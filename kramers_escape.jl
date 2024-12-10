using Pkg
Pkg.activate("ChaoticTipping")
using CSV 
using DataFrames
using DifferentialEquations
using LaTeXStrings
using Plots
using PrettyTables
using Random
include("parameters.jl")
include("utils.jl")

# Fill in details of precomputed lorenz data if needed
###############################################
num_precomputed = 0 #200
t_span_precomputed = (0.0, 1000000.0)
saveat_precomputed = 0.1
###############################################

# Select if to use SDE model (true) or chaotic ODE (false) 
noise_limit = true

# If ODE case: Select if to use precomputed lorenz trajectories
precomputed = false
# parameter for ODE: limit epsilon to 0 it converges to noise process # double well: ~ 0.1 AMOC ~1
epsilon = 1

# Set true for double well, false for AMOC 3box
double_well = false

# initial_condition (on/off) and param for 3 five_box, no hosing/T-forcing
paramset = 2
start_on_state = true
param = create_param(paramset)

# Number of realizations for each noise strength
num_ic = 100
# parameters for Kramers Law (noise/chaos)
# dw: [0.1, 0.2, 0.25, 0.3, 0.4, 0.5]/[0.25, 0.3, 0.4, 0.5] t:[90000, 5000, 1500, 500, 300, 200]
# 3box: 2 on/off: [1.2, 1.5, 2.0, 3.0, 4.0] 3 on: delta = [4.0, 5.0, 6.0, 7.5, 9.0] off [10, 12, 15, 19, 23] t:[100000, 55000, 30000, 20000, 15000].*1.5
delta = double_well ? [0.1, 0.2, 0.25, 0.3, 0.4, 0.5] : [1.2, 1.5, 2.0, 3.0, 4.0]
integration_len = double_well ? [90000, 5000, 1500, 500, 300, 200] : [100000, 55000, 30000, 20000, 15000].*1.5
# noise strengths
del_sc = sqrt.(2*delta) 
tipping_border = 1 .-delta
noise_len = length(delta)
# Estimation of standard deviation of limit of first lorenz coordinate noise
sigma = sqrt(60)
# δ is set seperatly
l_param = lorenz_parameters(28, 10, 8/3, epsilon, 0, sigma)
# initial_condition for double well and lorenz63 system
u0 = [-1.0]
model = precomputed ? dw_lorenz_precomputed! : dw_lorenz!
num_lorenz_trajs = num_ic

name = "dw"
sde_prob = SDEProblem(dw_drift!, one_d_diffusion!, u0, (0.0, 1000.0))
# change data for 3-box if necessary
if (double_well) 
    if !noise_limit && !precomputed ic_lorenz = sample_lorenz_attractor(num_ic*noise_len)[1] end
else
    num_lorenz_trajs *= 2
    if !noise_limit && !precomputed ic_lorenz = sample_lorenz_attractor(2*num_ic*noise_len)[1] end
    u0 = start_on_state ? param.on_state : param.off_state
    sde_prob = SDEProblem(three_box!, two_d_diffusion!, u0, (0.0, 1000.0), noise_rate_prototype = zeros(2, 2)) 
    model = precomputed ? three_box_lorenz_precomputed! : three_box_lorenz!
    tipping_border =  0.8 *(start_on_state ? param.q_off : param.q_on) .+ 0.0*delta
    name = "3b"
end
ode_prob = ODEProblem(model, u0, (0.0, 1000.0))

current_dir = @__DIR__
if !noise_limit && precomputed
    if (num_lorenz_trajs * noise_len > num_precomputed) && precomputed
        error("There are not enough precomputed Lorenz trajectories.")
    end

    # Choose precomputed trajectories
    rand_trajs = randperm(num_precomputed)[1:num_lorenz_trajs * noise_len]
    max_index = ceil(Int,(1/saveat_precomputed)*maximum(integration_len)/epsilon^2)+1

    # Load precomputed trajectories
    lorenz_trajs = Vector{}()
    for j in rand_trajs
        traj = CSV.read(joinpath(current_dir, "lorenz/trajectory_$j.csv"), DataFrame, limit=max_index, select=[2])[:,  1]
        push!(lorenz_trajs, traj)
    end

    for l in integration_len
        if (l/epsilon^2) >= t_span_precomputed[2] error("Precomputation to short for integration length ", l) end
    end

    println("trajectories loaded")
end
#H = -0.04
#Tf = 0


# switch noise strength, integration_length and lorenz trajectory/initial_condition
function prob_func_p(prob, i, repeat)
    # switch noise after every num_ic problem and lorenz traj for every realization
    noise_index = div(i-1, num_ic) + 1
    new_p = double_well ? (del_sc[noise_index], l_param, saveat_precomputed, lorenz_trajs[i]) : (del_sc[noise_index], l_param, saveat_precomputed, lorenz_trajs[i], lorenz_trajs[noise_len * num_ic+i], param)
    new_tspan = (0.0, integration_len[noise_index])
    remake(prob, p=new_p, tspan=new_tspan)
end
function prob_func_n(prob, i, repeat)
    noise_index = div(i-1, num_ic) + 1
    new_p = (del_sc[noise_index], param) # new_p = (del_sc[noise_index], param, forcing_parameters(true, H, H, H, 0, 1, 1, 1), forcing_parameters(true, Tf, Tf, Tf, 0, 1, 1, 1))
    new_tspan = (0.0, integration_len[noise_index]) 
    remake(prob, p=new_p, tspan=new_tspan)
end

function prob_func(prob, i, repeat)
    noise_index = div(i-1, num_ic) + 1
    new_p = (del_sc[noise_index], l_param, param)
    new_tspan = (0.0, integration_len[noise_index])
    new_u0 = double_well ? vcat(u0, ic_lorenz[:, i]) : vcat(u0, ic_lorenz[:, i], ic_lorenz[:, noise_len * num_ic+i])
    remake(prob, p=new_p, tspan=new_tspan, u0 = new_u0)
end

if noise_limit
    ensemble_prob = EnsembleProblem(sde_prob, prob_func = prob_func_n)
    ensemble_solution = solve(ensemble_prob, EM(), EnsembleThreads(), trajectories=noise_len * num_ic, dt = double_well ? 0.1 : 1)
else
    ensemble_prob = EnsembleProblem(ode_prob, prob_func = precomputed ? prob_func_p : prob_func)
    ensemble_solution = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=noise_len * num_ic, maxiters = 1e7)
end

# detect, which trajectories tipped at which time
tipped = Array{Any}(undef, noise_len, num_ic)
e_sol = Array{Any}(undef, noise_len, num_ic)
for i in 1:noise_len
    for j in 1:num_ic
        e_sol[i, j] = ensemble_solution[num_ic*(i-1)+j]
        # detect tipping:
            # double well: trajectory is above 1-δ
            # 3-box: check b = (q(t) <= tipping boarder). If started in 
                # on_state -> b == true -> tipped
                # off_state -> b == true -> not yet tipped
        tipped[i, j] = double_well ? e_sol[i, j].t[coalesce(findfirst(e_sol[i, j][1, :] .>= tipping_border[i]), 1)] : 
            e_sol[i, j].t[coalesce(findfirst(start_on_state .== (q0(e_sol[i, j][1, :], param) .<= tipping_border[i])), 1)] 
    end
end

for i in 1:noise_len
    if (sum(tipped .!= 0, dims=2)[i]) == 0 error("No trajectories tipped for δ = " * string(delta[i])) end
    if (sum(tipped .!= 0, dims=2)[i]) < num_ic @warn("Not all trajectories tipped for δ = " * string(delta[i]), ". Average result may be too low") end
end

# calculate mean
estimated_tipping = (sum(tipped, dims=2)./sum(tipped .!= 0, dims=2))[:, 1] 

# Assuming, drift term is a double well V = (x-1)^2(x+1)^2 with communication height 1 for plotting
pot_h = 1 
y_0 = log(2*pi/(sqrt(32)))

if(!double_well)
    if paramset == 2 && start_on_state 
        pot_h = 4.3 # 1.5 6.5
        y_0 = 6.2
    elseif paramset == 2 && !start_on_state 
        pot_h = 3.7 # 14 1.1
        y_0 = 6.9
    elseif paramset == 3 && start_on_state 
        pot_h = 13.3
        y_0 = 6.5
    elseif paramset == 3 && !start_on_state 
        pot_h = 40
        y_0 = 6.4
    else
        pot_h = 1
        y_0 = 0
    end
end
x_ax = (0, 5/(4*delta[1]))

# Plot 1/δ vs log tipping times
title = double_well ? "log transition times " * (noise_limit ? "" : (", " * L"\epsilon = " * string(epsilon))) : 
    "log transition times " * param.name  *  (noise_limit ? ", " : (", " * L"\epsilon = " * string(epsilon) * ",\n")) * "noise x" * string(param.noise_scaling)
p_detail = scatter(1 ./delta, log.(estimated_tipping), xlabel=L"\frac{1}{\delta}", ylabel="ln first transition time", size=(600, 600), 
    label="log avg τ", title=title, xlims = x_ax) # aspect_ratio=1,

plot!(p_detail, 1 ./ delta, pot_h  ./ delta .+ y_0, label="asymptotics", lw=2, linecolor=:red, linestyle=:dash)
if !noise_limit plot!(p_detail, 1 ./ delta, pot_h  ./ delta .+ y_0.-2, label="coef_error", lw=2, linecolor=:green, linestyle=:dash) end
if !noise_limit plot!(p_detail, 1 ./ delta, pot_h  ./ delta .+ y_0.+1, lw=2, label = missing, linecolor=:green, linestyle=:dash) end
for i in 1:noise_len
    vals = tipped[i, :][tipped[i, :] .!=0]
    scatter!(p_detail, 0. *vals .+ (1 /delta[i]), log.(vals), markershape=:cross, markersize=1, markercolor=:black, label=missing)
end
#p_detail = plot(1 ./epsilon, log.(estimated_tipping))
on_str = start_on_state ? "_on" : "_off"
appendix = noise_limit ? "_n" : (precomputed ? "_p" : "")
file_name = double_well ? "Kramers/kramers_" * name * (noise_limit ? "" : "_l_e" * string(epsilon)) * "_n" * string(num_ic) * appendix * ".png" : 
    "Kramers/kramers_" * name * "_" * string(paramset) * on_str * (noise_limit ? "" : "_l_e" * string(epsilon)) * "_n" * string(num_ic) *  appendix * ".png"
savefig(p_detail, joinpath(current_dir, file_name))
display(p_detail)

# Print results
header = ["δ", "sqrt{2δ}", "tipping boarder", "# tipped", "amount tipped", "average tipping time δ"]
data = hcat(
    delta,
    del_sc,
    tipping_border,
    count(tipped .!= 0, dims=2)[:, 1],
    string.(100 .*sum(tipped .!= 0, dims=2)[:, 1] ./ num_ic) .* "%",
    estimated_tipping
)
if !noise_limit println("epsilon = ", epsilon) end
pretty_table(data, header=header)

# Test plot for double well
#=
# plot a trajectory, to see if it has tipped
p_test = plot(label=missing,)
for j in 1:1
     plot!(p_test, e_sol[2, j], label = "1")
end
for j in 1:1
     plot!(p_test, e_sol[1, j], label = "2")
end
# x-axis is taken from the time span of last plot
display(p_test)
=#


# Test plot for 3-box
#=
i = length(delta)
j = 1
a1 = e_sol[i, 1][1, :]
a2 = e_sol[i, 1][2, :]
b1 = e_sol[j, 2][1, :]
b2 = e_sol[j, 2][2, :]

t1 = range(0, integration_len[i], length(a1))
t2 = range(0, integration_len[j], length(b1))
# Plot trajectories
p1 = plot(t1, a1, ylabel = L"\tilde S_N", label = missing, legend=:topright, color =:black)
hline!([param.on_state[1]], line=:dash, color=:blue, label="on_state")
hline!([param.off_state[1]], line=:dash, color=:red, label="off_state")
p2 = plot(t2, b1, label = missing, color =:darkgreen)
hline!([param.on_state[1]], line=:dash, color=:blue, label = missing)
hline!([param.off_state[1]], line=:dash, color=:red, label = missing)

p3 = plot(t1, a2, ylabel = L"\tilde S_T", label = missing, legend=:topright, color =:black)
hline!([param.on_state[2]], color=:blue, label="on_state")
hline!([param.off_state[2]], line=:dash, color=:red, label="off_state")
p4 = plot(t2, b2, label = missing, color =:darkgreen)
hline!([param.on_state[2]], color=:blue, label = missing)
hline!([param.off_state[2]], line=:dash, color=:red, label = missing)

# Plot q
#q_values = [q(S_N, (param, T_param), 0) for S_N in a1]
q_values = [q0(a1[i], param) for i in 1:length(a1)]
q_values_2 = [q0(b1[i], param) for i in 1:length(b1)]
p5 = plot(t1, q_values.*1e-6, ylabel=L"q\ [Sv]", label = missing, legend=:topright, color =:black)
hline!([param.q_on]*1e-6, line=:dash, color=:blue, label="on_state")
hline!([param.q_off]*1e-6, line=:dash, color=:red, label="off_state")
p6 = plot(t2, q_values_2.*1e-6, label = missing, color =:darkgreen)
hline!([param.q_on]*1e-6, line=:dash, color=:blue, label = missing)
hline!([param.q_off]*1e-6, line=:dash, color=:red, label = missing)

p_d = plot(p1, p2, p3, p4, p5, p6, layout=(:,1), size=(800, 700))
savefig(p_d, "Kramers_traj.png")
display(p_d) 
println("Green: δ = ", string(delta[j]), ", Black: δ = ", string(delta[i]))
=#