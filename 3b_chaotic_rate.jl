using Pkg
Pkg.activate("ChaoticTipping")
using ColorTypes
using DelimitedFiles
using DifferentialEquations
using LaTeXStrings
using Plots
include("parameters.jl")
include("utils.jl")

# number initial conditions
num_ic = 15 
# Define 3-box parameters
paramset = 2
param = create_param(paramset)

# Define time span
initial_time = 0.0
final_time = 3000
num_points = 10000

initial_condition = param.on_state

# lorenz forcing parameters: 
forcing_l = true
coupled_l = true # 2d lorenz forcing
epsilon = 3
delta = 1.5 #0.5 # ->  sqrt(2δ) = 1
sigma = sqrt(60)

# ignores all settings below and plots the chaotic tipping window for hosing
c_tipping_window = false

# temperature forcing parameters:
forcing_T = false
T_pert = 1.4 # 2 for paramset = 2, 5 for paramset = 3
T_end = 1

# hosing parameters
forcing_h = true
H_pert = 0.35
H_end = -0.04

# general forcing parameters
t_mid = 200 # Determines how wide the forcing should vary
t_f = 0 #100
up = false

r_min = abs(T_pert-0)/(2*(t_mid-t_f))
t_var_max = 2*(t_mid-t_f)

T_param_all = Vector{}()
h_param_all = Vector{}()
for k in 1:num_ic
    t_var = t_var_max - (k-1)*(t_var_max-1)/num_ic
    T_p = missing
    if up
        fall = 500
        T_p = forcing_parameters(forcing_T, 0, T_end, T_pert, t_mid-0.5*t_var, t_var, t_mid-0.5*t_var, fall)
        H_p = forcing_parameters(forcing_h, 0, H_end, H_pert, t_mid-0.5*t_var, t_var, t_mid-0.5*t_var, fall)
    else
        rise = 20
        pert = 80
        T_p = forcing_parameters(forcing_T, 0, T_end, T_pert, t_f, rise, pert, t_var)
        H_p = forcing_parameters(forcing_h, 0, H_end, H_pert, t_f, rise, pert, t_var)
    end
    if c_tipping_window
        T_p = forcing_parameters(false, 0, 0, 0, 0, 0, 0, 0)
        H_p = forcing_parameters(true, 0, 0, 0.35, 0, final_time, 2000, 0)
    end
    push!(T_param_all, T_p)
    push!(h_param_all, H_p)
end

l_param = lorenz_parameters(28, 10, 8/3, epsilon, forcing_l ? delta : 0, sigma)
# Define model and initial condition dependent on settings on chaotic forcing
model = three_box
dim = 2
if forcing_l ic_lorenz, p = coupled_l ? sample_lorenz_attractor(2*num_ic) : sample_lorenz_attractor(num_ic) end
if (forcing_l) 
    model = three_box_chaotic
    ic = [vcat(initial_condition, ic_lorenz[:, k]) for k in 1:num_ic]
    initial_condition = ic
    dim += 3
    if (coupled_l)
        model = three_box_coupled_chaos
        ic = [vcat(initial_condition[k], ic_lorenz[:, num_ic + k]) for k in 1:num_ic]
        initial_condition = ic
        dim += 3
    end
else
    ic = [initial_condition for k in 1:num_ic]
    initial_condition = ic
end

# Calculate trajectories
t_span = (initial_time, final_time)
dt=(final_time - initial_time) / (num_points - 1)

function prob_func(prob, i, repeat)
    new_u0 = initial_condition[i]
    new_p = (param, h_param_all[i], T_param_all[i], l_param)
    remake(prob, p=new_p, u0 = new_u0)
end

prob = ODEProblem(model, initial_condition[1], t_span)
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
ensemble_solution = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=num_ic, maxiters = 1e7, saveat=dt)

# Plot results
t_range = range(initial_time, final_time, num_points)
colors = color_gradient(colorant"yellow", colorant"purple", num_ic)

# Plot trajectories
p1 = plot(ylabel = L"\tilde S_N")
p2 = plot(ylabel = L"\tilde S_T")
p3 = plot(ylabel=L"q\ [Sv]")
p4 = plot(ylabel=L"H_{pert}")
p5 = plot(ylabel=L"T^{f}")
p6 = plot(ylabel=L"\mu q + T_0 + T^f")
if(c_tipping_window) 
    forcing_h = true
    forcing_T = false
end
for k in 1:num_ic
    trajk = [[u[i] for u in ensemble_solution[k].u] for i in 1:2]
    u1, u2 = trajk
    plot!(p1, t_range, u1, label = missing, linecolor=colors[k])
    plot!(p2, t_range, u2, label = missing, linecolor=colors[k])
    q_values = [q(u1[i], (param, T_param_all[k]), t_range[i]) for i in 1:num_points]
    plot!(p3, t_range, q_values.*1e-6, label = missing, linecolor=colors[k])
    h_values = forcing_h ? [forcing_pwl(t, h_param_all[k]) for t in t_range] : 0 .*t_range
    plot!(p4, t_range, h_values, label = missing, linecolor=colors[k], xlabel="Time (t) [Years]")
    T_values = forcing_T ? [forcing_pwl(t, T_param_all[k]) for t in t_range] : 0 .*t_range
    plot!(p5, t_range, T_values, label = missing, linecolor=colors[k])
    TN_values = q_values.*param.mu .+param.T_0 .+ T_values
    plot!(p6, t_range, TN_values, label = missing, linecolor=colors[k])
end
if(c_tipping_window)
    h = [forcing_pwl(t, h_param_all[1]) for t in t_range]
    bif_index = argmin(abs.(h .- (paramset == 2 ?  0.2134 :  0.3895)))
    bif_time = initial_time + bif_index*dt
    if paramset in (2, 3) vline!(p4, [bif_time], label="Hopf bifurcation", line=:solid, color=:green) end
end

hline!(p1, [param.on_state[1]], line=:dash, color=:blue, label = missing)#label="on_state")
hline!(p1, [param.off_state[1]], line=:dash, color=:red, label = missing)#label="off_state")
hline!(p2, [param.on_state[2]], line=:dash, color=:blue, label = missing)#label="on_state")
hline!(p2, [param.off_state[2]], line=:dash, color=:red, label = missing)#label="off_state")
hline!(p3, [param.q_on]*1e-6, line=:dash, color=:blue, label = missing)#label="on_state")
hline!(p3, [param.q_off]*1e-6, line=:dash, color=:red, label = missing)#label="off_state"

if forcing_h && forcing_T
    p_detail = plot(p1, p2, p3, p6, p5, p4, layout=(6,1), size=(800, 700))
elseif forcing_h && !forcing_T
    p_detail = plot(p1, p2, p3, p4, layout=(4,1), size=(800, 500))
elseif !forcing_h && forcing_T
    xlabel!(p5, "Time (t) [Years]")
    p_detail = plot(p1, p2, p3, p6, p5, layout=(5,1), size=(800, 600))
elseif !forcing_h && !forcing_T
    xlabel!(p3, "Time (t) [Years]")
    p_detail = plot(p1, p2, p3, layout=(3,1), size=(800, 400))
end
current_dir = @__DIR__
savefig(p_detail, joinpath(current_dir, "trajs/" * string(paramset) * "_chaotic_forcing_rate_detail" * (up ? "_up" : "") * ".png"))
display(p_detail)

# Plot 2d phase space
p_phase = plot(xlabel=L"\tilde S_N", ylabel=L"\tilde S_T")
for k in 1:num_ic
    trajk = [[u[i] for u in ensemble_solution[k].u] for i in 1:2]
    u1, u2 = trajk
    t_var = floor(Int64, (t_var_max - (k-1)*(t_var_max-1)/num_ic))
    plot!(u1, u2,  linecolor=colors[k], label = missing) #label = string(t_var))
end
scatter!([param.on_state[1]], [param.on_state[2]], color=:blue, label="on_state")
scatter!([param.off_state[1]], [param.off_state[2]], color=:red, label="off_state")
savefig(p_phase, joinpath(current_dir, "trajs/" * string(paramset) * "_chaotic_forcing_rate_traj" * (up ? "_up" : "") * ".png"))
#display(p_phase)