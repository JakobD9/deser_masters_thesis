using Pkg
Pkg.activate("ChaoticTipping")
using DelimitedFiles
using DifferentialEquations
using LaTeXStrings
using Plots
include("parameters.jl")
include("utils.jl")

# Define parameters
paramset = 2
param = create_param(paramset)
initial_condition = param.on_state


# Define time span
initial_time = 0.0
final_time = 1500
num_points = 10000

# lorenz forcing parameters:
forcing_l = true
coupled_l = true
epsilon = 1 # 0.2
delta = 1.6 #0.5 # ->  prefactor sqrt(2δ)
sigma = sqrt(60)

# temperature forcing parameters:
forcing_T = false
rate_TN = 0.0075
rate_TN_2 = 0.05
TN_pert = 3.05
use_hosing_rise_fall = true

# hosing parameters:
forcing_h = true
H_0 = 0
H_end = 0
H_pert = 0.15 #0.3706
t_rise = 50
t_pert = 100
t_fall = 500
t_fall_2 = 100

# Generell forcing params
t_f = 0

h_param = forcing_parameters(forcing_h, H_0, H_0, H_pert, t_f, t_rise, t_pert, t_fall)
h_param_2 = forcing_parameters(forcing_h, H_0, H_0, H_pert, t_f, t_rise, t_pert, t_fall_2)
# Temperature forcing upramp only
T_param = forcing_parameters(forcing_T, 0, 0, TN_pert, t_f, abs(TN_pert-0)/rate_TN, final_time, 0)
T_param_2 = forcing_parameters(forcing_T, 0, 0, TN_pert, t_f, abs(TN_pert-0)/rate_TN_2, final_time, 0)
# Overwrite for up/downramp forcing
if use_hosing_rise_fall
    T_param = forcing_parameters(forcing_T, H_0, H_0, TN_pert, t_f, t_rise, t_pert, t_fall)
    T_param_2 = forcing_parameters(forcing_T, H_0, H_0, TN_pert, t_f, t_rise, t_pert, t_fall_2)
end
l_param = lorenz_parameters(28, 10, 8/3, epsilon, forcing_l ? delta : 0, sigma)

# Set model up
dim = 2
if forcing_l ic_lorenz, p = coupled_l ? sample_lorenz_attractor(4) : sample_lorenz_attractor(2) end
model = three_box
if (forcing_l) 
    model = three_box_chaotic
    ic = [vcat(initial_condition, ic_lorenz[:, k]) for k in 1:2]
    initial_condition = ic
    dim += 3
    if (coupled_l)
        model = three_box_coupled_chaos
        ic = [vcat(initial_condition[k], ic_lorenz[:, 2 + k]) for k in 1:2]
        initial_condition = ic
        dim += 3
    end
else
    ic = [initial_condition for k in 1:2]
    initial_condition = ic
end

t_span = (initial_time, final_time)

dt=(final_time - initial_time) / (num_points - 1)

# Calculate solutions
prob = ODEProblem(model, initial_condition[1], t_span, (param, h_param, T_param, l_param))
sol = solve(prob, Tsit5(), saveat = dt)

prob2 = ODEProblem(model, initial_condition[2], t_span, (param, h_param_2, T_param_2, l_param))
sol2 = solve(prob2, Tsit5(), saveat = dt)

traja = [[u[i] for u in sol.u] for i in 1:dim]
trajb = [[u[i] for u in sol2.u] for i in 1:dim]


if dim == 2
    a1, a2 = traja 
    b1, b2 = trajb 
    a3 = a4 = a5 = a6 = a6 = a8 = zeros(num_points)
    b3 = b4 = b5 = b6 = b7 = b8 = zeros(num_points)
elseif dim == 5
    a1, a2, a3, a4, a5 = traja 
    b1, b2, b3, b4, b5 = trajb 
    a6 = a7 = a8 = zeros(num_points)
    b6 = b7 = b8 = zeros(num_points)
else
    a1, a2, a3, a4, a5, a6, a7, a8 = traja
    b1, b2, b3, b4, b5, b6, b7, b8 = trajb 
end

#icn = [a1[num_points], a2[num_points]]
#q_end = [q(a1[num_points], param)]

# Plot results
t_range = range(initial_time, final_time, num_points)

# Plot trajectories
p1 = plot(t_range, a1, ylabel = L"\tilde S_N", label = missing, legend=:topright, color=:black)
plot!(t_range, b1, label = missing, color=:darkgreen)
hline!([param.on_state[1]], line=:dash, color=:blue, label="on_state")
hline!([param.off_state[1]], line=:dash, color=:red, label="off_state")
p2 = plot(t_range, a2, ylabel = L"\tilde S_T", label = missing, legend=:topright, color=:black)
plot!(t_range, b2, label = missing, color=:darkgreen)
hline!([param.on_state[2]], color=:blue, label="on_state")
hline!([param.off_state[2]], line=:dash, color=:red, label="off_state")

# Plot q
q_values = [q(a1[i], (param, T_param), t_range[i]) for i in 1:num_points]
q_values_2 = [q(b1[i], (param, T_param_2), t_range[i]) for i in 1:num_points]
p3 = plot(t_range, q_values.*1e-6, ylabel=L"q\ [Sv]", label = missing, legend=:topright, color=:black)
plot!(t_range, q_values_2.*1e-6, label = missing, color=:darkgreen)
hline!([param.q_on]*1e-6, line=:dash, color=:blue, label="on_state")
hline!([param.q_off]*1e-6, line=:dash, color=:red, label="off_state")

# Plot hosing
h_values = forcing_h ? [forcing_pwl(t, h_param) for t in t_range] : 0 .*t_range
h_values_2 = forcing_h ? [forcing_pwl(t, h_param_2) for t in t_range] : 0 .*t_range
p4 = plot(t_range, h_values, ylabel=L"H(t)", label = missing, color=:black)
plot!(t_range, h_values_2, label = missing, color=:darkgreen)

# Plot Temperature forcing
T_values = forcing_T ? [forcing_pwl(t, T_param) for t in t_range] : 0 .*t_range
T_values_2 = forcing_T ? [forcing_pwl(t, T_param_2) for t in t_range] : 0 .*t_range
p5 = plot(t_range, T_values, ylabel=L"T^{f}", label = missing, color=:black)
plot!(t_range, T_values_2, label = missing, color=:darkgreen)

# Plot TN = \mu q + T_0
TN_values = q_values.*param.mu .+param.T_0 .+ T_values
TN_values_2 = q_values_2.*param.mu .+param.T_0 .+T_values_2
p6 = plot(t_range, TN_values, ylabel=L"\mu q + T_0 + T^f", label = missing, color=:black)
plot!(t_range, TN_values_2, label = missing, color=:darkgreen)
hline!([param.q_on]*param.mu .+param.T_0, line=:dash, color=:blue, label="on_state")
hline!([param.q_off]*param.mu .+param.T_0, line=:dash, color=:red, label="off_state")

# Plot lorenz forcing
p7 = plot(t_range, sqrt(2*l_param.delta)* 1/l_param.epsilon * 1/l_param.sigma * ((coupled_l ? param.B11 : 1)*a3 + 0*a6), label = missing, ylabel = L"f_0(y)_1", color=:black)
plot!(t_range, sqrt(2*l_param.delta)* 1/l_param.epsilon * 1/l_param.sigma * ((coupled_l ? param.B11 : 1)*b3 + 0*b6), label = missing, color=:darkgreen)

p8 = plot(t_range, sqrt(2*l_param.delta)* 1/l_param.epsilon * 1/l_param.sigma * ((coupled_l ? param.B21 : 0)*a3 + (coupled_l ? param.B22 : 0)*a6), label = missing, xlabel="Time (t) [Years]", ylabel = L"f_0(y)_2", color=:black)
plot!(t_range, sqrt(2*l_param.delta)* 1/l_param.epsilon * 1/l_param.sigma * ((coupled_l ? param.B21 : 0)*b3 + (coupled_l ? param.B22 : 0)*b6), label = missing, color=:darkgreen)

if forcing_h && forcing_T
    p_detail = plot(p1, p2, p3, p6, p5, p4, p7, p8, layout=(:,1), size=(800, 900))
elseif forcing_h && !forcing_T
    p_detail = plot(p1, p2, p3, p4, p7, p8, layout=(:,1), size=(800, 700))
elseif !forcing_h && forcing_T
    p_detail = plot(p1, p2, p3, p6, p5, p7, p8, layout=(:,1), size=(800, 800))
elseif !forcing_h && !forcing_T
    p_detail = plot(p1, p2, p3, p7, p8, layout=(:,1), size=(800, 600)) 
end
# Plots without chaotic forcing
#=
if forcing_h && forcing_T
    p_detail = plot(p1, p2, p3, p6, p5, p4,  layout=(:,1), size=(800, 700))
elseif forcing_h && !forcing_T
    p_detail = plot(p1, p2, p3, p4,  layout=(:,1), size=(800, 500))
elseif !forcing_h && forcing_T
    p_detail = plot(p1, p2, p3, p6, p5, layout=(:,1), size=(800, 600))
elseif !forcing_h && !forcing_T
    p_detail = plot(p1, p2, p3,  layout=(:,1), size=(800, 400))
end
=#
current_dir = @__DIR__
savefig(p_detail, joinpath(current_dir, "trajs/" * string(paramset) * (forcing_l ? "_chaotic" : "") * "_forcing_detail.png"))

display(p_detail)

# Plot 2d phase space
p_phase = plot(a1, a2, xlabel=L"\tilde S_N", ylabel=L"\tilde S_T", label = "traj_1", color=:black)
plot!(b1, b2, label = "traj_2", color=:darkgreen)
scatter!([param.on_state[1]], [param.on_state[2]], color=:blue, label="on_state")
scatter!([param.off_state[1]], [param.off_state[2]], color=:red, label="off_state")
savefig(p_phase, joinpath(current_dir, "trajs/" * string(paramset) * (forcing_l ? "_chaotic" : "") * "_forcing_traj.png"))
#display(p_phase)

# Estimate covariance matrices from trajectories
#=
println(cov(hcat(a1, a2)))
println(cov(hcat(b1, b2)))
l1 = sqrt(2*l_param.delta) * 1/l_param.epsilon * 1/l_param.sigma .* (param.B11.*a3 .+ 0 .*a6)
l2 = sqrt(2*l_param.delta)* 1/l_param.epsilon * 1/l_param.sigma .* (param.B21.*a3 .+ param.B22.*a6)
println(cov(hcat(l1, l2)))
=#

# Plot stochastic system
#=
u0 = initial_condition[1][1:2]
p=(sqrt(2*delta), param, h_param, T_param)
p2=(sqrt(2*delta), param, h_param_2, T_param_2)
sde_prob = SDEProblem(three_box_f!, two_d_diffusion!, u0, t_span, p, noise_rate_prototype = zeros(2, 2)) 
ssol = solve(sde_prob, EM(), dt=dt)

sde_prob2 = SDEProblem(three_box_f!, two_d_diffusion!, u0, t_span, p2, noise_rate_prototype = zeros(2, 2)) 
ssol2 = solve(sde_prob, EM(), dt=dt)

c1, c2 = [[u[i] for u in ssol.u] for i in 1:2]
d1, d2 = [[u[i] for u in ssol2.u] for i in 1:2]

# Plot trajectories
ps1 = plot(t_range, c1, ylabel = L"\tilde S_N", label = missing, legend=:topright, color=:black)
plot!(t_range, d1, label = missing, color=:darkgreen)
hline!([param.on_state[1]], line=:dash, color=:blue, label="on_state")
hline!([param.off_state[1]], line=:dash, color=:red, label="off_state")
ps2 = plot(t_range, c2, ylabel = L"\tilde S_T", label = missing, legend=:topright, color=:black)
plot!(t_range, d2, label = missing, color=:darkgreen)
hline!([param.on_state[2]], color=:blue, label="on_state")
hline!([param.off_state[2]], line=:dash, color=:red, label="off_state")

# Plot q
q_values = [q(c1[i], (param, T_param), t_range[i]) for i in 1:num_points]
q_values_2 = [q(d1[i], (param, T_param_2), t_range[i]) for i in 1:num_points]
ps3 = plot(t_range, q_values.*1e-6, ylabel=L"q\ [Sv]", label = missing, legend=:topright, color=:black)
plot!(t_range, q_values_2.*1e-6, label = missing, color=:darkgreen)
hline!([param.q_on]*1e-6, line=:dash, color=:blue, label="on_state")
hline!([param.q_off]*1e-6, line=:dash, color=:red, label="off_state")

if forcing_h && forcing_T
    ps_detail = plot(ps1, ps2, ps3, p5, p4,  layout=(:,1), size=(800, 600))
elseif forcing_h && !forcing_T
    ps_detail = plot(ps1, ps2, ps3, p4,  layout=(:,1), size=(800, 500))
elseif !forcing_h && forcing_T
    ps_detail = plot(ps1, ps2, ps3, p5, layout=(:,1), size=(800, 500))
elseif !forcing_h && !forcing_T
    ps_detail = plot(ps1, ps2, ps3,  layout=(:,1), size=(800, 400))
end

savefig(ps_detail, joinpath(current_dir, "trajs/" * string(paramset) * "_stochastic_forcing_detail.png"))
=#