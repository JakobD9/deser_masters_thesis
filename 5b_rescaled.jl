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

# Define time span
initial_time = 0.0
final_time = 5000
num_points = 10000
initial_condition = [param.S_N, param.S_T, param.S_S, param.S_IP]
initial_condition = param.on_state_5

# temperature forcing parameters:
forcing_T = true
T_0 = 0
T_pert = 4
T_end = -3

# hosing parameters:
forcing_h = false
H_0 = 0
H_end = -0.3
H_pert = 0.4

# Generell forcing params
t_f = 0
t_rise = 1100
t_pert = 500
t_fall = 2000
t_fall_2 = 1800

h_param = forcing_parameters(forcing_h, H_0, H_end, H_pert, t_f, t_rise, t_pert, t_fall)
h_param_2 = forcing_parameters(forcing_h, H_0, H_end, H_pert, t_f, t_rise, t_pert, t_fall_2)

T_param = forcing_parameters(forcing_T, T_0, T_end, T_pert, t_f, t_rise, t_pert, t_fall)
T_param_2 = forcing_parameters(forcing_T, T_0, T_end, T_pert, t_f, t_rise, t_pert, t_fall_2)

# Set model up
t_span = (initial_time, final_time)

dt=(final_time - initial_time) / (num_points - 1)

# Calculate solutions
prob = ODEProblem(five_box, initial_condition, t_span, (param, h_param, T_param, missing))
sol = solve(prob, Tsit5(), saveat = dt)
x1 = [u[1] for u in sol.u]
x2 = [u[2] for u in sol.u]
x3 = [u[3] for u in sol.u]
x4 = [u[4] for u in sol.u]

prob2 = ODEProblem(five_box, initial_condition, t_span, (param, h_param_2, T_param_2, missing))
sol2 = solve(prob2, Tsit5(), saveat = dt)
y1 = [u[1] for u in sol2.u]
y2 = [u[2] for u in sol2.u]
y3 = [u[3] for u in sol2.u]
y4 = [u[4] for u in sol2.u]

#icn = [x1[num_points], x2[num_points], x3[num_points], x4[num_points]]
#println(round.(icn, digits=7))
#q_end = [q(x1[num_points], x3[num_points], (param, missing), 0)]
#println(round.(q_end, digits=6)[1])

# Plot results
t_range = range(initial_time, final_time, num_points)
# Plot trajectories
p11 = plot(t_range, x1, ylabel = L"\tilde S_N", label = missing, legend=:topright)
plot!(t_range, y1, label = missing)
hline!([param.on_state_5[1]], line=:dash, color=:blue, label="on_state")
hline!([param.off_state_5[1]], line=:dash, color=:red, label="off_state")
p12 = plot(t_range, x2, ylabel = L"\tilde S_T", label = missing, legend=:topright)
plot!(t_range, y2, label = missing)
hline!([param.on_state_5[2]], line=:dash, color=:blue, label="on_state")
hline!([param.off_state_5[2]], line=:dash, color=:red, label="off_state")
p13 = plot(t_range, x3, ylabel = L"\tilde S_S", label = missing, legend=:topright)
plot!(t_range, y3, label = missing)
hline!([param.on_state_5[3]], line=:dash, color=:blue, label="on_state")
hline!([param.off_state_5[3]], line=:dash, color=:red, label="off_state")
p14 = plot(t_range, x4, ylabel = L"\tilde S_{IP}", label = missing, legend=:topright)
plot!(t_range, y4, label = missing)
hline!([param.on_state_5[4]], line=:dash, color=:blue, label="on_state")
hline!([param.off_state_5[4]], line=:dash, color=:red, label="off_state")
# Plot q
q_values = [q(x1[i], x3[i], (param, T_param), t_range[i]) for i in 1:num_points]
q_values_2 = [q(y1[i], y3[i], (param, T_param_2), t_range[i]) for i in 1:num_points]
p3 = plot(t_range, q_values.*1e-6, ylabel=L"q \ [Sv]", label = missing, legend=:topright)
plot!(t_range, q_values_2.*1e-6, label = missing)
hline!([param.q_on_5]*1e-6, line=:dash, color=:blue, label="on_state")
hline!([param.q_off_5]*1e-6, line=:dash, color=:red, label="off_state")
# Plot hosing
h_values = forcing_h ? [forcing_pwl(t, h_param) for t in t_range] : 0 .*t_range
h_values_2 = forcing_h ? [forcing_pwl(t, h_param_2) for t in t_range] : 0 .*t_range
p4 = plot(t_range, h_values, xlabel="Time (t) [Years]", ylabel=L"H(t)", label = missing, color=:blue)
plot!(t_range, h_values_2, label = missing, color=:red)
# Plot Temperature forcing
T_values = forcing_T ? [forcing_pwl(t, T_param) for t in t_range] : 0 .*t_range
T_values_2 = forcing_T ? [forcing_pwl(t, T_param_2) for t in t_range] : 0 .*t_range
p5 = plot(t_range, T_values, xlabel="Time (t) [Years]", ylabel=L"T^{f}", label = missing)
plot!(t_range, T_values_2, label = missing)
# Plot TN = \mu q + T_0
TN_values = q_values.*param.mu .+param.T_0 .+ T_values
TN_values_2 = q_values_2.*param.mu .+param.T_0 .+T_values_2
p6 = plot(t_range, TN_values, ylabel=L"\mu q + T_0 + T^f", label = missing)
plot!(t_range, TN_values_2, label = missing)
#hline!([param.q_on], line=:dash, color=:blue, label="on_state")

if forcing_h && forcing_T
    p_detail = plot(p11, p12, p13, p14, p3, p6, p5, p4, layout=(8,1), size=(800, 900))
elseif forcing_h && !forcing_T
    p_detail = plot(p11, p12, p13, p14, p3, p4, layout=(6,1), size=(800, 700))
elseif !forcing_h && forcing_T
    p_detail = plot(p11, p12, p13, p14, p3, p6, p5, layout=(7,1), size=(800, 800))
elseif !forcing_h && !forcing_T
    p_detail = plot(p11, p12, p13, p14, p3, p6, layout=(6,1), size=(800, 700))
end

current_dir = @__DIR__
savefig(p_detail, joinpath(current_dir, "trajs/" * string(paramset) * "_forcing_detail_5.png"))

display(p_detail)

# Plot 2d phase space
p_phase = plot(x1, x2, xlabel="S_N", ylabel="S_T", label = "traj_1")
plot!(y1, y2, label = "traj_2")
scatter!([param.on_state[1]], [param.on_state[2]], color=:blue, label="on_state")
scatter!([param.off_state[1]], [param.off_state[2]], color=:red, label="off_state")
savefig(p_phase, joinpath(current_dir, "trajs/" * string(paramset) * "_forcing_traj_5.png"))
display(p_phase)

SB_pw(S_N, S_T, S_S, S_IP, p) = (100 * (p.C - p.S_0*(p.V_B+p.V_N+p.V_T+p.V_IP+p.V_S)) .- (p.V_N .* S_N .+ p.V_T .* S_T .+ p.V_S .* S_S .+ p.V_IP .* S_IP) ) ./ p.V_B

x5 = SB_pw(x1, x2, x3, x4, param)
y5 = SB_pw(y1, y2, y3, y4, param)
S_Nx = x1./100 .+ param.S_0
S_Tx = x2./100 .+ param.S_0
S_Sx = x3./100 .+ param.S_0
S_IPx = x4./100 .+ param.S_0
S_Bx = x5./100 .+ param.S_0
S_Ny = y1./100 .+ param.S_0
S_Ty = y2./100 .+ param.S_0
S_Sy = y3./100 .+ param.S_0
S_IPy = y4./100 .+ param.S_0
S_By = y5./100 .+ param.S_0

p_scale = plot(t_range, S_Nx, ylabel = "Salinity", label = L"S_N", color = :red, legend=:topright)
plot!(t_range, S_Tx, label = L"S_T", color = :blue)
plot!(t_range, S_Sx, label = L"S_S", color = :turquoise)
plot!(t_range, S_IPx, label = L"S_{IP}", color = :orange)
plot!(t_range, S_Bx, label = L"S_B", color = :black)

savefig(p_scale, joinpath(current_dir, "trajs/" * string(paramset) * "_forcing_detail_backscale_5.png"))
display(p_scale)