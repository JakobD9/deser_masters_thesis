using Pkg
Pkg.activate("ChaoticTipping")
using DelimitedFiles
using DifferentialEquations
using LaTeXStrings
using LinearAlgebra
using Plots
include("parameters.jl")
include("utils.jl")

# Define parameters
paramset = 2
param = create_param(paramset)

# choose model function_ Either three_box or five_box
model = three_box

# Define time span
initial_time = 0.0
final_time = 5000
num_points = 1000

# Define phase space window
x_len = 200
y_len = 200
x_ax = [-0.4, 0.4] #paramset 3 [-0.5, 0.5] paramset 2 [-0.4, 0.4]  [-0.02, 0.01]
y_ax = [-0.1, 0.3] #paramset 3 [0, 0.75] paramset 2 [-0.1, 0.3] [0.08, 0.12]

#forcing parameters:
forcing_h = true
H_pert = -0.05 #0.2115

forcing_T = false
TN_pert = 2.9

# set run up
h_param = forcing_parameters(forcing_h, 0, 0, H_pert, 0, 0, final_time, 0)
T_param = forcing_parameters(forcing_T, 0, 0, TN_pert, 0, 0, final_time, 0)
l_param = lorenz_parameters(28, 10, 8/3, 0, 0, sqrt(60))

x_range = range(x_ax[1], stop=x_ax[2], length=x_len)
y_range = range(y_ax[1], stop=y_ax[2], length=y_len)

initial_conditions = [[x, y] for x in x_range, y in y_range]
if occursin("five", string(model))
    initial_conditions = [vcat(ic, [param.S_S, param.S_IP]) for ic in initial_conditions]
    println("five_detected")
end

t_span = (initial_time, final_time)
dt=(final_time - initial_time) / (num_points - 1)

function solve_ode(initial_condition)
    prob = ODEProblem(model, initial_condition, t_span, (param, h_param, T_param, l_param))
    sol = solve(prob, Tsit5(), dt=dt)
    return sol
end

# Solve the ODE for each initial condition
solutions = [solve_ode(ic) for ic in initial_conditions]
S_N = [solutions[i, j][1, :] for i in 1:x_len, j in 1:y_len]
S_T = [solutions[i, j][2, :] for i in 1:x_len, j in 1:y_len]
last_values = [solutions[i, j][end] for i in 1:x_len, j in 1:y_len]
test_ind = Int(num_points/200)
test_end_values = [solutions[i, j][end - test_ind] for i in 1:x_len, j in 1:y_len]
diff = norm.(last_values.-test_end_values)

println("Not convergent initial conditions: ", 100*sum(diff .>= 1e-3)/(x_len*y_len), "%")

last_S_N = [last_values[i, j][1, :][1] for i in 1:x_len, j in 1:y_len]
if occursin("five", string(model))
    last_S_S = [last_values[i, j][3, :][1] for i in 1:x_len, j in 1:y_len]
    q_end = q(last_S_N, last_S_S, (param, T_param), 0)
else
    q_end = q(last_S_N, (param, T_param), 0)
end

# decide if trajectory ended in on/off state
on_basin = q_end .>= 0

on_vectors = last_values[on_basin]
off_vectors = last_values[.!on_basin]
#println("fraction on state: ", length(on_vectors)/(x_len*y_len))
# max/min are w.r.t first component as of lexicographic ordering of vectors
if (length(on_vectors) > 0)
    println("variance on ", norm(maximum(on_vectors) - minimum(on_vectors)))
    avg_on = sum(on_vectors)/length(on_vectors)
    #println(avg_on)
end

if (length(off_vectors) > 0)
    println("variance off ", norm(maximum(off_vectors) - minimum(off_vectors)))
    avg_off = sum(off_vectors)/length(off_vectors)
    #println(avg_off)
end

# Define colors for true and false values
color_true = :blue
color_false = :red

# Convert ranges to matrices
X = repeat(x_range, 1, y_len)
Y = repeat(y_range', x_len, 1)

title = "Basins of attraction for " * (forcing_T ? L"T^{f} = " * string(TN_pert) : "") * 
    ((forcing_h && forcing_T) ? ", " : "" ) * (forcing_h ? L" H_{pert} = " * string(H_pert) : "")
p1 = scatter(title=title)

# Plot points based on the boolean array
scatter!(X[on_basin], Y[on_basin], color=color_true, markersize=3, label=missing)
scatter!(X[.!on_basin], Y[.!on_basin], color=color_false, markersize=3, label=missing)
if length(on_vectors) > 0  scatter!([avg_on[1]], [avg_on[2]], color = :green, label="on") end
if length(off_vectors) > 0  scatter!([avg_off[1]], [avg_off[2]], color = :yellow , label="off") end


xlabel!(L"\tilde S_N")
ylabel!(L"\tilde S_T")
current_dir = @__DIR__
savefig(p1, joinpath(current_dir, "basins/basin_" * (forcing_T ? "T" * string(TN_pert)*"_" : "") * 
     (forcing_h ? "h" * string(H_pert)*"_" : "") * "paramset_" * string(paramset) * ".png"))
display(p1)