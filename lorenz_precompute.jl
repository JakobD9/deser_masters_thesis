using Pkg
Pkg.activate("ChaoticTipping")
using CSV
using DataFrames
using DelimitedFiles
using DifferentialEquations
include("parameters.jl")
include("utils.jl")

# Saves time and all 3 components of lorenz trajectory. Attention, with settings below ~ 60GB data
# How many trajectories to precompute
num_ic = 5

# Time span and saving intervals (saves every t = k*saveat from k = t_span[1] to k = t_span[2]/saveat)
t_span = (0.0, 1000000.0)
saveat = 0.1

# sample initial conditions
x, p = sample_lorenz_attractor(num_ic)
initial_conditions = [Vector(x[:, i]) for i in 1:size(x, 2)]


l_param = lorenz_parameters(28, 10, 8/3, 0, 0, 0)

for i in 1:num_ic
    prob = ODEProblem(lorenz!, initial_conditions[i], t_span, l_param)
    sol = solve(prob, Tsit5(), saveat=saveat, maxiters = 1e8)

    t = sol.t 
    u1 = [u[1] for u in sol.u]
    u2 = [u[2] for u in sol.u]
    u3 = [u[3] for u in sol.u] 

    # Store the result in a DataFrame
    df = DataFrame(time=t, u1=u1, u2=u2, u3=u3)
    j = i # n+i if already n precomputed exist
    # Save trajectory to CSV file
    CSV.write("lorenz/trajectory_$(j).csv", df)
    
    println("Saved trajectory", i)
end

println("All trajectories are saved.")
