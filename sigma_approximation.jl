using Pkg
Pkg.activate("ChaoticTipping")
using DelimitedFiles
using DifferentialEquations
using Distributed
using Integrals
using LaTeXStrings
using Plots
using QuadGK
using SharedArrays
using Statistics
@everywhere using QuadGK, DifferentialEquations  # Ensure all workers have access to these packages, only for sampled = true

include("parameters.jl")
include("utils.jl")

sampled = false # true -> one integral and sample false -> double integral over one trajectory

num_ic = 20 # For sampling

initial_time = 0.0
final_time = 5000
num_points = 100*final_time
t_span = (initial_time, final_time)
dt=(final_time - initial_time) / (num_points - 1)
model = lorenz
initial_conditions, p = sample_lorenz_attractor(num_ic)

#####
# Don't change
vel = 1 # 1 is standard lorenz time
strength = 1 # 
l_param = lorenz_parameters(28, 10, 8/3, 1/sqrt(vel), 1/2 * strength^2, 0)
#####

#display(p)
u0 = 0
ode_prob = ODEProblem(lorenz!, u0, t_span, p=l_param)

function prob_func(prob, i, repeat)
    new_u0 = initial_conditions[:, i]
    remake(prob, u0 = new_u0)
end

ensemble_prob = EnsembleProblem(ode_prob, prob_func = prob_func)
ensemble_solution = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=num_ic, maxiters = 1e7, saveat=dt, p=l_param)


# Precompute values that are reused
y0_values =  [ensemble_solution[i](0)[1] for i in 1:num_ic]

# Define a new integrant function that avoids repeated computation
function integrant_precomputed(t, i, p)
    return ensemble_solution[i](t)[1] * y0_values[i]
end

if(sampled)
    res = SharedVector{Float64}(num_ic)
    # Calculate integral over every Lorenz trajectory
    @distributed for i in 1:num_ic
        prob = IntegralProblem((s, _) -> integrant_precomputed(s, i, p), 0, final_time)
        res[i] = solve(prob, QuadGKJL())[1]
    end
    
    println(res)
    #println("final_time: ", final_time)
    println("mean: ", mean(res))
    println("variance: ", var(res))
    println("max: ", maximum(res))
    
else
    f_0_proj(x, y, z) = x

    # Funtion for arbitrary f_0
    function integrant(t, i, f_0)
        sol = ensemble_solution[i]
        return f_0(sol(t)[1], sol(t)[2], sol(t)[3])
    end

    function integrant(t, s, i, p)
        # f_0 is projection to first coordinate 
        return ensemble_solution[i](t)[1]*ensemble_solution[i](t+s)[1]
    end
    T_end = 1000
    inner_time = final_time-T_end-1
    #println(inner_time)
    res = zeros(T_end+1) #res[T+1] gives integral up to T
    for T in 1:T_end
        println(T)
        integral, error = quadgk(s -> quadgk(t -> integrant(s, t, 1, 0), 0, inner_time)[1], T-1, T)
        res[T+1] = res[T]+integral
    end

    for T in 1:T_end
        res[T+1] = res[T+1]/T
    end
    plot(0:T_end, res, label=missing)

    sigmasq = cov(res[:, end])./t
    println("Approximated sigma: ", sigmasq)
end

current_dir = @__DIR__
#savefig(p_limit, joinpath(current_dir, "Lorenz_variance.png"))
#display(p_limit)
#println(res[1, end-10:end])

#=
# Check if integral is calculate correctly as sum 
prob = IntegralProblem((s, p) -> integrant(s, 1, p), 0.0, n_range[end]*t, p=0)
x = (1/sqrt(n_range[end])) * solve(prob, QuadGKJL())[1] # integral(s -> integrant(s, i), 0, n*t)
println("Steps: ", res[1, ind(n_range[end])])
println("Complete: ", x)
=#