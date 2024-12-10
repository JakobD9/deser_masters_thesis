using Pkg
Pkg.activate("ChaoticTipping")
using DelimitedFiles
using DifferentialEquations
using Integrals
using LaTeXStrings
using Plots
using Statistics
include("parameters.jl")
include("utils.jl")

num_ic = 500

initial_time = 0.0
final_time = 10000
t_span = (initial_time, final_time)
dt=0.01
initial_conditions, p = sample_lorenz_attractor(num_ic)

#display(p)
#####
# Don't change
vel = 1 # 1 is standard lorenz time
strength = 1 # 
l_param = lorenz_parameters(28, 10, 8/3, 1/sqrt(vel), 1/2 * strength^2, sqrt(60))
#####

u0 = 0
ode_prob = ODEProblem(lorenz!, u0, t_span, p=l_param)

function prob_func(prob, i, repeat)
    new_u0 = initial_conditions[:, i]
    remake(prob, u0 = new_u0)
end

ensemble_prob = EnsembleProblem(ode_prob, prob_func = prob_func)
ensemble_solution = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=num_ic, maxiters = 1e7, saveat=dt, p=l_param)

#f_0_proj(x, y, z) = x

# Funtion for arbitrary f_0
#function integrant(t, i, f_0)
    #sol = traj_all[i]
    #return f_0(sol(t)[1], sol(t)[2], sol(t)[3])
#end

t = 1.0 # fix a t > 0, only consider W_k(t) -> \sqrt(\Sigma)W_t 
k_step = 50
k_range = 1+k_step:k_step:floor(Int64, final_time/t)

res = zeros((num_ic, 1+length(k_range)))

ind(k) = 1 + floor.(Int, (k-0.5)/k_step)

# Calculate W_k for all initial conditions
for i in 1:num_ic
    println(i)
    prob = IntegralProblem((s, _) -> ensemble_solution[i](s)[1], 0.0, t)
    res[i, 1] =solve(prob, QuadGKJL())[1]
    for k in k_range
        # Integrate next intervall and add it to partial result
        prob = IntegralProblem((s, _) -> ensemble_solution[i](s)[1], (k-k_step)*t, k*t)
        res[i, ind(k)] = res[i, ind(k)-1] + solve(prob, QuadGKJL())[1] # integral(t -> integrant(t, i), 0, k*t)
    end
    for k in k_range
        res[i, ind(k)] /= sqrt(k)
    end
end

current_dir = @__DIR__
# read file
res = readdlm(joinpath(current_dir, "Green_Kubo", "green_kubo_est_500.csv"), ';')

sigmasq = cov(res[:, end])./t
println("Estimated sigma^2: ", sigmasq)

p_limit = plot(ylabel=L"W_k(" * string(t) * ")", xlabel = "k")
# plot(p_limit, 1:100:n_range[end], res[1, 1:100:end], label=missing)
for i in 1:num_ic
    plot!(p_limit, 1:k_step:k_range[end], res[i, :], label=missing)
end

savefig(p_limit, joinpath(current_dir, "Green_Kubo/Lorenz_variance_" * string(num_ic) * ".png"))
display(p_limit)
#println(res[1, end-10:end])

# Historgram
pr(x) = 1/sqrt(2pi*sigmasq) * exp(-x^2/(2*sigmasq))
bin_edges = range(-30, 30, length=20)
histogram(res[:, end], bins=bin_edges, normalize=:pdf)
plot!(pr, label="estimated", lw=3, color=:red)
savefig(joinpath(current_dir, "Green_Kubo/Lorenz_histo_" * string(num_ic) * ".png"))

filename = "Green_Kubo/green_kubo_est_" * string(num_ic) * ".csv"
writedlm(joinpath(current_dir, filename), res, ';')

sigma_t = zeros(1+length(k_range))
for i in 1:(1+length(k_range))
    sigma_t[i] = cov(res[:, i])./t
end

p_time = plot(1:k_step:k_range[end], sigma_t, ylabel=L"\hat \sigma^2", xlabel = L"k", label = missing)
savefig(joinpath(current_dir, "Green_Kubo/green_kubo_est_sequence_" * string(num_ic) * ".png"))

# Check if integral is calculate correctly as sum 
#=
prob = IntegralProblem((s, _) -> ensemble_solution[1](s)[1], 0.0, n_range[end]*t)
x = (1/sqrt(k_range[end])) * solve(prob, QuadGKJL())[1] # integral(s -> ensemble_solution[1](s)[1], 0, k*t)
println("Steps: ", res[1, ind(k_range[end])])
println("Complete: ", x)
=#