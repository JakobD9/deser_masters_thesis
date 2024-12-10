using Pkg
#=
Pkg.generate("ChaoticTipping")
Pkg.activate("ChaoticTipping")
Pkg.add("CSV")
Pkg.add("ColorTypes")
Pkg.add("DataFrames")
Pkg.add("DelimitedFiles")
Pkg.add("DifferentialEquations")
Pkg.add("Distributed")
Pkg.add("Integrals")
Pkg.add("LaTeXStrings")
Pkg.add("LinearAlgebra")
Pkg.add("NLsolve")
Pkg.add("PrettyTables")
Pkg.add("Plots")
Pkg.add("QuadGK")
Pkg.add("Random")
Pkg.add("Roots")
Pkg.add("SharedArrays")
Pkg.add("Statistics")
=#

using ColorTypes
using DelimitedFiles
using DifferentialEquations
using LaTeXStrings
using Plots
using StatsBase


# AMOC strength q
function q(S_N, p, t)
    p, t_p = p
    return p.lambda * (p.alpha * (p.T_S - p.T_0 - (t_p.forcing ? forcing_pwl(t, t_p) : 0)) .+ (p.beta / 100) .* (S_N .- p.S_S)) ./ (1 + p.lambda * p.alpha * p.mu)
end

function q(S_N, S_S, p, t)
    p, t_p = p
    return p.lambda * (p.alpha * (p.T_S - p.T_0 - (t_p.forcing ? forcing_pwl(t, t_p) : 0)) .+ (p.beta / 100) .* (S_N .- S_S)) ./ (1 + p.lambda * p.alpha * p.mu)
end

# q without T-forcing
function q0(S_N, p)
    return p.lambda .* (p.alpha .* (p.T_S - p.T_0) .+ (p.beta / 100) .* (S_N .- p.S_S)) ./ (1 + p.lambda * p.alpha * p.mu)
end

# out of place model formulation
function three_box(u, p, t)
    p, h_p, t_p, lp = p
    u1, u2 = u
    q_val = q(u1, (p, t_p), t)
    du1 = ifelse(q_val >= 0, f1p(u1, u2, p, h_p, t_p, t), f1n(u1, u2, p, h_p, t_p, t))
    du2 = ifelse(q_val >= 0, f2p(u1, u2, p, h_p, t_p, t), f2n(u1, u2, p, h_p, t_p, t))
    return [du1, du2]
end

# model function: Chaos only in S_N
function three_box_chaotic(u, p, t)
    p, h_p, t_p, l_p = p
    u1, u2, x1, x2, x3 = u
    q_val = q(u1, (p, t_p), t)

    du1 = sqrt(2*l_p.delta)* 1/l_p.epsilon * 1/l_p.sigma * x1 + 
        ifelse(q_val >= 0, f1p(u1, u2, p, h_p, t_p, t), f1n(u1, u2, p, h_p, t_p, t))

    du2 = 0*sqrt(2*l_p.delta)* 1/l_p.epsilon * 1/l_p.sigma * x1 + # only chaos in S_N
        ifelse(q_val >= 0, f2p(u1, u2, p, h_p, t_p, t), f2n(u1, u2, p, h_p, t_p, t))
    dx1 = l_p.mu*(x2-x1)
    dx2 = x1*(l_p.rho-x3) - x2
    dx3 = x1*x2 - l_p.beta*x3
    return [du1, du2, (1/l_p.epsilon^2)*dx1, (1/l_p.epsilon^2)*dx2, (1/l_p.epsilon^2)*dx3]
end

# chaos matrix B
function three_box_coupled_chaos(u, p, t)
    p, h_p, t_p, l_p = p
    u1, u2, x1, x2, x3, x4, x5, x6 = u
    q_val = q(u1, (p, t_p), t)
    du1 = sqrt(2*l_p.delta)* 1/l_p.epsilon * 1/l_p.sigma * (p.B11*x1 + 0*x4) + 
        ifelse(q_val >= 0, f1p(u1, u2, p, h_p, t_p, t), f1n(u1, u2, p, h_p, t_p, t))

    du2 = sqrt(2*l_p.delta)* 1/l_p.epsilon * 1/l_p.sigma * (p.B21*x1 + p.B22*x4) + 
        ifelse(q_val >= 0, f2p(u1, u2, p, h_p, t_p, t), f2n(u1, u2, p, h_p, t_p, t))
    dx1 = l_p.mu*(x2-x1)
    dx2 = x1*(l_p.rho-x3) - x2
    dx3 = x1*x2 - l_p.beta*x3
    dx4 = l_p.mu*(x5-x4)
    dx5 = x4*(l_p.rho-x6) - x5
    dx6 = x4*x5 - l_p.beta*x6
    return [du1, du2, (1/l_p.epsilon^2) *dx1, (1/l_p.epsilon^2)*dx2, (1/l_p.epsilon^2)*dx3, (1/l_p.epsilon^2)*dx4, (1/l_p.epsilon^2)*dx5, (1/l_p.epsilon^2)*dx6]
end

function lorenz(u, p, t)
    l_p = p
    vel = 1/l_p.epsilon^2
    x1, x2, x3 = u
    dx1 = l_p.mu*(x2-x1)
    dx2 = x1*(l_p.rho-x3) - x2
    dx3 = x1*x2 - l_p.beta*x3
    return [vel*dx1, vel*dx2, vel*dx3]
end

function lorenz!(du, u, p, t)
    l_p = p
    du[1] = l_p.mu * (u[2] - u[1])
    du[2] = u[1] * (l_p.rho - u[3]) - u[2]
    du[3] = u[1] * u[2] - l_p.beta * u[3]
end

function sample_lorenz_attractor(num_points, t_end=5000+100*num_points, integr_steps = 500000)
    initial_condition = [1.6,-2.6,26.6].+ randn(3).*0.0001
    t_start = 0
    dt_ic=(t_end - t_start) / (integr_steps - 1)
    times = t_start:dt_ic:t_end
    prob = ODEProblem(lorenz, initial_condition, (t_start, t_end), lorenz_parameters(28, 10, 8/3, 1, 0.5, 0))
    traj = solve(prob, Tsit5(), saveat=times, maxiters = 1e7)

    indices = floor.(Int64, sample(0.2*integr_steps:integr_steps, num_points, replace=false)) # Choose random indizes ofter the first fifth of the time span
    attr_points = traj[indices].+randn(3, num_points).*0.0001
    
    p1 = plot3d(traj[1, :], traj[2, :], traj[3, :], label="trajectory", title="Sampling of initial contitions", xlabel=L"y_1", ylabel=L"y_2", zlabel=L"y_3")
    x = attr_points[1, :]
    y = attr_points[2, :]
    z = attr_points[3, :]
    
    scatter3d!(p1, x, y, z, color = :red, label = "samples")

    x = sort(indices)
    y = zeros(Int64, num_points-1)
    for i in 1:num_points-1
        y[i] = x[i+1]-x[i]
    end
    println("Minimum time distance of sampled points: ", minimum(y)*dt_ic)
    return attr_points, p1
end

function five_box(u, p, t)
    p, h_p, t_p, lp = p
    u1, u2, u3, u4 = u
    q_val = q(u1, u3, (p, t_p), t)
    aq = abs(q_val)
    # S_N
    du1 = ifelse(q_val >= 0,
                 p.Y/p.V_N * (q_val * (u2 - u1) + p.K_N * (u2 - u1) - 100 * (p.F_N + p.A_N*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0),
                 p.Y/p.V_N * (aq * (SB(u1, u2, u3, u4, p) - u1) + p.K_N * (u2 - u1) - 100 * (p.F_N + p.A_N*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0)) 
    # S_T
    du2 = ifelse(q_val >= 0,
                 p.Y/p.V_T * (q_val * (p.gamma * u3 + (1 - p.gamma) * (u4) - u2) + p.K_S * (u3 - u2) + p.K_N * (u1 - u2) - 100 * (p.F_T + p.A_T*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0),
                 p.Y/p.V_T * (aq * (u1 - u2) + p.K_S * (u3 - u2) + p.K_N * (u1 - u2) - 100 * (p.F_T + p.A_T*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0))
    # S_S
    du3 = ifelse(q_val >= 0,
                 p.Y/p.V_S * (p.gamma * q_val * (SB(u1, u2, u3, u4, p) - u3) + p.K_IP * (u4 - u3) + p.K_S * (u2 - u3) + p.eta * (SB(u1, u2, u3, u4, p)-u3)- 100 * (p.F_S + p.A_S*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0),
                 p.Y/p.V_S * (p.gamma * aq * (u2 - u3) + p.K_IP * (u4 - u3) + p.K_S * (u2 - u3) + p.eta * (SB(u1, u2, u3, u4, p)-u3)- 100 * (p.F_S + p.A_S*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0)) 
    # S_IP
    du4 = ifelse(q_val >= 0,
                 p.Y/p.V_IP * ((1-p.gamma)*q_val * (SB(u1, u2, u3, u4, p) - u4) + p.K_IP * (u3 - u4) - 100 * (p.F_IP + p.A_IP*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0),
                 p.Y/p.V_IP * ((1-p.gamma)*aq * (u2- u4) + p.K_IP * (u3 - u4) - 100 * (p.F_IP + p.A_IP*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0)) 
    return [du1, du2, du3, du4]
end

function five_box_chaotic(u, p, t) # not implemented
    error("Not implemented")
    return NaN
end

# Calculate S_IP from other variables
SIP(S_N, S_T, p) = (100 * (p.C - p.S_0*(p.V_B+p.V_N+p.V_T+p.V_IP+p.V_S)) - (p.V_N* S_N + p.V_T* S_T+ p.V_S *p.S_S + p.V_B* p.S_B) ) / p.V_IP

# Calculate S_B from other variables (5-box)
SB(S_N, S_T, S_S, S_IP, p) = (100 * (p.C - p.S_0*(p.V_B+p.V_N+p.V_T+p.V_IP+p.V_S)) - (p.V_N* S_N + p.V_T* S_T+ p.V_S *S_S + p.V_IP* S_IP) ) / p.V_B

# Define SDE functions as drift and diffusion term
# double well potential
function dw_drift!(du,u,p,t)
    du[1] = -4u[1].^3 .+4 .*u[1]
end

# constant diffusion, i.e. sqrt(2δ)
function one_d_diffusion!(du,u,p,t)
    du[1] = p[1]
end

function two_d_diffusion!(du, u, p, t)
    du[1, 1] = p[2].B11 * p[1]
    du[1, 2] = 0
    du[2, 1] = p[2].B21 * p[1]
    du[2, 2] = p[2].B22 * p[1]
end

# unused, potentially usefull when handing proc = CorrelatedWienerProcess!(Γ, 0.0, zeros(2)) to sde problem
function two_d_diffusion2!(du, u, p, t)
    du[1] = p[1]
    du[2] = p[1]
end

# Define double well with lorenz trajectory as parameter
function dw_lorenz_precomputed!(du, u, p, t)
    # Assuming p[1] = sqrt(2δ), p[2]  = lorenz parameter, p[3] = saveat, p[4] = 1st coordinate of lorenz traj
    # lorenz data has step size saveat -> time / saveat and speed 1/epsilon^2 -> time / epsilon^2
    l_p = p[2]
    index = Int(floor((1/p[3])*t / l_p.epsilon^2)) + 1
    du[1] = -4 * u[1]^3 + 4 * u[1] + p[1]/l_p.sigma * (1/l_p.epsilon) * p[4][index]
end

# Double well potential with coupled lorenz63 system
function dw_lorenz!(du, u, p, t)
    u, x1, x2, x3 = u
    l_p = p[2]
    du[1] = -4 * u^3 + 4 * u + p[1]/l_p.sigma * (1/l_p.epsilon) *  x1
    du[2] = (1/l_p.epsilon^2)*(l_p.mu * (x2 - x1))
    du[3] = (1/l_p.epsilon^2)*(x1 * (l_p.rho - x3) - x2)
    du[4] = (1/l_p.epsilon^2)*(x1 * x2 - l_p.beta * x3)
end

# 3-box model as deterministic drift
function three_box!(du, u, p, t)
    p = p[2]
    u1, u2 = u
    q_val = q0(u1, p)
    aq = abs(q_val)
    du[1] = ifelse(q_val >= 0, f1p(u1, u2, p), f1n(u1, u2, p))
    du[2] = ifelse(q_val >= 0, f2p(u1, u2, p), f2n(u1, u2, p))
end

# 3-box model with two coupled lorenz63 system
function three_box_lorenz!(du, u, p, t)
    del_sc = p[1]
    l_p = p[2]
    p = p[3]
    u1, u2, x1, x2, x3, x4, x5, x6 = u
    q_val = q0(u1, p)
    du1 = du[1] = (del_sc * 1/l_p.epsilon * 1/l_p.sigma * (p.B11*x1 + 0*x4) + 
        ifelse(q_val >= 0, f1p(u1, u2, p), f1n(u1, u2, p)))
    du2 = du[2] = (del_sc * 1/l_p.epsilon * 1/l_p.sigma * (p.B21*x1 + p.B22*x4) + 
        ifelse(q_val >= 0, f2p(u1, u2, p), f2n(u1, u2, p)))
    du[3] = l_p.mu*(x2-x1)/l_p.epsilon^2
    du[4] = (x1*(l_p.rho-x3) - x2)/l_p.epsilon^2
    du[5] = (x1*x2 - l_p.beta*x3)/l_p.epsilon^2
    du[6] = l_p.mu*(x5-x4)/l_p.epsilon^2
    du[7] = (x4*(l_p.rho-x6) - x5)/l_p.epsilon^2
    du[8] = (x4*x5 - l_p.beta*x6)/l_p.epsilon^2
end

# Define 3-box model with lorenz trajectories as parameter
function three_box_lorenz_precomputed!(du, u, p, t)
    # Assuming p[1] = sqrt(2δ), p[2]  = lorenz parameter, p[3] = saveat, p[4] = 1st coordinate of lorenz traj1, p[5] = 1st coordinate of lorenz traj2, p6 = param
    # lorenz data has step size saveat -> time / saveat and speed 1/epsilon^2 -> time / epsilon^2
    del_sc = p[1]
    l_p = p[2]
    x1 = p[4]
    x4 = p[5]
    index = Int(floor((1/p[3])*t / l_p.epsilon^2)) + 1
    p = p[6]
    u1, u2 = u
    q_val = q0(u1, p)
    du[1] = (del_sc * 1/l_p.epsilon * 1/l_p.sigma * (p.B11*x1[index]) + 
        ifelse(q_val >= 0, f1p(u1, u2, p), f1n(u1, u2, p)))

    du[2] = (del_sc * 1/l_p.epsilon * 1/l_p.sigma * (p.B21*x1[index] + p.B22*x4[index]) + 
        ifelse(q_val >= 0, f2p(u1, u2, p), f2n(u1, u2, p)))
end

# homogeneous AMOC equations

f1p(u1, u2, p) = p.Y / p.V_N * (q0(u1, p) * (u2 - u1) + p.K_N * (u2 - u1) - 100 * p.F_N * p.S_0)

f2p(u1, u2, p) = p.Y/p.V_T * (q0(u1, p) * (p.gamma * p.S_S + (1 - p.gamma) * (SIP(u1, u2, p)) - u2) + p.K_S * (p.S_S - u2) + p.K_N * (u1 - u2) - 100 * p.F_T * p.S_0)

f1n(u1, u2, p) = p.Y / p.V_N * (-q0(u1, p) * (p.S_B - u1) + p.K_N * (u2 - u1) - 100 * p.F_N * p.S_0)

f2n(u1, u2, p) = p.Y / p.V_T * (-q0(u1, p) * (u1 - u2) + p.K_S * (p.S_S - u2) + p.K_N * (u1 - u2) - 100 * p.F_T * p.S_0)

function fp(x, p)
    u1, u2 = x
    return [f1p(u1, u2, p), f2p(u1, u2, p)]
end

function fn(x, p)
    u1, u2 = x
    return [f1n(u1, u2, p), f2n(u1, u2, p)]
end

function f(x, p)
    u1, u2 = x
    return q0(u1, p) >= 0 ? fp(x, p) : fn(x, p)
end

function xtf(x, p)
    u1, u2 = x
    return q0(u1, p) >= 0 ? u1* f1p(u1, u2, p) + u2* f2p(u1, u2, p) : u1*f1n(u1, u2, p) +u2*f2n(u1, u2, p)
end

# inhomogeneous AMOC equations

f1p(u1, u2, p, h_p, t_p, t) = p.Y/p.V_N * (q(u1, (p, t_p), t) * (u2 - u1) + p.K_N * (u2 - u1) - 100 * (p.F_N + p.A_N*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0)

f2p(u1, u2, p, h_p, t_p, t) = p.Y/p.V_T * (q(u1, (p, t_p), t) * (p.gamma * p.S_S + (1 - p.gamma) * (SIP(u1, u2, p)) - u2) + p.K_S * (p.S_S - u2) 
    + p.K_N * (u1 - u2) - 100 * (p.F_T + p.A_T*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0)

f1n(u1, u2, p, h_p, t_p, t) = p.Y/p.V_N * (-q(u1, (p, t_p), t) * (p.S_B - u1) + p.K_N * (u2 - u1) - 100 * (p.F_N + p.A_N*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0)

f2n(u1, u2, p, h_p, t_p, t) = p.Y/p.V_T * (-q(u1, (p, t_p), t) * (u1 - u2) + p.K_S * (p.S_S - u2) + p.K_N * (u1 - u2) 
    - 100 * (p.F_T + p.A_T*(h_p.forcing ? forcing_pwl(t, h_p) : 0)) * p.S_0)

# forced 3-box model as deterministic drift
function three_box_f!(du, u, p, t)
    h_p = p[3]
    t_p = p[4]
    p = p[2]
    u1, u2 = u
    q_val = q(u1, (p, t_p), t)
    aq = abs(q_val)
    du[1] = ifelse(q_val >= 0, f1p(u1, u2, p, h_p, t_p, t), f1n(u1, u2, p, h_p, t_p, t))
    du[2] = ifelse(q_val >= 0, f2p(u1, u2, p, h_p, t_p, t), f2n(u1, u2, p, h_p, t_p, t))
end

# time dependent piecewise linear forcing (hosing or temperature ascent)
function forcing_pwl(t, f_p)
    t = t-f_p.t_f
    r_rise = (- f_p.initial_forcing + f_p.height_forcing) / (f_p.t_rise+0.00001)
    r_fall = (f_p.final_forcing - f_p.height_forcing) / (f_p.t_fall+0.00001)
    a(t) = r_rise * t + f_p.initial_forcing
    b(t) = r_fall * (t - f_p.t_rise - f_p.t_pert - f_p.t_fall) + f_p.final_forcing
    
    if 0 <= t <= f_p.t_rise
        f_pwl = a(t)
    elseif f_p.t_rise <= t <= f_p.t_rise + f_p.t_pert
        f_pwl = f_p.height_forcing
    elseif f_p.t_rise + f_p.t_pert <= t <= f_p.t_rise + f_p.t_pert + f_p.t_fall
        f_pwl = b(t)
    elseif 0 > t
        f_pwl = f_p.initial_forcing
    else
        f_pwl = f_p.final_forcing
    end
    return f_pwl
end

function color_gradient(color1, color2, num_colors)
    return [RGB(
        (1 - i / num_colors) * red(color1) + (i / num_colors) * red(color2),
        (1 - i / num_colors) * green(color1) + (i / num_colors) * green(color2),
        (1 - i / num_colors) * blue(color1) + (i / num_colors) * blue(color2)
    ) for i in 1:num_colors]
end

function solve_ode(model, initial_condition, t_span, dt, param)
    prob = ODEProblem(model, initial_condition, t_span, param)
    sol = solve(prob, Tsit5(), saveat=dt)
    return sol
end

# returns the first argument that is neither missing, nor nothing
function coalesce(args...)
    for arg in args
        if !ismissing(arg) && !isnothing(arg)
            return arg
        end
    end
    return nothing
end

# returns n points on a circle of radius r
function tuple_circle(r, n)
    points = Vector{}()
    for i in 0:(n-1)
        angle = 2 * π * i / n
        x = r * cos(angle)
        y = r * sin(angle)
        push!(points, (x, y))
    end
    return points
end