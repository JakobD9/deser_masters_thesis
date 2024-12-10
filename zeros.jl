using Pkg
Pkg.activate("ChaoticTipping")
using NLsolve
using Roots
using Plots
include("parameters.jl")
include("utils.jl")

# Define parameters
paramset = 8
param = create_param(paramset)

x_crit = find_zero(x -> q0(x, param), 0)
println("x_crit: ", x_crit)

sol_fp = nlsolve(x -> fp(x, param), [x_crit + 0.5, 0.1])
println("Solution for fp (u1 > x_crit): ", sol_fp.zero[1] > x_crit ? sol_fp.zero : "Not in range")

sol_fn = nlsolve(x -> fn(x, param), [x_crit - 0.5, 0.1])
println("Solution for fn (u1 < x_crit): ", sol_fn.zero[1] < x_crit ? sol_fn.zero : "Not in range")

# Create contour plots
u1_range = -1:0.01:1
u2_range = -1:0.01:1

# Contour plots for f1p, f2p where q(u1) >= 0
c1 = contour(u1_range, u2_range, (u1, u2) -> f1p(u1, u2, param), levels=20, title="f1p Contour")
c2 = contour(u1_range, u2_range, (u1, u2) -> f2p(u1, u2, param), levels=20, title="f2p Contour")
plot(c1, c2, layout=(2,1), size=(800, 800))

# Contour plots for f1n, f2n where q(u1) < 0
c3 = contour(u1_range, u2_range, (u1, u2) -> f1n(u1, u2, param), levels=20, title="f1n Contour")
c4 = contour(u1_range, u2_range, (u1, u2) -> f2n(u1, u2, param), levels=20, title="f2n Contour")
plot(c3, c4, layout=(2,1), size=(800, 800))


print("estimation plus: ")
println(fp(sol_fp.zero, param))
print("estimation minus: ")
println(fn(sol_fn.zero, param))
