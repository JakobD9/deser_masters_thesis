using Pkg
Pkg.activate("ChaoticTipping")
using LaTeXStrings
using Plots
include("utils.jl")
include("parameters.jl")

paramset = 2
param = create_param(paramset)
# Choose range for R
r_range = 1:0.01:50

res = zeros(length(r_range))
for i in 1:lastindex(r_range)
    radius = r_range[i]
    num_points = floor(Int64, 20*radius)
    tuples = tuple_circle(radius, num_points)

    xtAMOC = zeros(num_points)
    for i in 1:num_points
        xtAMOC[i] = xtf(tuples[i], param)
    end
    res[i] = maximum(xtAMOC)
end

pic = scatter(r_range, res, markersize=1, markercolor=:black, title="Dissipativity condition AMOC", xlabel= L"R", ylabel="max " * L"(x, y)f(x, y)" * " on " * L"S_R(0)", label = missing)
current_dir = @__DIR__
savefig(pic, joinpath(current_dir, "other_pics",  "AMOC_dissipativity.png"))

# Picture of sampled points on the lorenz attractor
#=
x, p = sample_lorenz_attractor(50)
display(p)
savefig(p, joinpath(current_dir, "other_pics", "lorenz_sample.png"))
=#

# Picture of double well potential
#=
p = 1.3
x = range(-1.5, 1.5, 1000)
f(x) = (x-1)^2*(x+1)^2 -p*x
y = f.(x)

points_x = [-1, 0, 1]
points_y = f.(points_x)
labels = [L"\bar x_1", L" x_s", L"\bar x_2"]

# Create the plot
p2 = plot(x, y, label=L"V", legend=:top, xlabel=L"x", title="Double well potential")
scatter!(points_x, points_y, color=:red, label=missing)
annotate!([(px, py+0.3, l) for (px, py, l) in zip(points_x, points_y, labels)])
savefig(p2, joinpath(current_dir, "other_pics", "double_well.png"))
=#