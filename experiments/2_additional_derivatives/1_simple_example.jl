using ProbNumDiffEq
using ProbNumDiffEq: stack
using Statistics
using LinearAlgebra
using CairoMakie
using LaTeXStrings
using ForwardDiff

include("../theme.jl")
DIR = @__DIR__

update_theme!(Axis=(xgridvisible=false, ygridvisible=false))

# Problem definition
logistic(du, u, p, t) = (@. du = p * u * (1 - u))
analytic(u0, p, t) = [exp(p[1] * t) / (1 / u0[1] - 1 + exp(p[1] * t))]
u0 = [1e-2]
tspan = (0.0, 3.0)
p = [3.0]
prob = ODEProblem(ODEFunction(logistic, analytic=analytic), u0, tspan, p)
prob = ProbNumDiffEq.remake_prob_with_jac(prob)

# True solution and derivatives
true_solution(t) = analytic(u0, p, t)[1]
true_derivative(t) =
    p[1] * u0[1] * exp(p[1] * t) * (1 - u0[1]) / (1 + u0[1] * (exp(p[1] * t) - 1))^2
true_second_derivative(t) = ForwardDiff.derivative(true_derivative, t)

# Plot them
prange = tspan[1]:0.1:tspan[2]
fig = Figure(
    # resolution = (600, 150),
    resolution=(250, 300),
    figure_padding=5,
)

ax1 =
    fig[1, 1] = Axis(fig, ylabel=L"Y^{(0)}", xticksvisible=false, xticklabelsvisible=false)
ax2 =
    fig[2, 1] =
        Axis(fig, ylabel=L"Y^{(0)} - y", xticksvisible=false, xticklabelsvisible=false)
ax3 =
    fig[3, 1] = Axis(
        fig,
        ylabel=L"Y^{(1)} - \dot{y}",
        xticksvisible=false,
        xticklabelsvisible=false,
    )
ax4 = fig[4, 1] = Axis(fig, ylabel=L"Y^{(2)} - \ddot{y}", xlabel=L"t")

lines!(ax1, 0:0.1:3, true_solution.(0:0.1:3), linestyle=:dash, color=:black)
hlines!(ax2, [0], linestyle=:dash, color=:black)
hlines!(ax3, [0], linestyle=:dash, color=:black)
hlines!(ax4, [0], linestyle=:dash, color=:black)

# Solve with and without the additional update!
interp(sol, times) = StructArray([
    sol.interp(t, sol.t, sol.x_filt, sol.x_smooth, sol.diffusions) for t in times
])
function plot_sol(
    sol,
    c;
    color=COLORS[c],
    marker=[:circle, :utriangle][c],
    axes=(ax1, ax2, ax3, ax4),
    times=tspan[1]:0.025:tspan[2],
)
    x = interp(sol, times)
    means = stack(x.μ)
    stds = sqrt.(stack(diag.(x.Σ)))

    true_vals_dense =
        [true_solution.(times) true_derivative.(times) true_second_derivative.(times)]
    true_vals =
        [true_solution.(sol.t) true_derivative.(sol.t) true_second_derivative.(sol.t)]

    ax1, ax2, ax3, ax4 = axes

    lines!(ax1, times, means[:, 1], linewidth=2.5, color=(color, 0.7))
    scatter!(
        ax1,
        sol.t,
        stack(sol.x_smooth.μ)[:, 1],
        color=color,
        markercolor=color,
        size=5,
        markersize=c == 1 ? 10 : 8,
        marker=marker,
        strokewidth=0.5,
    )

    for (i, ax) in enumerate((ax2, ax3, ax4))
        lines!(
            ax,
            times,
            means[:, i] - true_vals_dense[:, i],
            linewidth=2.5,
            color=(color, 0.7),
        )
        scatter!(
            ax,
            sol.t,
            stack(sol.x_smooth.μ)[:, i] - true_vals[:, i],
            color=color,
            markercolor=color,
            size=5,
            markersize=c == 1 ? 10 : 8,
            marker=marker,
            strokewidth=0.5,
        )
        band!(
            ax,
            times,
            means[:, i] .- 1.96 .* stds[:, i] - true_vals_dense[:, i],
            means[:, i] .+ 1.96 .* stds[:, i] - true_vals_dense[:, i],
            color=(color, 0.2),
        )
    end
end

N = 8
order = 3
sol1 = solve(
    prob,
    EK1(order=order, diffusionmodel=:fixed),
    adaptive=false,
    tstops=range(tspan...; length=N),
)
sol2 = solve(
    prob,
    EK1(order=order, diffusionmodel=:fixed, fdb_improved=3),
    adaptive=false,
    tstops=range(tspan...; length=N),
)
plot_sol(sol1, 1)
plot_sol(sol2, 2)

e1 = [
    LineElement(color=COLORS[1], linestyle=:solid),
    MarkerElement(color=COLORS[1], marker=:circle, strokecolor=:black, strokewidth=0.5),
]
e2 = [
    LineElement(color=COLORS[2], linestyle=:solid),
    MarkerElement(color=COLORS[2], marker=:utriangle, strokecolor=:black, strokewidth=0.5),
]
e3 = LineElement(color=:black, linestyle=:dash)

leg =
    fig[0, :] = Legend(
        fig,
        [e1, e2],
        ["Without ÿ", "With ÿ"],
        orientation=:horizontal,
        tellwidth=false,
        tellheight=true,
        framevisible=false,
    )
trim!(fig.layout)
fig

ax1.xticks = [0, 3]
ax2.xticks = [0, 3]
ax3.xticks = [0, 3]
ax4.xticks = [0, 3]

ylims!(ax2, -0.0125, 0.0125)
ax2.yticks = [-0.01, 0.01]
ylims!(ax3, -0.025, 0.025)
ax3.yticks = [-0.02, 0.02]
ylims!(ax4, -0.25, 0.25)
ax4.yticks = [-0.2, 0.2]
ylims!(ax1, -0.2, 1.2)
ax1.yticks = [0, 1]

rowgap!(fig.layout, 10)
rowgap!(fig.layout, 1, 0)
colgap!(fig.layout, 15)
trim!(fig.layout)

# If we'd want to save just this figure, uncomment the following line
# save(joinpath(DIR, "fig4.pdf"), fig, pt_per_unit=1)
