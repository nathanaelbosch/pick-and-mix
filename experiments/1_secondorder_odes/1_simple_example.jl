using ProbNumDiffEq
using ProbNumDiffEq: stack
using DifferentialEquations
import Plots: RGB
using Statistics
using LinearAlgebra
# using GLMakie
using CairoMakie


DIR = @__DIR__


COLORS = parse.(
    RGB, ["#107D79", "#FF9933", "#1F77B4", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"]
)



#Parameters
ω = 1

#Initial Conditions
x₀ = [0.0]
dx₀ = [π/2]
tspan = (0.0, 2π)

ϕ = atan((dx₀[1]/ω)/x₀[1])
A = √(x₀[1]^2 + dx₀[1]^2)

#Define the problem
function harmonicoscillator1(du,u,ω,t)
    du[2] = u[1]
    du[1] = -ω^2 * u[2]
end
function harmonicoscillator2(ddu,du,u,ω,t)
    ddu .= -ω^2 * u
end

#Pass to solvers
prob1 = ODEProblem(harmonicoscillator1, [dx₀..., x₀...], tspan, ω)
prob2 = SecondOrderODEProblem(harmonicoscillator2, dx₀, x₀, tspan, ω)

# sol = solve(prob1, Tsit5())
# sol = solve(prob2, DPRKN6())

#Plot
# Plots.plot(sol, vars=[2,1], linewidth=2, title ="Simple Harmonic Oscillator", xaxis = "Time", yaxis = "Elongation", label = ["x" "dx"])
# Plots.plot!(t->A*cos(ω*t-ϕ), lw=3, ls=:dash, label="Analytical Solution x")
# Plots.plot!(t->-A*ω*sin(ω*t-ϕ), lw=3, ls=:dash, label="Analytical Solution dx")





true_solution(t) = A*cos(ω*t-ϕ)
true_derivative(t) = -A*ω*sin(ω*t-ϕ)
true_second_derivative(t) = -A*ω^2*cos(ω*t-ϕ)




set_theme!()
theme = Theme(
    fontsize=10,
    Axis=(
        xgridvisible=false,
        ygridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        spinewidth=0.7,
        xtickwidth=0.7,
        ytickwidth=0.7,
        xticksize=2,
        yticksize=2,
        xticklabelsize=9,
        yticklabelsize=9,
    ),
    Legend=(
        labelsize=9,
        framevisible=false,
        patchsize=(13,13),
    )
)
set_theme!(theme)

fig = Figure(
    resolution = (600, 150),
    # resolution = (250, 150),
    figure_padding = 5,
)
ax1 = fig[1, 1] = Axis(fig, title=L"Y^{(0)}", xlabel=L"t")
ax2 = fig[1, 2] = Axis(fig, title=L"Y^{(0)} - y", xlabel=L"t")
ax3 = fig[1, 3] = Axis(fig, title=L"Y^{(1)} - \dot{y}", xlabel=L"t")
ax4 = fig[1, 4] = Axis(fig, title=L"Y^{(2)} - \ddot{y}", xlabel=L"t")

lines!(ax1, 0:0.1:2π, true_solution.(0:0.1:2π), linestyle=:dash, color=:black)
hlines!(ax2, [0], linestyle=:dash, color=:black)
hlines!(ax3, [0], linestyle=:dash, color=:black)
hlines!(ax4, [0], linestyle=:dash, color=:black)




interp(sol, times) = StructArray([
    sol.interp(t, sol.t, sol.x_filt, sol.x_smooth, sol.diffusions) for t in times])
function plot_sol(
    sol, c, firstorder;
    times = tspan[1]:0.025:tspan[2]
    )

    x = interp(sol, times)
    means = stack(x.μ)
    stds = sqrt.(stack(diag.(x.Σ)))
    discrete_means = stack(sol.x_smooth.μ)

    d = length(sol.u[1])
    q = sol.interp.q
    D = d*(q+1)
    if firstorder
        means = means[:, [5, 1, 2, 3]]
        stds = stds[:, [5, 1, 2, 3]]
        discrete_means = discrete_means[:, [5, 1, 2, 3]]
    end

    true_vals_dense = [true_solution.(times) true_derivative.(times) true_second_derivative.(times)]
    true_vals = [true_solution.(sol.t) true_derivative.(sol.t) true_second_derivative.(sol.t)]

    linealpha = 0.8

    lines!(
        fig[1, 1], times, means[:, 1],
        linewidth=2.5,
        color=(COLORS[c], linealpha),
    )
    scatter!(
        fig[1, 1], sol.t, discrete_means[:, 1],
        color=COLORS[c], markercolor=COLORS[c], size=5,
        markersize=firstorder ? 12 : 8,
        marker=firstorder ? :circle : :utriangle,
        strokewidth=0.5,
    )

    for i in 1:3
        lines!(
            fig[1, i+1], times, means[:, i] - true_vals_dense[:, i],
            color=(COLORS[c],linealpha),
            linewidth=2.5,
        )
        scatter!(
            fig[1, i+1], sol.t, discrete_means[:, i] - true_vals[:, i],
            color=COLORS[c], markercolor=COLORS[c], size=5,
            markersize=firstorder ? 10 : 8,
            marker=firstorder ? :circle : :utriangle,
            strokewidth=0.5,
        )
        band!(fig[1, i+1], times,
              means[:, i] .- 1.96 .* stds[:, i] - true_vals_dense[:, i],
              means[:, i] .+ 1.96 .* stds[:, i] - true_vals_dense[:, i],
              color=(COLORS[c], 0.2))
    end
end



# PN:
N = 7
order = 4
sol1 = solve(prob1, EK1(order=order-1, diffusionmodel=:fixed),
             adaptive=false, tstops=range(tspan...; length=N))
sol2 = solve(prob2, EK1(order=order, diffusionmodel=:fixed),
             adaptive=false, tstops=range(tspan...; length=N))

function plot_samples(sol, c, firstorder; N=5)
    alpha = 0.6
    width = 0.5

    means, times = ProbNumDiffEq.dense_sample_states(sol, N)

    d = length(sol.u[1])
    q = sol.interp.q
    D = d*(q+1)
    if firstorder
        # means = means[:, D÷2+1:end, :]
        means = means[:, [5, 1, 2, 3], :]
    end

    true_vals_dense = [true_solution.(times) true_derivative.(times) true_second_derivative.(times)]
    true_vals = [true_solution.(sol.t) true_derivative.(sol.t) true_second_derivative.(sol.t)]

    for j in 1:N
        lines!(
        fig[1, 1], times, means[:, 1, j],
        linewidth=width,
        color=(COLORS[c], alpha),
        )

        for i in 1:3
            lines!(
                fig[1, i+1], times, means[:, i, j] - true_vals_dense[:, i],
                color=(COLORS[c], alpha),
                linewidth=width,
            )
        end
    end
end

# plot_samples(sol1, 1, true)
# plot_samples(sol2, 2, false)
plot_sol(sol1, 1, true)
plot_sol(sol2, 2, false)


ax1.xticks = ([0, 2π], ["0", "2π"])
ax2.xticks = ([0, 2π], ["0", "2π"])
ax3.xticks = ([0, 2π], ["0", "2π"])
ax4.xticks = ([0, 2π], ["0", "2π"])
# hidexdecorations!(ax1, grid = false)
# hidexdecorations!(ax2, grid = false)
# hidexdecorations!(ax3, grid = false)

ylims!(ax1, -2.2, 2.2)
ax1.yticks = [-2, 2]
ylims!(ax2, -0.03, 0.03)
ax2.yticks = [-0.02, 0.02]
ylims!(ax3, -0.04, 0.04)
ax3.yticks = [-0.03, 0.03]
ylims!(ax4, -0.12, 0.12)
ax4.yticks = [-0.1, 0.1]




e1 = [LineElement(color = COLORS[1], linestyle = :solid),
      MarkerElement(color = COLORS[1], marker=:circle, strokecolor = :black, strokewidth=0.5)]
e2 = [LineElement(color = COLORS[2], linestyle = :solid),
      MarkerElement(color = COLORS[2], marker=:utriangle, strokecolor = :black, strokewidth=0.5)]
e3 = LineElement(color = :black, linestyle = :dash)

leg = fig[0, :] = Legend(
    fig, [e1, e2,
          # e3
          ],
    ["Solved as a first-order ODE", "Solved as a second-order ODE",
     # "true solution"
     ],
    orientation=:horizontal,
    # tellwidth=false,
    tellheight=true,
    # framevisible=false,
    valign=:bottom,
)
rowgap!(fig.layout, 0)
colgap!(fig.layout, 15)
trim!(fig.layout)

save(joinpath(DIR, "fig2_secondorder_example.pdf"), fig, pt_per_unit=1)
