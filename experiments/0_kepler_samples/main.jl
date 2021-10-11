using ProbNumDiffEq
using ProbNumDiffEq: X_A_Xt, stack, remake_prob_with_jac
using DifferentialEquations
using LinearAlgebra
using CairoMakie

include("../theme.jl")
DIR = @__DIR__


##########################################################################################
# Problem definition
function twobody(du, u, p, t)
    R3 = norm(u[1:2])^3
    @. du[3:4] = -u[1:2] / R3
    @. du[1:2] = u[3:4]
end
u0, du0 = [0.4, 0.0], [0.0, 2.0]
# tspan = (0, 3 * 2π)
# tspan = (0, 7/8 * 2π)
# tspan = (0, 9999/10000 * 2π)
tspan = (0, 99/100 * 2π)
# tspan = (0, 2π)
prob = ODEProblem(twobody, [u0...; du0...], tspan)



function twobody2(ddu, du, u, p, t)
    R3 = norm(u)^3
    @. ddu = -u / R3
end
prob2 = SecondOrderODEProblem(twobody2, du0, u0, tspan)

# appxsol = solve(prob, Vern9(), abstol=1e-12, reltol=1e-12)



##########################################################################################
# Constant quantities
H(p, q) = norm(p)^2/2 - inv(norm(q))
L(p, q) = q[1]*p[2] - p[1]*q[2]
H(u) = H(u[3:4], u[1:2])
L(u) = L(u[3:4], u[1:2])
H2(u) = H(u[1:2], u[3:4])
L2(u) = L(u[1:2], u[3:4])


g1(u) = [H(u) - H(du0, u0); L(u) - L(du0, u0)]
g2(u) = [H2(u) - H(du0, u0); L2(u) - L(du0, u0)]









# dt = 1//5
# sol1 = solve(prob, EK1(order=2, diffusionmodel=:fixed), adaptive=false, dt=dt//2);
# sol2 = solve(prob, EK1(order=2, diffusionmodel=:fixed, manifold=g1), adaptive=false, dt=dt);
# sol3 = solve(prob2, EK1(order=3, diffusionmodel=:fixed), adaptive=false, dt=dt//2);
# sol4 = solve(prob2, EK1(order=3, diffusionmodel=:fixed, manifold=g2), adaptive=false, dt=dt);
# # sol4 = solve(prob2, EK1(order=3, diffusionmodel=:fixed), adaptive=false, dt=dt,
# #              callback=ProbNumDiffEq.ManifoldUpdate(g2; save_positions=(false, false)));
# # sol1 = solve(prob, EK1(order=2), abstol=1e-3, reltol=1e-1)
# # sol2 = solve(prob2, EK1(order=3), abstol=1e-3, reltol=1e-1,
# #              callback=ProbNumDiffEq.ManifoldUpdate(g; save_positions=(false, false)));
# # N = 100
# # samples1, times = ProbNumDiffEq.dense_sample(sol1, N)
# # samples2, times = ProbNumDiffEq.dense_sample(sol2, N)


dt = 1//12
sol1 = solve(prob, EK1(order=1, diffusionmodel=:fixed), adaptive=false, dt=dt//2);
sol2 = solve(prob, EK1(order=1, diffusionmodel=:fixed, manifold=g1), adaptive=false, dt=dt);
sol3 = solve(prob2, EK1(order=2, diffusionmodel=:fixed), adaptive=false, dt=dt//2);
sol4 = solve(prob2, EK1(order=2, diffusionmodel=:fixed, manifold=g2), adaptive=false, dt=dt);




fig = Figure(
    # resolution=(600,170)
    resolution=(250, 220),
    figure_padding=5,
)
kwargs = (
    xgridvisible=false,
    ygridvisible=false,
    rightspinevisible=false,
    topspinevisible=false,
    leftspinevisible=false,
    bottomspinevisible=false,
    # rightspinevisible=true,
    # topspinevisible=true,
    # leftspinevisible=true,
    # bottomspinevisible=true,
    aspect = DataAspect(),
    xticksvisible=false,
    xticklabelsvisible=false,
    yticksvisible=false,
    yticklabelsvisible=false,
    limits=((-2.5, 0.5), (-1.2, 1.2)),
    xlabelsize=8,
    ylabelsize=8,
    xaxisposition=:top,
)
fig[1, 1] = Axis(fig; xlabel="First-order ODE", ylabel="Conventional", kwargs...)
fig[1, 2] = Axis(fig; xlabel="Second-order ODE", kwargs...)
fig[2, 1] = Axis(fig; ylabel="Conserved energy", kwargs...)
fig[2, 2] = Axis(fig; kwargs...)


for (i, sol) in enumerate((sol1, sol2))
    N = 100
    samples1, times = ProbNumDiffEq.dense_sample(sol, N)

    lines!(
        fig[i, 1],
        [u[1] for u in appxsol.u],
        [u[2] for u in appxsol.u],
        color=:black,
        linestyle=:dash,
        linewidth=1,
    )
    for j in 1:N
        lines!(fig[i, 1], samples1[:, 1, j], samples1[:, 2, j],
               # color=(COLORS[1], 1.0),
               color=(:gray, 0.5),
               label="", alpha=0.8, linewidth=0.1)
    end
    ds1 = sol(times).u.μ
    lines!(
        fig[i, 1],
        [u[1] for u in ds1],
        [u[2] for u in ds1],
        color=(COLORS[i], 0.8),
    )
    scatter!(
        fig[i, 1],
        [u[1] for u in sol.u],
        [u[2] for u in sol.u],
        color=COLORS[i],
        markersize=2,
    )
    scatter!(fig[i, 1], [sol.u[end][1]], [sol.u[end][2]], color=COLORS[i], markersize=8)
    scatter!(fig[i, 1], [0], [0], color=:gray, markersize=15)
end


for (i, sol) in enumerate((sol3, sol4))
    N = 100
    samples2, times = ProbNumDiffEq.dense_sample(sol, N)
    lines!(
        fig[i, 2],
        [u[1] for u in appxsol.u],
        [u[2] for u in appxsol.u],
        color=:black,
        linestyle=:dash,
        linewidth=1,
    )
    for j in 1:N
        lines!(fig[i, 2], samples2[:, 1, j], samples2[:, 2, j],
               # color=(COLORS[2], 1.0),
               color=(:gray, 0.5),
               label="", alpha=0.8, linewidth=0.1)
    end
    ds2 = sol(times).u.μ
    lines!(
        fig[i, 2],
        [u[3] for u in ds2],
        [u[4] for u in ds2],
        color=(COLORS[i+2], 0.8),
    )
    scatter!(
        fig[i, 2],
        [u[3] for u in sol.u],
        [u[4] for u in sol.u],
        color=COLORS[i+2],
        strokewidth=0.5,
        markersize=2,
    )
    scatter!(fig[i, 2], [sol.u[end][3]], [sol.u[end][4]], color=COLORS[i+2], markersize=8)
    scatter!(fig[i, 2], [0], [0], color=:gray, markersize=15)
end
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
# Label(fig[1,1][1,1,TopLeft()], "A", font=noto_sans_bold, textsize=10)
# Label(fig[1,2][1,1,TopLeft()], "B", font=noto_sans_bold, textsize=10)
# Label(fig[2,1][1,1,TopLeft()], "C", font=noto_sans_bold, textsize=10)
# Label(fig[2,2][1,1,TopLeft()], "D", font=noto_sans_bold, textsize=10)
# Label(fig[1,1][1,1,Top()], "First-order ODE", textsize=8)
# Label(fig[1,2][1,1,Top()], "Second-order ODE", textsize=8)
# Label(fig[2,1][1,1,Top()], "Physical conservation", textsize=8)
# Label(fig[2,2][1,1,Top()], "Fully informed", textsize=8)

rowgap!(fig.layout, -5)
colgap!(fig.layout, -5)
trim!(fig.layout)


save(joinpath(DIR, "figure1.pdf"), fig, pt_per_unit=1)
