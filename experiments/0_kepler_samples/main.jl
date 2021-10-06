using ProbNumDiffEq
using ProbNumDiffEq: X_A_Xt, stack, remake_prob_with_jac
using DifferentialEquations
using LinearAlgebra
using Plots

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

appxsol = solve(prob, Vern9(), abstol=1e-12, reltol=1e-12)



##########################################################################################
# Constant quantities
H(p, q) = norm(p)^2/2 - inv(norm(q))
L(p, q) = q[1]*p[2] - p[1]*q[2]
H(u) = H(u[3:4], u[1:2])
L(u) = L(u[3:4], u[1:2])
H2(u) = H(u[1:2], u[3:4])
L2(u) = L(u[1:2], u[3:4])


g(u) = [H2(u) - H(du0, u0); L2(u) - L(du0, u0)]









dt = 1//20
sol1 = solve(prob, EK1(order=3, diffusionmodel=:fixed), adaptive=false, dt=dt)
# sol2 = solve(prob2, EK1(order=2, diffusionmodel=:fixed), adaptive=false, dt=dt)
sol2 = solve(prob2, EK1(order=4, diffusionmodel=:fixed), adaptive=false, dt=dt,
             callback=ProbNumDiffEq.ManifoldUpdate(g; save_positions=(false, false)));
# sol1 = solve(prob, EK1(order=1))
# sol2 = solve(prob2, EK1(order=2))


Plots.plot(sol1, vars=(1,2))
Plots.plot!(sol2, vars=(3,4))


# Plots.plot(sol1, vars=(1,2), xlims=(0.35, 0.45), ylims=(-0.05, 0.05))
# Plots.plot!(sol2, vars=(3,4))


N = 100
samples1, _ = ProbNumDiffEq.dense_sample(sol1, N; density=10000)
# Plots.plot()
for i in 1:N
    Plots.plot!(samples1[:, 1, i], samples1[:, 2, i], label="", color=1, alpha=0.8, linewidth=0.1)
end
Plots.plot!()

samples2, _ = ProbNumDiffEq.dense_sample(sol2, N; density=10000)
# Plots.plot()
for i in 1:N
    Plots.plot!(samples2[:, 1, i], samples2[:, 2, i], label="", color=2, alpha=0.8, linewidth=0.1)
end
Plots.plot!()

Plots.plot!(xlims=(0.398, 0.402), ylims=(-0.02, 0.01))
# Plots.plot!(xlims=(0.35, 0.45), ylims=(-0.1, 0.1))


Plots.plot!(xlims=(-1.6005, -1.5995), ylims=(-0.01, 0.01))
# Plots.plot!(xlims=(-1.61, -1.59), ylims=(-0.05, 0.05))

##########################################################################################
# Can we see something in the residuals?
s1, s2 = stack(sol1.u), stack([u[:] for u in sol2.u])
as1, as2 = stack(appxsol.(sol1.t)), stack(appxsol.(sol2.t))
N = 10
samples1 = ProbNumDiffEq.sample(sol1, N)
samples2 = ProbNumDiffEq.sample(sol2, N)

Plots.plot(s1[:, 1] - as1[:, 1], s1[:, 2] - as1[:, 2])
Plots.plot!(s2[:, 3] - as2[:, 1], s2[:, 4] - as2[:, 2])

for i in 1:N
    Plots.plot!(samples1[:, 1, i] - as1[:, 1], samples1[:, 2, i] - as1[:, 2], color=1, label="", linewidth=0.1)
end
for i in 1:N
    Plots.plot!(samples2[:, 1, i] - as2[:, 1], samples2[:, 2, i] - as2[:, 2], color=2, label="", linewidth=0.1)
end
Plots.plot!()



##########################################################################################
# Utility
function scale_solution!(sol, scale, diffscale=scale)
    sol.diffusions .*= diffscale
    [copy!(s.Σ, ProbNumDiffEq.apply_diffusion(s.Σ, scale)) for s in sol.x_filt]
    [copy!(s.Σ, ProbNumDiffEq.apply_diffusion(s.Σ, scale)) for s in sol.x_smooth]
end
function plot_samples!(p, sol, N=3; color=:black, width=1, alpha=0.1, label="", kwargs...)
    samples = ProbNumDiffEq.sample(sol, N);
    for i in 1:N
        Plots.plot!(samples[:, 1, i], samples[:, 2, i],
              color=color, width=width, alpha=alpha, label=label, kwargs...)
    end
    return p, samples
end

N = 80
alpha = 0.5
width = 0.25
meanwidth = 2


##########################################################################################
# How do samples look WITHOUT the additional update?
@warn "Make sure to modify the solver code to not use rescaling after the solve!"
sol = solve(prob, EK0(order=1, diffusionmodel=:fixed), adaptive=false, dt=1//100)
scale_solution!(sol, 10)

p = Plots.scatter([0], [0], color=:gray, markersize=10, label="", aspect_ratio=:equal)
p, samples = plot_samples!(p, sol, N; color=:gray, width=width, alpha=alpha)
p = Plots.plot!([u[1] for u in sol.u], [u[2] for u in sol.u], color=COLORS[1], width=meanwidth,)
Plots.plot!(ylims=(-1, 1), xlims=(-2, 1), legend=false, ticks=false, axis=false)
Plots.scatter!(sol[end][1:1], sol[end][2:2], color=COLORS[1], markersize=5, label="")
Plots.plot!(size=(200,200/3*2))
savefig(joinpath(DIR, "samples_vanilla.pdf"))






##########################################################################################
# How do samples look WITH the additional update?
g(u) = [H(u) - H(du0, u0); L(u) - L(du0, u0)]

@warn "Make sure to modify the solver code to not use rescaling after the solve!"
sol_m = solve(prob, EK0(order=1, diffusionmodel=:fixed, manifold=g), adaptive=false, dt=1//70)
scale_solution!(sol_m, 20)

p_m = scatter([0], [0], color=:gray, markersize=10, label="", aspect_ratio=:equal)
p_m, samples_m = plot_samples!(p_m, sol_m, N; color=:gray, width=width, alpha=alpha)
p_m = plot!([u[1] for u in sol_m.u], [u[2] for u in sol_m.u], color=COLORS[2], width=meanwidth,)
plot!(ylims=(-1, 1), xlims=(-2, 1), legend=false, ticks=false, axis=false)
scatter!(sol_m[end][1:1], sol_m[end][2:2], color=COLORS[2], markersize=5, label="")
plot!(size=(200,200/3*2))
savefig(joinpath(DIR, "samples_informed.pdf"))
