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
tspan = (0, 99/100 * 2π)
prob = ODEProblem(twobody, [u0...; du0...], tspan)

function twobody2(ddu, du, u, p, t)
    R3 = norm(u)^3
    @. ddu = -u / R3
end
prob2 = SecondOrderODEProblem(twobody2, du0, u0, tspan)




fig = Figure(resolution=(300, 115), figure_padding=0)

for j in 1:2
    fig[1, j] = Axis(fig,
                     xgridvisible=false,
                     ygridvisible=false,
                     rightspinevisible=false,
                     topspinevisible=false,
                     leftspinevisible=false,
                     bottomspinevisible=false,
                     aspect = DataAspect(),
                     xticksvisible=false,
                     xticklabelsvisible=false,
                     yticksvisible=false,
                     yticklabelsvisible=false,
                     limits=((-1.9, 0.7), (-1.0, 1.0))
                     )
end





##########################################################################################
# Utility
function scale_solution!(sol, scale, diffscale=scale)
    sol.diffusions .*= diffscale
    [copy!(s.Σ, ProbNumDiffEq.apply_diffusion(s.Σ, scale)) for s in sol.x_filt]
    [copy!(s.Σ, ProbNumDiffEq.apply_diffusion(s.Σ, scale)) for s in sol.x_smooth]
end
function plot_samples!(ax, sol, N=3; color=:black, width=1, alpha=0.1, label="", kwargs...)
    samples = ProbNumDiffEq.sample(sol, N);
    for i in 1:N
        lines!(ax, samples[:, 1, i], samples[:, 2, i],
              color=(color, alpha), linewidth=width, label=label, kwargs...)
    end
    return samples
end

N = 80
alpha = 0.5
width = 0.1
meanwidth = 2


ALG = EK0
ORDER, SCALE = 2, 100
# ORDER, SCALE = 2, 500
DT = 1//50
RESCALE = 10

@info "Experiment information:" ALG ORDER DT SCALE informed_solver_scale=SCALE*RESCALE



##########################################################################################
# How do samples look WITHOUT the additional update?
# @warn "Make sure to modify the solver code to not use rescaling after the solve!"
# sol = solve(prob, ALG(order=ORDER, diffusionmodel=:fixed), adaptive=false, dt=DT)
sol = solve(prob2, ALG(order=ORDER+1, diffusionmodel=:fixed), adaptive=false, dt=DT)
# sol = solve(prob, ALG(order=ORDER))
@info "Uninformed nf" sol.destats.nf
scale_solution!(sol, SCALE)

scatter!(fig[1,1], [0], [0], color=:gray, markersize=10, label="", aspect_ratio=:equal)
# plot!(appxsol, vars=(3,4), #label="True solution",
#       color=:black)
samples = plot_samples!(fig[1,1], sol, N;
                        # color=COLORS[1],
                        color=:gray,
                        width=width,
                        alpha=alpha
                        )
p = lines!(fig[1,1],
           # [u[1] for u in sol.u], [u[2] for u in sol.u],
           [u[3] for u in sol.u], [u[4] for u in sol.u],
          color=COLORS[1],
          # color=:black,
          # label="EK0",
          linewidth=meanwidth,
          # marker=:o, markersize=1
          )
scatter!(fig[1,1],
         # sol[end][1:1], sol[end][2:2],
         sol[end][3:3], sol[end][4:4],
         color=COLORS[1], markersize=5, label="")
# plot!(size=(200,200/3*2))
# savefig(joinpath(FILEDIR, "samples_vanilla.pdf"))




##########################################################################################
# Constant quantities
H(p, q) = norm(p)^2/2 - inv(norm(q))
L(p, q) = q[1]*p[2] - p[1]*q[2]
H(u) = H(u[3:4], u[1:2])
L(u) = L(u[3:4], u[1:2])
H2(u) = H(u[1:2], u[3:4])
L2(u) = L(u[1:2], u[3:4])
# plot(appxsol.t, H2.(appxsol.u) .- H(du0, u0))
# plot!(appxsol.t, L2.(appxsol.u) .- L(du0, u0))
# plot!(sol_euler.t, H2.(sol_euler.u) .- H(du0, u0))
# plot!(sol_euler.t, L2.(sol_euler.u) .- L(du0, u0))
# plot!(sol_ek0.t, H.(sol_ek0.u) .- H(du0, u0))
# plot!(sol_ek0.t, L.(sol_ek0.u) .- L(du0, u0))



# How do samples look WITH the additional update?
g1(u) = [H(u) - H(du0, u0); L(u) - L(du0, u0)]
g2(u) = [H2(u) - H(du0, u0); L2(u) - L(du0, u0)]

# sol_m = solve(prob, ALG(order=ORDER),
#               callback=ProbNumDiffEq.ManifoldUpdate(g1; save_positions=(false, false), maxiters=1));
# sol_m = solve(prob, ALG(order=ORDER, diffusionmodel=:fixed), adaptive=false, dt=DT,
#               callback=ProbNumDiffEq.ManifoldUpdate(g1; save_positions=(false, false)));
sol_m = solve(prob2, ALG(order=ORDER+1, diffusionmodel=:fixed), adaptive=false, dt=DT,
              callback=ProbNumDiffEq.ManifoldUpdate(g2; save_positions=(false, false), maxiters=1));
@info "informed nf" sol_m.destats.nf
scale_solution!(sol_m, SCALE*RESCALE)

scatter!(fig[1,2],[0], [0], color=:gray, markersize=10, label="", aspect_ratio=:equal)
# plot!(appxsol, vars=(3,4), #label="True solution",
#       color=:black)
samples_m = plot_samples!(fig[1,2], sol_m, N;
                               # color=COLORS[2],
                               color=:gray,
                               width=width, alpha=alpha)
lines!(fig[1,2],
       # [u[1] for u in sol_m.u], [u[2] for u in sol_m.u],
       [u[3] for u in sol_m.u], [u[4] for u in sol_m.u],
       color=COLORS[2], linewidth=meanwidth,
       )
# scatter!(fig[1,2], sol_m[end][1:1], sol_m[end][2:2], color=COLORS[2], markersize=5, label="")
scatter!(fig[1,2], sol_m[end][3:3], sol_m[end][4:4], color=COLORS[2], markersize=5, label="")
# plot!(size=(200,200/3*2))


rowgap!(fig.layout, 0)
colgap!(fig.layout, 0)
trim!(fig.layout)


save(joinpath(DIR, "inflated_samples.pdf"), fig, pt_per_unit=1)
