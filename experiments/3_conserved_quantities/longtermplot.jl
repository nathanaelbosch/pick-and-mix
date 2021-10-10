using ProbNumDiffEq
using DifferentialEquations
using Statistics
using LinearAlgebra
using UnPack
using CairoMakie


include("../theme.jl")
DIR = @__DIR__


LINEWIDTH = 0.1


######################################################################################
# Problem setup
######################################################################################
u0, du0 = [0.0, 0.1], [0.5, 0.0]
initial = [u0...; du0...]
# tspan = (0,30.)
tspan = (0, 1e3)
V(x,y) = 1//2 * (x^2 + y^2 + 2x^2*y - 2//3 * y^3)  # Potential
E(dx,dy,x,y) = V(x,y) + 1//2 * (dx^2 + dy^2);  # Total energy of the system
function Hénon_Heiles(du,u,p,t)
    x  = u[3]
    y  = u[4]
    dx = u[1]
    dy = u[2]
    du[3] = dx
    du[4] = dy
    du[1] = -x - 2x*y
    du[2] = y^2 - y -x^2
end
prob = ODEProblem(Hénon_Heiles, initial, tspan)
prob = ProbNumDiffEq.remake_prob_with_jac(prob)
E(u) = E(u...)
g(u) = [E(u) - E(initial)]


function hh2(ddu, du,u,p,t)
    ddu[1] = - u[1] - 2*u[1]*u[2]
    ddu[2] = u[2]^2 - u[2] - u[1]^2
end
prob2 = SecondOrderODEProblem(hh2, du0, u0, tspan)


# Appxsol
appxsol = solve(prob2, KahanLi8(), dt=1//100)

# EK1
sol1 = solve(prob2, EK1(order=3), abstol=1e-3, reltol=5e-3)
# sol2 = solve(prob2, EK1(order=3), abstol=1e-3, reltol=5e-3,
#              callback=ProbNumDiffEq.ManifoldUpdate(
#                  g; save_positions=(false, true), maxiters=1))
sol2 = solve(prob2, EK1(order=3, manifold=g), abstol=1e-3, reltol=5e-3)


fig = Figure(
    resolution=(600,200),
    # resolution=(6000,2000),
    figure_padding=5,
)
ax1 = fig[1, 1] = Axis(fig, aspect = DataAspect())
ax2 = fig[1, 2] = Axis(fig, aspect = DataAspect())
ax3 = fig[1, 3] = Axis(fig, aspect = DataAspect())
hidedecorations!(ax1)
hidedecorations!(ax2)
hidedecorations!(ax3)
hidespines!(ax1)
hidespines!(ax2)
hidespines!(ax3)

lines!(ax1,
       [u[3] for u in appxsol.u],
       [u[4] for u in appxsol.u],
       linewidth=LINEWIDTH,
       color=:black,
       )
lines!(ax2,
       [u[3] for u in sol1.u],
       [u[4] for u in sol1.u],
       linewidth=LINEWIDTH,
       color=PN_COLORS[1],
       )
lines!(ax3,
       [u[3] for u in sol2.u],
       [u[4] for u in sol2.u],
       linewidth=LINEWIDTH,
       color=PN_COLORS[2],
       )

colgap!(fig.layout, 0)
trim!(fig.layout)

save(joinpath(DIR, "longtermplot.pdf"), fig, pt_per_unit=1)
save(joinpath(DIR, "longtermplot.png"), fig, pt_per_unit=1)
