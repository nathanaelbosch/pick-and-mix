using ProbNumDiffEq
using CairoMakie
using ColorSchemes


include("../theme.jl")
DIR = @__DIR__


######################################################################################
# Problem
######################################################################################
function rober_mm(du,u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    du[1] = -k₁*y₁ + k₃*y₂*y₃
    du[2] =  k₁*y₁ - k₃*y₂*y₃ - k₂*y₂^2
    du[3] =  y₁ + y₂ + y₃ - 1
    nothing
end

M = [1. 0  0
     0  1. 0
     0  0  0]
f_mm = ODEFunction(rober_mm, mass_matrix=M)
u0 = [1.0,0.0,0.0]
# tspan = (0.0, 1e2)
tspan = (0.0, 1e5)
p = (0.04,3e7,1e4)
prob_mm = ODEProblem(f_mm, u0, tspan, p)


######################################################################################
# Solve
######################################################################################
sol = solve(prob_mm, EK1())
# appxsol = solve(remake(prob_mm, u0=big.(prob_mm.u0), tspan=big.(prob_mm.tspan), p=big.(prob_mm.p)),
#                 RadauIIA5(), reltol=1e-18, abstol=1e-18)


######################################################################################
# Plot
######################################################################################
# fig = Figure(
#     resolution=(100,100),
#     figure_padding=5,
# )
# ax1 = fig[1, 1] # = Axis(fig, xscale=log10, xticklabelsvisible=false)
# ax2 = fig[2, 1] # = Axis(fig, xscale=log10, xticklabelsvisible=false)
# ax3 = fig[3, 1] # = Axis(fig, xscale=log10, xlabel=L"t")

xlims!(p0ax1, 1e-5, 1e5)
xlims!(p0ax2, 1e-5, 1e5)
xlims!(p0ax3, 1e-5, 1e5)

lines!(p0ax1, sol.t, [u[1] for u in sol.u])
lines!(p0ax2, sol.t, [u[2] for u in sol.u])
lines!(p0ax3, sol.t, [u[3] for u in sol.u])

p0ax1.yticks = [0, 1]
p0ax2.yticks = ([0, 3e-5], ["0", "3⋅10⁻⁵"])
p0ax3.yticks = [0, 1]

rowgap!(fig.layout, 5)

save(joinpath(DIR, "rober.pdf"), fig, pt_per_unit=1)
