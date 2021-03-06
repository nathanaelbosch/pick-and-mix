using DifferentialEquations
using ProbNumDiffEq
using LinearAlgebra
using Statistics
using DiffEqDevTools
using ModelingToolkit
using JLD

include("../workprecision.jl")
DIR = @__DIR__


######################################################################################
# Problem
######################################################################################
function rober_mm(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end

M = [
    1.0 0 0
    0 1.0 0
    0 0 0
]
f_mm = ODEFunction(rober_mm, mass_matrix=M)
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 1e2)
p = (0.04, 3e7, 1e4)
prob_mm = ODEProblem(f_mm, u0, tspan, p)

function rober_ode(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    return nothing
end
prob_ode = ODEProblem(rober_ode, u0, tspan, p)

solve_ref(prob) = solve(
    remake(prob, u0=big.(prob.u0), tspan=big.(prob.tspan), p=big.(prob.p)),
    RadauIIA5(),
    reltol=1e-18,
    abstol=1e-18,
)
appxsol_mm = appxsol = solve_ref(prob_mm);
appxsol_ode = solve_ref(prob_ode);

# Conservation law:
g(u) = [u[1] + u[2] + u[3] - 1]
solve(prob_ode, EK1(manifold=g))


######################################################################################
# Work-Precision Diagrams
######################################################################################
wps = Dict()

_prob, _appxsol = prob_mm, appxsol_mm

abstols = 1.0 ./ 10.0 .^ (5:9)
reltols = 1.0 ./ 10.0 .^ (2:6)

wps["Rodas5"] =
    MyWorkPrecision(_prob, Rodas5(), abstols ./ 1e3, reltols ./ 1e3; appxsol=_appxsol)
wps["Rodas4"] =
    MyWorkPrecision(_prob, Rodas4(), abstols ./ 1e3, reltols ./ 1e3; appxsol=_appxsol)
wps["Rosenbrock23"] =
    MyWorkPrecision(_prob, Rosenbrock23(), abstols ./ 1e3, reltols ./ 1e3; appxsol=_appxsol)
wps["Radau"] =
    MyWorkPrecision(_prob, RadauIIA5(), abstols ./ 1e3, reltols ./ 1e3; appxsol=_appxsol)
wps["QNDF"] =
    MyWorkPrecision(_prob, QNDF(), abstols ./ 1e3, reltols ./ 1e3; appxsol=_appxsol)

for o in (2, 3, 5)
    wps["EK1($o)"] = MyWorkPrecision(
        _prob,
        EK1(order=o, smooth=false),
        abstols,
        reltols;
        appxsol=_appxsol,
        dense=false,
        save_everystep=false,
    )
end

save(joinpath(DIR, "rober_wps.jld"), "wps", wps)
