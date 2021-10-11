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
function pendulum!(du, u, p, t)
    x, dx, y, dy, T = u
    g, L = p
    du[1] = dx; du[2] = T*x
    du[3] = dy; du[4] = T*y - g
    du[5] = x^2 + y^2 - L^2
end
pendulum_fun! = ODEFunction(pendulum!, mass_matrix=Matrix(Diagonal([1,1,1,1,0])))
u0 = [1.0, 0, 0, 0, 0]; p = [9.8, 1]; tspan = (0, 10.0)
pendulum_prob = ODEProblem(pendulum_fun!, u0, tspan, p)
# solve(pendulum_prob,Rodas4())


# Transform it
traced_sys = modelingtoolkitize(pendulum_prob)
pendulum_sys = structural_simplify(dae_index_lowering(traced_sys))
prob = ODEProblem(pendulum_sys, Pair[], tspan)
bigprob = remake(prob, u0=big.(prob.u0), tspan=big.(prob.tspan), p=big.(prob.p))
sol = solve(prob, Rodas4())
# appxsol = solve(bigprob, RadauIIA5(), abstol=1e-14, reltol=1e-14)
appxsol = solve(bigprob, Rodas4(), abstol=1e-16, reltol=1e-16)


######################################################################################
# Work-Precision Diagrams
######################################################################################

abstols = 1.0 ./ 10.0 .^ (5:12)
reltols = 1.0 ./ 10.0 .^ (2:9)

wps = Dict()
wps["Rodas5"] = MyWorkPrecision(prob, Rodas5(), abstols ./ 100, reltols ./ 100; appxsol=appxsol)
wps["Rosenbrock23"] = MyWorkPrecision(prob, Rosenbrock23(), abstols[1:end-1], reltols[1:end-1]; appxsol=appxsol)
wps["RadauIIA5"] = MyWorkPrecision(prob, RadauIIA5(), abstols ./ 100, reltols ./ 100; appxsol=appxsol)
wps["QNDF"] = MyWorkPrecision(prob, QNDF(), abstols ./ 100, reltols ./ 100; appxsol=appxsol)
for (i, o) in enumerate((2, 3, 5))
    wps["EK1($o)"] = MyWorkPrecision(
        prob, EK1(order=o, smooth=false),
        abstols, reltols; appxsol=appxsol, dense=false, save_everystep=false)
    # wps["EK1($o)-index4"] = MyWorkPrecision(
    #     pendulum_prob, EK1(order=o, smooth=false),
    #     abstols, reltols; appxsol=appxsol, dense=false, save_everystep=false)
end

save(joinpath(DIR, "pendulum_wps.jld"), "wps", wps)
