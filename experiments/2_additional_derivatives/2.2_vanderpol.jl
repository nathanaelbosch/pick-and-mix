using ProbNumDiffEq
using DifferentialEquations
using SciPyDiffEq
using Statistics
using LinearAlgebra
using UnPack
using DiffEqDevTools
using JLD


include("../workprecision.jl")
DIR = @__DIR__



######################################################################################
# Problem setup
######################################################################################
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_lotkavoltera, prob_ode_fitzhughnagumo, prob_ode_threebody, prob_ode_vanderpol, prob_ode_vanstiff
prob = remake(prob_ode_vanstiff, tspan=(0.0, 10.0))
appxsol = solve(remake(prob, u0=big.(prob.u0)), RadauIIA5(), abstol=1e-16, reltol=1e-16, maxiters=1e6)



######################################################################################
# Solve
######################################################################################
orders = [3,5]
fdbs = [3]

abstols = 1.0 ./ 10.0 .^ (6:11)
reltols = 1.0 ./ 10.0 .^ (3:8)

wps = Dict()
for o in orders
    wps["$o"] = MyWorkPrecision(prob, EK1(order=o), abstols, reltols)
    for f in fdbs
        wps["$o,$f"] = MyWorkPrecision(prob, EK1FDB(order=o, jac_quality=f), abstols, reltols)
    end
end

# wps["Tsit5"] = MyWorkPrecision(prob, Tsit5(), abstols, reltols; save_everystep=false, dense=false)
wps["RadauIIA5"] = MyWorkPrecision(prob, RadauIIA5(), abstols ./ 10, reltols ./ 10; save_everystep=false, dense=false)
wps["Radau-scipy"] = MyWorkPrecision(prob, SciPyDiffEq.Radau(), abstols ./ 10, reltols ./ 10; save_everystep=false, dense=false)
wps["LSODA-scipy"] = MyWorkPrecision(prob, SciPyDiffEq.LSODA(), abstols ./ 10, reltols ./ 10; save_everystep=false, dense=false)

save(joinpath(DIR, "vdp_wps.jld"), "wps", wps)
