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
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems;
importodeproblems();
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_vanstiff
prob = remake(prob_ode_vanstiff, tspan=(0.0, 10.0))
prob = ProbNumDiffEq.remake_prob_with_jac(prob)
appxsol = solve(
    remake(prob, u0=big.(prob.u0)),
    RadauIIA5(),
    abstol=1e-16,
    reltol=1e-16,
    maxiters=1e6,
)

######################################################################################
# Solve
######################################################################################
orders = [3, 5]
fdbs = [3]

abstols = 1.0 ./ 10.0 .^ (6:11)
reltols = 1.0 ./ 10.0 .^ (3:8)

wps = Dict()
for o in orders
    wps["$o"] = MyWorkPrecision(
        prob,
        EK1(order=o, smooth=false),
        abstols,
        reltols;
        dense=false,
        save_everystep=false,
    )
    for f in fdbs
        wps["$o,$f"] = MyWorkPrecision(
            prob,
            EK1(order=o, smooth=false, fdb_improved=f),
            abstols,
            reltols;
            dense=false,
            save_everystep=false,
        )
    end
end

wps["RadauIIA5"] = MyWorkPrecision(
    prob,
    RadauIIA5(),
    abstols ./ 10,
    reltols ./ 10;
    save_everystep=false,
    dense=false,
)

save(joinpath(DIR, "vdp_wps.jld"), "wps", wps)
