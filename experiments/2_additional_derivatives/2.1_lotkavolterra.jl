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
function lotkavolterra(du, u, p, t)
    x, y = u
    a, b, c, d = p
    du[1] = a*x - b*x*y
    du[2] = d*x*y - c*y
end
u0 = [1.0; 1.0]
tspan = (0.0, 7.0)
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(lotkavolterra, u0, tspan, p)
prob = ProbNumDiffEq.remake_prob_with_jac(prob)
appxsol = solve(remake(prob, u0=big.(prob.u0)), Vern9(), abstol=1e-20, reltol=1e-20, maxiters=1e6)


######################################################################################
# Solve
######################################################################################
orders = [3,5]
fdbs = [1,2,3]

abstols = 1.0 ./ 10.0 .^ (6:11)
reltols = 1.0 ./ 10.0 .^ (3:8)

wps = Dict()
for o in orders
    wps["$o"] = MyWorkPrecision(prob, EK1(order=o, smooth=false), abstols, reltols; dense=false, save_everystep=false)
    for f in fdbs
        wps["$o,$f"] = MyWorkPrecision(prob, EK1FDB(order=o, smooth=false, jac_quality=f), abstols,
                                       reltols; dense=false, save_everystep=false)
    end
end

wps["Tsit5"] = MyWorkPrecision(prob, Tsit5(), abstols ./ 100, reltols ./ 100; save_everystep=false, dense=false)
wps["RadauIIA5"] = MyWorkPrecision(prob, RadauIIA5(), abstols ./ 100, reltols ./ 100; save_everystep=false, dense=false)
# wps["RK45-scipy"] = MyWorkPrecision(prob, SciPyDiffEq.RK45(), abstols ./ 100, reltols ./ 100; save_everystep=false, dense=false)

save(joinpath(DIR, "lv_wps.jld"), "wps", wps)
