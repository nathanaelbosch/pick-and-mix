using ProbNumDiffEq
using ProbNumDiffEq: stack
using OrdinaryDiffEq
using SciPyDiffEq
using Plots: RGB
using BenchmarkTools
using DiffEqDevTools
using LinearAlgebra
using JLD

include("../workprecision.jl")
DIR = @__DIR__

@fastmath function pleiades(du, u, p, t)
    v = view(u, 1:7)   # x
    w = view(u, 8:14)  # y
    x = view(u, 15:21) # x′
    y = view(u, 22:28) # y′
    du[15:21] .= v
    du[22:28] .= w
    @inbounds @simd ivdep for i in 1:14
        du[i] = zero(eltype(u))
    end
    @inbounds @simd ivdep for i in 1:7
        @inbounds @simd ivdep for j in 1:7
            if i != j
                r = ((x[i] - x[j])^2 + (y[i] - y[j])^2)^(3 / 2)
                du[i] += j * (x[j] - x[i]) / r
                du[7+i] += j * (y[j] - y[i]) / r
            end
        end
    end
end
x0 = [3.0, 3.0, -1.0, -3.0, 2.0, -2.0, 2.0]
y0 = [3.0, -3.0, 2.0, 0, 0, -4.0, 4.0]
dx0 = [0, 0, 0, 0, 0, 1.75, -1.5]
dy0 = [0, 0, 0, -1.25, 1, 0, 0]
u0 = [dx0; dy0; x0; y0]
tspan = (0.0, 3.0)
prob1 = ODEProblem(pleiades, u0, tspan)
appxsol = solve(remake(prob1, u0=big.(prob1.u0)), Vern9(), abstol=1e-20, reltol=1e-20)

function pleiades2(ddu, du, u, p, t)
    x = view(u, 1:7)
    y = view(u, 8:14)
    for i in 1:14
        ddu[i] = zero(eltype(u))
    end
    for i in 1:7, j in 1:7
        if i != j
            r = ((x[i] - x[j])^2 + (y[i] - y[j])^2)^(3 / 2)
            ddu[i] += j * (x[j] - x[i]) / r
            ddu[7+i] += j * (y[j] - y[i]) / r
        end
    end
end
u0 = [x0; y0]
du0 = [dx0; dy0]
prob2 = SecondOrderODEProblem(pleiades2, du0, u0, tspan)


######################################################################################
# Work-Precision Diagrams
######################################################################################
wps = Dict()

abstols = 1.0 ./ 10.0 .^ (6:11)
reltols = 1.0 ./ 10.0 .^ (3:8)

wps["Tsit5"] = MyWorkPrecision(
    prob1,
    Tsit5(),
    abstols ./ 10,
    reltols ./ 10;
    dense=false,
    save_everystep=false,
)
wps["RadauIIA5"] = MyWorkPrecision(
    prob1,
    RadauIIA5(),
    abstols ./ 10,
    reltols ./ 10;
    dense=false,
    save_everystep=false,
)
wps["Vern6"] = MyWorkPrecision(
    prob1,
    Vern6(),
    abstols ./ 10,
    reltols ./ 10;
    dense=false,
    save_everystep=false,
)
wps["DPRKN6"] =
    MyWorkPrecision(prob2, DPRKN6(), abstols, reltols; dense=false, save_everystep=false)
wps["RK45-SciPy"] = MyWorkPrecision(
    prob1,
    SciPyDiffEq.RK45(),
    abstols,
    reltols;
    dense=false,
    save_everystep=false,
)
wps["LSODA-SciPy"] = MyWorkPrecision(
    prob1,
    SciPyDiffEq.LSODA(),
    abstols,
    reltols;
    dense=false,
    save_everystep=false,
)

smooth = false
wps["ODE1;o=3;EK0"] = MyWorkPrecision(
    prob1,
    EK0(order=2, smooth=smooth),
    abstols,
    reltols;
    dense=smooth,
    save_everystep=false,
)
wps["ODE2;o=3;EK0"] = MyWorkPrecision(
    prob2,
    EK0(order=3, smooth=smooth),
    abstols,
    reltols;
    dense=smooth,
    save_everystep=false,
)
wps["ODE1;o=5;EK1"] = MyWorkPrecision(
    prob1,
    EK1(order=4, smooth=smooth),
    abstols,
    reltols;
    dense=smooth,
    save_everystep=false,
)
wps["ODE2;o=5;EK1"] = MyWorkPrecision(
    prob2,
    EK1(order=5, smooth=smooth),
    abstols,
    reltols;
    dense=smooth,
    save_everystep=false,
)

save(joinpath(DIR, "workprecisiondata.jld"), "wps", wps)
