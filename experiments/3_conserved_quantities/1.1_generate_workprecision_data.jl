using ProbNumDiffEq
using DifferentialEquations
import Plots: RGB
using Statistics
using LinearAlgebra
using UnPack
using DiffEqDevTools
using JLD

DIR = @__DIR__
include("../workprecision.jl")
# Here we need a slightly different function in order to save the error that the
# solution has in the conserved quantity:
function MyWorkPrecision(
    prob,
    alg,
    abstols,
    reltols,
    args...;
    dts=nothing,
    appxsol=appxsol,
    dense_errors=false,
    timeseries_errors=false,
    kwargs...,
)
    results = []
    for (i, (atol, rtol)) in enumerate(zip(abstols, reltols))
        @info "info" alg atol rtol

        sol_call() =
            isnothing(dts) ?
            solve(prob, alg, abstol=atol, reltol=rtol, args...; kwargs...) :
            solve(
                prob,
                alg,
                abstol=atol,
                reltol=rtol,
                adaptive=false,
                tstops=prob.tspan[1]:dts[i]:prob.tspan[2],
                args...;
                kwargs...,
            )

        sol = sol_call()
        errsol =
            sol isa ProbNumDiffEq.ProbODESolution ?
            appxtrue(
                mean(sol),
                appxsol;
                dense_errors=dense_errors,
                timeseries_errors=timeseries_errors,
            ) :
            appxtrue(
                sol,
                appxsol;
                dense_errors=dense_errors,
                timeseries_errors=timeseries_errors,
            )

        # Time
        tbest = 9999999999
        for i in 1:3
            t = @elapsed sol_call()
            tbest = min(tbest, t)
        end

        r = Dict(
            :final => errsol.errors[:final],
            :time => tbest,
            :nf => errsol.destats.nf,
            :njacs => errsol.destats.njacs,
            :e_final => E(sol[end]) - E(sol[1]),
        )

        push!(results, r)
    end
    return results
end

######################################################################################
# Problem setup
######################################################################################
u0, du0 = [0.0, 0.1], [0.5, 0.0]
initial = [du0...; u0...]
tspan = (0, 1e2)
@info "timespan:" tspan
V(x, y) = 1 // 2 * (x^2 + y^2 + 2x^2 * y - 2 // 3 * y^3)  # Potential
E(dx, dy, x, y) = V(x, y) + 1 // 2 * (dx^2 + dy^2);  # Total energy of the system
E(u) = E(u...)
g(u, E=E, initial=initial) = [E(u) - E(initial)]

function Hénon_Heiles(du, u, p, t)
    dx = u[1]
    dy = u[2]
    x = u[3]
    y = u[4]
    du[1] = -x - 2x * y
    du[2] = y^2 - y - x^2
    du[3] = dx
    return du[4] = dy
end
prob = ODEProblem(Hénon_Heiles, initial, tspan)
prob = ProbNumDiffEq.remake_prob_with_jac(prob)

function hh2(ddu, du, u, p, t)
    ddu[1] = -u[1] - 2 * u[1] * u[2]
    return ddu[2] = u[2]^2 - u[2] - u[1]^2
end
prob2 = SecondOrderODEProblem(hh2, du0, u0, tspan)

# Compute reference solution
appxsol = solve(
    SecondOrderODEProblem(hh2, big.(du0), big.(u0), tspan),
    Vern9(),
    abstol=1e-20,
    reltol=1e-20,
)

######################################################################################
# Compute the data for the work-precision diagrams
######################################################################################
abstols = 1.0 ./ 10.0 .^ (3:11)
reltols = 1.0 ./ 10.0 .^ (0:8)

wps = Dict()
wps["Tsit5"] = MyWorkPrecision(
    prob,
    Tsit5(),
    abstols ./ 1e3,
    reltols ./ 1e3,
    maxiters=1e6,
    dense=false,
    save_everystep=false,
)
wps["Tsit5 w/ g"] = MyWorkPrecision(
    prob,
    Tsit5(),
    abstols ./ 1e3,
    reltols ./ 1e3;
    maxiters=1e6,
    dense=false,
    save_everystep=false,
    callback=ManifoldProjection((resid, u, p, t) -> (resid .= g(u))),
)
wps["DPRKN6"] = MyWorkPrecision(
    prob2,
    DPRKN6(),
    abstols,
    reltols,
    maxiters=1e6,
    dense=false,
    save_everystep=false,
)
dts = 1.0 ./ 10.0 .^ (-0.25:0.25:1.0)
wps["KahanLi8"] = MyWorkPrecision(
    prob2,
    KahanLi8(),
    dts,
    dts;
    dts=dts,
    maxiters=1e6,
    dense=false,
    save_everystep=false,
)
wps["SymplecticEuler"] = MyWorkPrecision(
    prob2,
    SymplecticEuler(),
    dts,
    dts;
    dts=dts,
    maxiters=1e6,
    dense=false,
    save_everystep=false,
)

for o in (3, 5, 8)
    wps["EK1($o)"] = MyWorkPrecision(
        prob2,
        EK1(order=o, smooth=false),
        maxiters=1e6,
        abstols,
        reltols;
        dense=false,
        save_everystep=false,
    )
    wps["EK1($o) w/ g"] = MyWorkPrecision(
        prob2,
        EK1(order=o, smooth=false, manifold=g),
        maxiters=1e6,
        abstols,
        reltols;
        dense=false,
        save_everystep=false,
    )
end

save(joinpath(DIR, "henonheiles.jld"), "wps", wps)
