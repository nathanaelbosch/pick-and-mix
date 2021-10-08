using ProbNumDiffEq
using DifferentialEquations
using Statistics
using LinearAlgebra
using UnPack
using DiffEqDevTools
using JLD


DIR = @__DIR__


######################################################################################
# Problem setup
######################################################################################
logistic(du, u, p, t) = (@. du = p * u * (1 - u))
analytic(u0, p, t) = [exp(p[1]*t) / (1/u0[1] - 1 + exp(p[1]*t))]
u0 = [1e-2]
tspan = (0.0, 3.0)
p = [3.0]
prob = ODEProblem(ODEFunction(logistic, analytic=analytic), u0, tspan, p)
prob = ProbNumDiffEq.remake_prob_with_jac(prob)
appxsol = solve(remake(prob, u0=big.(u0)), Vern9(), abstol=1e-20, reltol=1e-20)

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
# Simple work-precision
######################################################################################
function chi2(gaussian_estimate, actual_value)
    μ, Σ = gaussian_estimate
    diff = μ - actual_value
    R = qr(Σ.squareroot').R

    @assert R'R ≈ Matrix(Σ)

    chi2_pinv = diff' * pinv(Matrix(Σ)) * diff
    return chi2_pinv
end
function MyWorkPrecision(prob, alg, abstols, reltols, args...;
                         dts=nothing, appxsol=appxsol, kwargs...)
    results = []
    for (i, (atol, rtol)) in enumerate(zip(abstols, reltols))
        @info "info" alg atol rtol

        sol_call() = isnothing(dts) ?
            solve(prob, alg, abstol=atol, reltol=rtol, args...; kwargs...) :
            solve(prob, alg, abstol=atol, reltol=rtol,
                  adaptive=false, tstops=prob.tspan[1]:dts[i]:prob.tspan[2],
                  args...; kwargs...)

        sol = sol_call()
        errsol = sol isa ProbNumDiffEq.ProbODESolution ?
            appxtrue(mean(sol), appxsol; dense_errors=true) : appxtrue(sol, appxsol)

        # Time
        tbest = 100
        for i in 1:5
            t = @elapsed sol_call()
            tbest = min(tbest, t)
        end

        r = Dict(
            :final => errsol.errors[:final],
            :l2 => errsol.errors[:l2],
            :L2 => errsol.errors[:L2],
            :time => tbest,
            :nf => errsol.destats.nf,
            :njacs => errsol.destats.njacs,
            # :chi2_final => chi2(sol.pu[end], appxsol.(sol.t[end])),
            # :final_std => std(sol.pu[end]),
        )
        if sol isa ProbNumDiffEq.ProbODESolution
            r[:chi2_final] = chi2(sol.pu[end], appxsol.(sol.t[end]))[1]
            # r[:final_std] = std(sol.pu[end])
        end

        push!(results, r)
    end
    return results
end





orders = [3,5]
fdbs = [1,2,3]


abstols = 1.0 ./ 10.0 .^ (6:11)
reltols = 1.0 ./ 10.0 .^ (3:8)


# dts = 1.0 ./ 10.0 .^ (0.0:0.5:2)
# abstols = dts
# reltols = dts
alg_kwargs = (
    diffusionmodel=:fixed,
)
wps = Dict()
for o in orders
    wps["$o"] = MyWorkPrecision(
        prob, EK1(; alg_kwargs..., order=o),
        abstols, reltols)
    for f in fdbs
        wps["$o,$f"] = MyWorkPrecision(
            prob, EK1FDB(; alg_kwargs..., order=o, jac_quality=f),
            abstols, reltols)
    end
end



wps["Tsit5"] = MyWorkPrecision(prob, Tsit5(), abstols, reltols;
                               # dts=dts,
                               save_everystep=false, dense=false)
wps["RadauIIA5"] = MyWorkPrecision(prob, RadauIIA5(), abstols, reltols;
                                   # dts=dts,
                                   save_everystep=false, dense=false)


save(joinpath(DIR, "logistic_workprecisiondata.jld"), "wps", wps)
