using ProbNumDiffEq
using DifferentialEquations
using Statistics
using LinearAlgebra
using UnPack
using DiffEqDevTools
using JLD

using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_lotkavoltera, prob_ode_fitzhughnagumo, prob_ode_threebody, prob_ode_vanderpol, prob_ode_vanstiff

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
appxsol = TestSolution(solve(remake(prob, u0=big.(u0)), Vern9(), abstol=1e-20, reltol=1e-20))

# Problem 1: Lotka-Volterra
function lotkavolterra(du, u, p, t)
    x, y = u
    a, b, c, d = p
    du[1] = a*x - b*x*y
    du[2] = d*x*y - c*y
end
u0 = [1.0; 1.0]
tspan = (0.0, 7.0)
p = [1.5,1.0,3.0,1.0]
prob1 = ODEProblem(lotkavolterra, u0, tspan, p)



# Problem 2: Van-der-Pol
prob2 = remake(prob_ode_vanstiff, tspan=(0.0, 10.0))


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

        sol_call() = solve(prob, alg, abstol=atol, reltol=rtol, args...; kwargs...)

        sol = sol_call()
        errsol = appxtrue(sol, appxsol)

        # Time
        tbest = 100
        for i in 1:5
            t = @elapsed sol_call()
            tbest = min(tbest, t)
        end

        r = Dict(
            :final => errsol.errors[:final],
            # :l2 => errsol.errors[:l2],
            # :L2 => errsol.errors[:L2],
            :time => tbest,
            :nf => isnothing(errsol.destats) ? nothing : errsol.destats.nf,
            :njacs => isnothing(errsol.destats) ? nothing : errsol.destats.njacs,
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



# prob = prob_ode_threebody
prob = prob_ode_fitzhughnagumo
# prob = prob_ode_vanderpol
# prob = prob_ode_vanstiff
# prob = ProbNumDiffEq.remake_prob_with_jac(prob1)
# appxsol = TestSolution(solve(remake(prob, u0=big.(prob.u0)), Vern9(), abstol=1e-20, reltol=1e-20))
prob = ProbNumDiffEq.remake_prob_with_jac(prob)
appxsol = TestSolution(solve(remake(prob, u0=big.(prob.u0)), RadauIIA5(), abstol=1e-20, reltol=1e-20))


orders = [3,5]
abstols = 1.0 ./ 10.0 .^ (6:11)
reltols = 1.0 ./ 10.0 .^ (3:8)
wps = Dict()
for o in orders
    wps["EK1($o)"] = MyWorkPrecision(prob, EK1(order=o, smooth=false), abstols, reltols,
                                     save_everystep=false, dense=false)
    wps["EK1FDB($o)"] = MyWorkPrecision(prob, EK1FDB(order=o, smooth=false), abstols, reltols,
                                        save_everystep=false, dense=false)
end
# wps["Tsit5"] = MyWorkPrecision(prob, Tsit5(), abstols ./ 100, reltols ./ 100;
#                                save_everystep=false, dense=false)
wps["RadauIIA5"] = MyWorkPrecision(prob, RadauIIA5(), abstols ./ 100, reltols ./ 100;
                                   save_everystep=false, dense=false)
# wps["RK45 (SciPy)"] = MyWorkPrecision(prob, SciPyDiffEq.RK45(), abstols ./ 100, reltols ./ 100;
#                                save_everystep=false, dense=false)
wps["LSODA (SciPy)"] = MyWorkPrecision(prob, SciPyDiffEq.LSODA(), abstols ./ 100, reltols ./ 100;
                                      save_everystep=false, dense=false)
