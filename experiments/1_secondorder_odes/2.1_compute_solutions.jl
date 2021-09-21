using ProbNumDiffEq
using ProbNumDiffEq: stack
using OrdinaryDiffEq
using Plots: RGB
using BenchmarkTools
using DiffEqDevTools
using LinearAlgebra


DIR = @__DIR__


function pleiades(du,u,p,t)
    v = view(u,1:7)   # x
    w = view(u,8:14)  # y
    x = view(u,15:21) # x′
    y = view(u,22:28) # y′
    du[15:21] .= v
    du[22:28].= w
    for i in 1:14
        du[i] = zero(eltype(u))
    end
    for i=1:7,j=1:7
        if i != j
            r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
            du[i] += j*(x[j] - x[i])/r
            du[7+i] += j*(y[j] - y[i])/r
        end
    end
end
x0 = [3.0,3.0,-1.0,-3.0,2.0,-2.0,2.0]
y0 = [3.0,-3.0,2.0,0,0,-4.0,4.0]
dx0 = [0,0,0,0,0,1.75,-1.5]
dy0 = [0,0,0,-1.25,1,0,0]
u0 = [dx0; dy0; x0; y0]
tspan = (0.0, 3.0)
prob1 = ODEProblem(pleiades, u0, tspan)
appxsol = solve(remake(prob1, u0=big.(prob1.u0)), Vern9(), abstol=1e-20, reltol=1e-20)


function pleiades2(ddu, du, u, p, t)
    x = view(u,1:7)
    y = view(u,8:14)
    for i in 1:14
        ddu[i] = zero(eltype(u))
    end
    for i=1:7,j=1:7
        if i != j
            r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
            ddu[i] += j*(x[j] - x[i])/r
            ddu[7+i] += j*(y[j] - y[i])/r
        end
    end
end
u0 = [x0; y0]
du0 = [dx0; dy0]
prob2 = SecondOrderODEProblem(pleiades2, du0, u0, tspan)



######################################################################################
# WP
######################################################################################
function chi2(gaussian_estimate, actual_value)
    μ, Σ = gaussian_estimate
    diff = μ - actual_value
    # R = qr(Σ.squareroot').R

    # @assert R'R ≈ Matrix(Σ)

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
            appxtrue(mean(sol), appxsol;
                     dense_errors=false,
                     timeseries_errors=false,
                     ) : appxtrue(sol, appxsol;
                                  dense_errors=false,
                                  timeseries_errors=false,
                                  )

        # Time
        # b = @benchmark $sol_call() samples=5 seconds=100
        # b = @benchmark $sol_call()
        tbest = 9999999999999999
        for i in 1:3
            t = @elapsed sol_call()
            tbest = min(tbest, t)
        end


        r = Dict(
            :final => errsol.errors[:final],
            # :l2 => errsol.errors[:l2],
            # :L2 => errsol.errors[:L2],
            :time => tbest,
            # :time => minimum(b).time / 1e9,
            # :memory => minimum(b).memory / 2^30,
            :nf => errsol.destats.nf,
            :njacs => errsol.destats.njacs,
        )
        if sol isa ProbNumDiffEq.ProbODESolution
            r[:chi2_final] = chi2(sol.pu[end], appxsol.(sol.t[end]))[1]
            # r[:final_std] = std(sol.pu[end])
        end


        push!(results, r)
    end
    return results
end


wps = Dict()



abstols = 1.0 ./ 10.0 .^ (6:12)
reltols = 1.0 ./ 10.0 .^ (3:9)


wps["Tsit5"] = MyWorkPrecision(prob1, Tsit5(), abstols ./ 10, reltols ./ 10; dense=false)
wps["RadauIIA5"] = MyWorkPrecision(prob1, RadauIIA5(), abstols ./ 10, reltols ./ 10; dense=false)
wps["Vern6"] = MyWorkPrecision(prob1, Vern6(), abstols ./ 10, reltols ./ 10; dense=false)
wps["DPRKN6"] = MyWorkPrecision(prob2, DPRKN6(), abstols, reltols; dense=false)


# order = 3
smooth = false
# for o in [5,4]
#     wps["ODE1;o=$o;EK1"] = MyWorkPrecision(
#         prob1, EK1(order=o-1, smooth=smooth, initial_derivatives=dfs1),
#         abstols, reltols; dense=smooth, save_everystep=false)
#     wps["ODE2;o=$o;EK1"] = MyWorkPrecision(
#         prob2, EK1(order=o, smooth=smooth, initial_derivatives=dfs2),
#         abstols, reltols; dense=smooth, save_everystep=false)
# end
wps["ODE1;o=4;EK0"] = MyWorkPrecision(
    prob1, EK0(order=3, smooth=smooth),
    abstols, reltols; dense=smooth, save_everystep=false)
wps["ODE2;o=4;EK0"] = MyWorkPrecision(
    prob2, EK0(order=4, smooth=smooth),
    abstols, reltols; dense=smooth, save_everystep=false)
wps["ODE1;o=5;EK1"] = MyWorkPrecision(
    prob1, EK1(order=4, smooth=smooth),
    abstols, reltols; dense=smooth, save_everystep=false)
wps["ODE2;o=5;EK1"] = MyWorkPrecision(
    prob2, EK1(order=5, smooth=smooth),
    abstols, reltols; dense=smooth, save_everystep=false)




using JLD
save(joinpath(DIR, "workprecisiondata.jld"), "wps", wps)
