using DifferentialEquations
using ProbNumDiffEq
using Plots
using LinearAlgebra
using Statistics
using DiffEqDevTools
using ModelingToolkit


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
solve(pendulum_prob,Rodas4())




traced_sys = modelingtoolkitize(pendulum_prob)
pendulum_sys = structural_simplify(dae_index_lowering(traced_sys))
prob = ODEProblem(pendulum_sys, Pair[], tspan)
bigprob = remake(prob, u0=big.(prob.u0), tspan=big.(prob.tspan), p=big.(prob.p))
sol = solve(prob, Rodas4())
# appxsol = solve(bigprob, RadauIIA5(), abstol=1e-14, reltol=1e-14)
appxsol = solve(bigprob, Rodas4(), abstol=1e-16, reltol=1e-16)

Plots.plot(sol, vars=states(traced_sys))







###########################################################################################
# Work-Precision diagrams

function chi2(gaussian_estimate, actual_value)
    μ, Σ = gaussian_estimate
    diff = μ - actual_value
    R = qr(Σ.squareroot').R
    @assert R'R ≈ Matrix(Σ)
    chi2_pinv = diff' * pinv(Matrix(Σ)) * diff
    return chi2_pinv
end
function MyWorkPrecision(prob, alg, abstols, reltols, args...; appxsol, kwargs...)
    @info "$alg"
    results = []
    for (atol, rtol) in zip(abstols, reltols)
        @info "info" atol rtol

        sol_call() = solve(prob, alg, abstol=atol, reltol=rtol, args...; kwargs...)

        sol = sol_call()
        errsol = sol isa ProbNumDiffEq.ProbODESolution ?
            appxtrue(sol, TestSolution(appxsol);
                     # dense_errors=false,
                     # timeseries_errors=false,
                     ) : appxtrue(sol, appxsol)

        # Time
        tbest = 100
        for i in 1:3
            t = @elapsed sol_call()
            tbest = min(tbest, t)
        end

        r = Dict(
            :final => errsol.errors[:final],
            # :l2 => errsol.errors[:l2],
            # :L2 => errsol.errors[:L2],
            :time => tbest,
            :nf => errsol.destats.nf,
            :njacs => errsol.destats.njacs,
            # :chi2_final => chi2(sol.pu[end], appxsol.(sol.t[end])),
        )

        push!(results, r)
    end
    return results
end


# _prob = bigprob
_prob = prob
wps = Dict()

abstols = 1.0 ./ 10.0 .^ (5:12)
reltols = 1.0 ./ 10.0 .^ (2:9)
wps["Rodas5"] = MyWorkPrecision(_prob, Rodas5(), abstols, reltols; appxsol=appxsol)
# wps["Rosenbrock23"] = MyWorkPrecision(_prob, Rosenbrock23(), abstols1, reltols1; appxsol=appxsol)
# wps["RadauIIA5"] = MyWorkPrecision(_prob, RadauIIA5(), abstols1, reltols1; appxsol=appxsol)
# wps["QNDF"] = MyWorkPrecision(_prob, QNDF(), abstols1, reltols1; appxsol=appxsol)

for (i, o) in enumerate((3, 5, 8))
    wps["EK1($o)"] = MyWorkPrecision(
        _prob, EK1(order=o, smooth=false),
        abstols, reltols; appxsol=appxsol, dense=false, save_everystep=false)
    # wps["EK1($o)-index4"] = MyWorkPrecision(
    #     pendulum_prob, EK1(order=o, smooth=false),
    #     abstols, reltols; appxsol=appxsol, dense=false, save_everystep=false)
end



x, y = :nf, :final
# x, y = :nf, :l2
p = Plots.plot(xscale=:log10, yscale=:log10, xlabel="#evaluations", ylabel="#Final error",
    # ylims=(1e-20, 1e0),
         legend=:outerright,
)
for (i, label) in enumerate(sort(collect(keys(wps))))
    wp = wps[label]
    # color = occursin("EK1", label) ? :auto : :gray
    # color = (i+1) ÷ 2
    # linestyle = (i-1)%2 == 0 ? :dash : :solid
    Plots.plot!(p, [r[x] for r in wp], [r[y] for r in wp],
          # color=color,
          # linestyle=linestyle,
          marker=:auto, markerstrokewidth=0.5,
          linewidth=4, linealpha=0.6,
          label=label)
end
title!(p, "Index-1 Pendulum DAE (MM)")
p


# Plots.plot!(p1, legend=false)
# Plots.plot!(p2, xlabel="runtime")
# Plots.plot!(p1, p2, title=false)
# savefig("pendulum_mm_dae.pdf")
