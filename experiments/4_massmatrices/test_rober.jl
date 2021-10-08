using DifferentialEquations
using ProbNumDiffEq
using Plots

using LinearAlgebra
using Statistics
using DiffEqDevTools


function rober_mm(du,u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    du[1] = -k₁*y₁ + k₃*y₂*y₃
    du[2] =  k₁*y₁ - k₃*y₂*y₃ - k₂*y₂^2
    du[3] =  y₁ + y₂ + y₃ - 1
    nothing
end

M = [1. 0  0
     0  1. 0
     0  0  0]
f_mm = ODEFunction(rober_mm, mass_matrix=M)
u0 = [1.0,0.0,0.0]
tspan = (0.0, 1e2)
# tspan = (0.0, 1e1)
p = (0.04,3e7,1e4)
prob_mm = ODEProblem(f_mm, u0, tspan, p)
# bigprob_mm = remake(prob_mm, u0=big.(prob_mm.u0), tspan=big.(prob_mm.tspan),
# p=big.(prob_mm.p))

function rober_ode(du,u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    du[1] = -k₁*y₁+k₃*y₂*y₃
    du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
    du[3] =  k₂*y₂^2
    nothing
end
prob_ode = ODEProblem(rober_ode, u0, tspan, p)


solve_ref(prob) = solve(remake(prob, u0=big.(prob.u0), tspan=big.(prob.tspan), p=big.(prob.p)),
                        RadauIIA5(), reltol=1e-18, abstol=1e-18)
appxsol_mm = appxsol = solve_ref(prob_mm);
appxsol_ode = solve_ref(prob_ode);


Plots.plot(appxsol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))

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
                     # dense_errors=false
                     ) : appxtrue(sol, appxsol)

        # Time
        tbest = 100
        for i in 1:3
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
        )

        push!(results, r)
    end
    return results
end


# _prob = prob_mm
# _prob = prob_ode
wps = Dict()

for (_l, _prob, _appxsol) in (
    ("DAE", prob_mm, appxsol_mm),
    ("ODE", prob_ode, appxsol_ode)
    )


    abstols2 = abstols1 = 1.0 ./ 10.0 .^ (5:9)
    reltols2 = reltols1 = 1.0 ./ 10.0 .^ (2:6)
    # wps["Rodas5 ($_l)"] = MyWorkPrecision(_prob, Rodas5(), abstols1, reltols1; appxsol=_appxsol)
    wps["R23 ($_l)"] = MyWorkPrecision(_prob, Rosenbrock23(), abstols1, reltols1; appxsol=_appxsol)
    wps["Radau ($_l)"] = MyWorkPrecision(_prob, RadauIIA5(), abstols1, reltols1; appxsol=_appxsol)
    # wps["QNDF ($_l)"] = MyWorkPrecision(_prob, QNDF(), abstols1, reltols1; appxsol=_appxsol)


    # abstols2 = 1.0 ./ 10.0 .^ (4:9)
    # reltols2 = 1.0 ./ 10.0 .^ (1:6)
    for o in (3, 5)
        @info "EK1($o)"
        wps["EK1($o) ($_l)"] = MyWorkPrecision(_prob, EK1(order=o,
                                                          # smooth=false
                                                          ), abstols2, reltols2;
                                               appxsol=_appxsol,
                                               # dense=false, save_everystep=false
                                               )
    end
end




# x, y = :time, :final
x, y = :nf, :final
p = Plots.plot(xscale=:log10, yscale=:log10, xlabel="#evaluations", ylabel="#Final error",
    # ylims=(1e-20, 1e0),
         legend=:outerright,
)
for (i, label) in enumerate(sort(collect(keys(wps))))
    wp = wps[label]
    # color = occursin("EK1", label) ? :auto : :gray
    color = (i+1) ÷ 2
    linestyle = (i-1)%2 == 0 ? :dash : :solid
    linewidth = (i-1)%2 == 0 ? 2 : 4
    linealpha = (i-1)%2 == 0 ? 0.6 : 1
    Plots.plot!(p, [r[x] for r in wp], [r[y] for r in wp],
          color=color,
          linestyle=linestyle,
          marker=:auto, markerstrokewidth=0.5,
          linewidth=linewidth, linealpha=linealpha,
          label=label)
end
p





Plots.plot!(p1, legend=false)
Plots.plot!(p2, xlabel="runtime")
Plots.plot!(p1, p2, title="")
savefig("rober_mm_dae.pdf")














######################################################
# Bigproblem test - Why do things not get more exact??
# bp = remake(prob_mm, u0=big.(prob_mm.u0), tspan=big.(prob_mm.tspan), p=big.(prob_mm.p))
# solve(bp, EK1(order=3))
