############################################################################################
# Plotting
using CairoMakie
using JLD
using Statistics

include("../theme.jl")
DIR = @__DIR__

d = load(joinpath(DIR, "workprecisiondata.jld"))
wps = d["wps"]

ODE1_o4_times = [w[:time] for w in wps["ODE1;o=3;EK0"]]
ODE2_o4_times = [w[:time] for w in wps["ODE2;o=3;EK0"]]
ODE1_o5_times = [w[:time] for w in wps["ODE1;o=5;EK1"]]
ODE2_o5_times = [w[:time] for w in wps["ODE2;o=5;EK1"]]

@info "runtimes and improvements" ODE1_o4_times ./ ODE2_o4_times mean(ODE1_o4_times ./ ODE2_o4_times) ODE1_o5_times ./ ODE2_o5_times mean(ODE1_o5_times ./ ODE2_o5_times)

using OrdinaryDiffEq
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






fig = Figure(
    resolution=(600,150),
    figure_padding = 1,
)

ax = fig[1, 1] = Axis(
    fig, title="Solution Trajectories",
    xlabel="x", ylabel="y",
    xgridvisible=false,
    ygridvisible=false,
    xticksvisible=false,
    yticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    aspect = AxisAspect(1),
)
ax2 = fig[1, 2] = Axis(fig, yscale=log10, xscale=log10)
ax3 = fig[1, 3] = Axis(fig, yscale=log10, xscale=log10)
ax2.xlabel = "Number of evaluations"
ax2.ylabel = "Final error"
ax3.xlabel = "Runtime [s]"
ax3.ylabel = "Final error"





for i in 1:7
    lines!(ax,
           [u[14+i] for u in appxsol.u],
           [u[21+i] for u in appxsol.u],
           color=:gray,
           )
    scatter!(ax,
             [appxsol.u[end][14+i]],
             [appxsol.u[end][21+i]],
             markersize=5,
             strokewidth=0.5,
             color=:gray,
           )
end



npn_keys = (
    # "Tsit5",
    "Vern6",
    "DPRKN6",
    "RadauIIA5",
    # "RK45-SciPy",
    # "LSODA-SciPy",
)
npn_labels = (
    # "Runge-Kutta (Tsit5)",
    "Runge-Kutta (Vern6)",
    "R-K-Nyström (DPRKN6)",
    "Implicit R-K (RadauIIA5)",
    # "Runge-Kutta (SciPy)",
    # "LSODA (SciPy)",
)
npn_sclines = []
for (j,(k,l)) in enumerate(zip(npn_keys, npn_labels))
    for (i,x) in enumerate((:nf, :time))
        wp = wps[k]
        @info "??" k wp isnothing(wp[1][x])
        if isnothing(wp[1][x])
            continue
        end
        scl = scatterlines!(
            fig[1,i+1],
            [r[x] for r in wp], [r[:final] for r in wp],
            color=(:gray, 0.5),
            markercolor=:gray,
            markersize=5,
            linewidth=3,
            strokewidth=0.5,
            marker=[:utriangle, :dtriangle, :rect, :rtriangle, :pentagon][j],
        )
        (i == 2) && push!(npn_sclines, scl)
    end
end


pn_labels =pn_keys = (
    "ODE1;o=3;EK0",
    "ODE2;o=3;EK0",
    "ODE1;o=5;EK1",
    "ODE2;o=5;EK1",
)
pn_labels = (
    # "first-order ODE; IWP(2)",
    "1st order  (EK0, IWP(2))",
    "2nd order (EK0, IWP(3))",
    "1st order  (EK1, IWP(4))",
    "2nd order (EK1, IWP(5))",
)
pn_sclines = []
for (j,(k,l)) in enumerate(zip(pn_keys, pn_labels))
    for (i,x) in enumerate((:nf, :time))
        wp = wps[k]
        scl = scatterlines!(
            fig[1,i+1],
            [r[x] for r in wp], [r[:final] for r in wp],
            # color=(COLORS[j],0.5),
            color=([COLORS[4], COLORS[3], COLORS[4], COLORS[3]][j], 0.5),
            markercolor=[COLORS[4], COLORS[3], COLORS[4], COLORS[3]][j],
            # marker=[:pentagon, :diamond, :star5, :star4][2i+j-2],
            marker=[:pentagon, :star5, :diamond, :star4][j],
            strokewidth=0.5, linewidth=3, markersize=10,
        )
        (i == 1) && push!(pn_sclines, scl)
    end
end



leg = fig[:, end+1] = Legend(
    fig,
    [pn_sclines..., npn_sclines...],
    [pn_labels..., npn_labels...],
    # orientation=:horizontal,
    # tellwidth=false,
    # tellheight=true,
)




ax2.xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"])
# xlims!(ax2, 10^(2.8), 10^(4.05))

ax2.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
# ylims!(ax2, 1e-10, 1e0)



ax3.xticks = ([1e-2, 1e0, 1e2], ["10⁻²", "10⁰", "10²"])
ax3.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
# ylims!(ax3, 1e0, 1e3)


ax3.ylabelvisible=false
ax3.yticklabelsvisible=false


save(joinpath(DIR, "fig3_secondorder_workprecision.pdf"), fig, pt_per_unit=1)
