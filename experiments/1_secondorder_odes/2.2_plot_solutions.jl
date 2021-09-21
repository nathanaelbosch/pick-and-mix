############################################################################################
# Plotting
using CairoMakie
using Plots: RGB
using JLD

DIR = @__DIR__

d = load(joinpath(DIR, "workprecisiondata.jld"))
wps = d["wps"]


COLORS = parse.(
    RGB, ["#107D79", "#FF9933", "#1F77B4", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"]
)


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






set_theme!()
theme = Theme(
    fontsize=10,
    Axis=(
        xgridvisible=false,
        ygridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        spinewidth=0.7,
        xtickwidth=0.7,
        ytickwidth=0.7,
        xticksize=2,
        yticksize=2,
        xticklabelsize=9,
        yticklabelsize=9,
    ),
    Legend=(
        labelsize=9,
        framevisible=false,
    )
)
set_theme!(theme)






fig = Figure(
    # resolution=(600,200)
    resolution=(600,180)
)



ax = fig[1, 1] = Axis(
    fig, title="Solution Trajectories",
    xlabel="x", ylabel="y",
    xticksvisible=false,
    yticksvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
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
           # color=COLORS[2+i],
           color=:gray,
           )
    scatter!(ax,
             [appxsol.u[end][14+i]],
             [appxsol.u[end][21+i]],
             # color=COLORS[2+i],
             markersize=5,
             strokewidth=0.5,
             color=:gray,
           )
end



npn_keys = (
    "Tsit5",
    "DPRKN6",
    "RadauIIA5",
    "Vern6",
)
npn_labels = (
    # "Runge-Kutta (Tsit5)",
    "Runge-Kutta (Vern6)",
    "R-K-Nyström (DPRKN6)",
    "Implicit R-K (RadauIIA5)",
)
npn_sclines = []
for (j,(k,l)) in enumerate(zip(npn_keys, npn_labels))
    for (i,x) in enumerate((:nf, :time))
        wp = wps[k]
        scl = scatterlines!(
            fig[1,i+1],
            [r[x] for r in wp], [r[:final] for r in wp],
            color=(:gray, 0.5),
            markersize=4,
            linewidth=3,
            strokewidth=0.5,
            marker=[:utriangle, :dtriangle, :rect, :pentagon][j],
        )
        (i == 1) && push!(npn_sclines, scl)
    end
end


pn_keys = (
    "ODE1;o=4;EK0",
    "ODE2;o=4;EK0",
    "ODE1;o=5;EK1",
    "ODE2;o=5;EK1",
)
pn_labels = (
    # "first-order ODE; IWP(2)",
    "1st order  (IWP(3), EK0)",
    "2nd order (IWP(4), EK0)",
    "1st order  (IWP(4), EK1)",
    "2nd order (IWP(5), EK1)",
)
pn_sclines = []
for (j,(k,l)) in enumerate(zip(pn_keys, pn_labels))
    for (i,x) in enumerate((:nf, :time))
        wp = wps[k]
        scl = scatterlines!(
            fig[1,i+1],
            [r[x] for r in wp], [r[:final] for r in wp],
            # color=(COLORS[j],0.5),
            color=([COLORS[1], COLORS[2], COLORS[1], COLORS[2]][j], 0.5),
            markercolor=[COLORS[1], COLORS[2], COLORS[1], COLORS[2]][j],
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
    patchsize=(10,10),
)




ax2.xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"])
# xlims!(ax2, 10^(2.8), 10^(4.05))

ax2.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
# ylims!(ax2, 1e-10, 1e0)



ax3.xticks = ([1e-2, 1e0, 1e2], ["10⁻²", "10⁰", "10²"])
# ax3.yticks = ([1e0, 1e1, 1e2], ["10⁰", "10¹", "10²"])
# ax3.yticks = ([1e0, 1e-5, 1e-10], ["10⁰", "10⁻⁵", "10⁻¹⁰"])
ax3.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
# ylims!(ax3, 1e0, 1e3)


ax3.ylabelvisible=false
ax3.yticklabelsvisible=false


save(joinpath(DIR, "fig3_secondorder_workprecision.pdf"), fig, pt_per_unit=1)
