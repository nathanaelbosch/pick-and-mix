using CairoMakie
import Plots: RGB
using JLD

include("../theme.jl")
DIR = @__DIR__

d = load(joinpath(DIR, "logistic_workprecisiondata.jld"))
wps = d["wps"]


fig = Figure(
    resolution=(600,170)
)


##########################################################################################
# Logistic
ax1 = fig[1, 1] = Axis(fig, yscale=log10, xscale=log10)
ax2 = fig[1, 2] = Axis(fig, yscale=log10, xscale=log10)
x, y = :nf, :final
npn_keys = ("Tsit5", "RadauIIA5")
npn_labels = ("Runge-Kutta (Tsit5)", "Implicit RK (RadauIIA5)")
npn_sclines = []
for i in (1,2)
for (j,(k,l)) in enumerate(zip(npn_keys, npn_labels))
    wp = wps[k]
    scl = scatterlines!(
        fig[1,i],
        [r[x] for r in wp[1:end-1]], [r[y] for r in wp[1:end-1]],
        color=(:gray, 0.5),
        markersize=4,
        linewidth=3,
        strokewidth=0.5,
        marker=[:utriangle, :dtriangle, :rect, :pentagon][j],
    )
    (i==1) && push!(npn_sclines, scl)
end
end

orders = [3,5]
fdbs = [0,1,2,3]
pn_sclines = []
for (i, o) in enumerate(orders),
    (j, f) in enumerate(fdbs)

    k = f == 0 ? "$o" : "$o,$f"
    @info "?" k

    wp = wps[k]

    scl = scatterlines!(
        fig[1,i],
        [r[x] for r in wp], [r[y] for r in wp],
        color=(COLORS[j], 0.5),
        markercolor=COLORS[j],
        # linestyle=[:solid, :dash, :dot][j],
        marker=[:circle, :diamond, :star5, :triangle_down][j],
        linewidth=3,
        strokewidth=0.5,
        markersize=10,
    )
    (i == 1) && push!(pn_sclines, scl)
end




leg = fig[:, end+1] = Legend(
    fig,
    [pn_sclines..., npn_sclines...],
    [
        "no ÿ information",
        "ÿ with approx. Jac. (1)",
        "ÿ with approx. Jac. (2)",
        "ÿ with exact Jac.",
        "Runge-Kutta (Tsit5)",
        "Implicit R-K (RadauIIA5)",
    ],
    patchsize=(10,10),
)




trim!(fig.layout)

# save("test.pdf", fig, pt_per_unit=1)
save(joinpath(DIR, "logistic_workprecision.pdf"), fig, pt_per_unit=1)
