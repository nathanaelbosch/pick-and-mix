using CairoMakie
using JLD
using Statistics

include("../theme.jl")
DIR = @__DIR__



fig = Figure(
    resolution=(600,180)
)

ax1 = fig[1, 1] = Axis(fig, yscale=log10, xscale=log10)
ax1 = fig[1, 2] = Axis(fig, yscale=log10, xscale=log10)



npn_keys = (
    # "Tsit5",
    "RadauIIA5",
    # "RK45 (SciPy)",
    "LSODA (SciPy)",
)
npn_labels = npn_keys
npn_sclines = []
for (j,(k,l)) in enumerate(zip(npn_keys, npn_labels))
    for (i,x) in enumerate((:nf, :time))
        wp = wps[k]
        # @info "??" k wp isnothing(wp[1][x])
        if isnothing(wp[1][x])
            continue
        end
        scl = scatterlines!(
            fig[1,i],
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


pn_keys = (
    "EK1(3)",
    "EK1FDB(3)",
    "EK1(5)",
    "EK1FDB(5)",
)
pn_labels = (
    "EK1(3)",
    "EK1FDB(3)",
    "EK1(5)",
    "EK1FDB(5)",
)
pn_sclines = []
for (j,(k,l)) in enumerate(zip(pn_keys, pn_labels))
    for (i,x) in enumerate((:nf, :time))
        wp = wps[k]
        scl = scatterlines!(
            fig[1,i],
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



# leg = fig[:, end+1] = Legend(
#     fig,
#     [pn_sclines..., npn_sclines...],
#     [pn_labels..., npn_labels...],
#     # orientation=:horizontal,
#     # tellwidth=false,
#     # tellheight=true,
# )


save(joinpath(DIR, "fdb_workprecision.pdf"), fig, pt_per_unit=1)
