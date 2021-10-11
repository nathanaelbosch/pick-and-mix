using CairoMakie
using JLD
using Statistics
using ColorSchemes

include("../theme.jl")
DIR = @__DIR__

d = load(joinpath(DIR, "henonheiles.jld"))
wps = d["wps"]


# COLORS = ColorSchemes.BrBG_4[[2,1,3,4]]
COLORS = ColorSchemes.BrBG_4[[3,4,2,1]]
MARKERS = Dict(
    "Tsit5" => :dtriangle,
    "Tsit5 w/ g" => :utriangle,
    "DPRKN6" => :ltriangle,
    "SymplecticEuler" => :rtriangle,
    "KahanLi8" => :rect,
    "EK1(3)" => :diamond,
    "EK1(3) w/ g" => :star4,
    "EK1(5)" => :pentagon,
    "EK1(5) w/ g" => :star5,
    "EK1(8)" => :hexagon,
    "EK1(8) w/ g" => :star8,
)


fig = Figure(
    resolution=(600,120),
    figure_padding=5,
)

ax1 = fig[1, 1] = Axis(fig, yscale=log10, xscale=log10)
ax2 = fig[1, 2] = Axis(fig, yscale=log10, xscale=log10)
ax3 = fig[1, 3] = Axis(fig, yscale=log10, xscale=log10)
axes = [ax1, ax2, ax3]

# [
# ],
npn_keys, npn_labels = zip([
    # ("Tsit5", "Runge-Kutta (Tsit5)"),
    ("Tsit5 w/ g", "RK + projection. (Tsit5)"),
    ("DPRKN6", "RK-Nystöm (DPRKN6)"),
    # ("SymplecticEuler",
    ("KahanLi8", "Symplectic (KahanLi8)"),
]...)

pn_keys, pn_labels = zip([
    ("EK1(3)", "EK1(3) (only ODE)"),
    # ("EK1(5)", "EK1(5) (only ODE)"),
    ("EK1(8)", "EK1(8) (only ODE)"),
    ("EK1(3) w/ g", "EK1(3) (with energy)"),
    # ("EK1(5) w/ g", "EK1(5) (with energy)"),
    ("EK1(8) w/ g", "EK1(8) (with energy)"),
]...)

npn_sclines = Dict()
pn_sclines = Dict()
for (r, (x,y)) in enumerate(((:nf, :final), (:time, :final), (:nf, :e_final)))

    for (i, k) in enumerate(npn_keys)
        @info "ahm" k
        wp = wps[k]
        scline = scatterlines!(
            axes[r],
            [r[x] for r in wp],
            [max(abs(r[y]), 1e-50) for r in wp],
            marker=MARKERS[k],
            markersize=4,
            color=(:gray, 0.5),
            markercolor=:gray,
            linewidth=3,
            strokewidth=0.5,
        )
        npn_sclines[k] = scline
    end

    # pn_keys = ("EK1(5)",)
    for (i, k) in enumerate(pn_keys)
        @info "ahm" k
        wp = wps[k]
        scline = scatterlines!(
            axes[r],
            [r[x] for r in wp],
            [max(abs(r[y]), 1e-30) for r in wp],
            color=(COLORS[i], 0.5),
            markercolor=COLORS[i],
            marker=MARKERS[k],
            strokewidth=0.5,
            linewidth=3,
            # markersize=10,
            markersize=8,
        )
        pn_sclines[k] = scline
    end

end

leg = Legend(
    fig[:, end+1],
    [[pn_sclines[k] for k in pn_keys]..., [npn_sclines[k] for k in npn_keys]...],
    [pn_labels..., npn_labels...],
    labelsize=8,
    # tellheight=true,
    rowgap=2,
)

ax3.xlabel="Number of evaluations"
ax1.xlabel="Number of evaluations"
ax2.xlabel="Runtime [s]"
# ax2.xlabel="Runtime"
ax3.ylabel="Change in energy"
ax1.ylabel="Final error"

ax3.xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"])
ax1.xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"])
# hidexdecorations!(ax1, grid=false, minorgrid=false)

# ylims!(ax2, 1e-15, 10^(-0.5))
ylims!(ax3, 10^-18, 2e-0)
ax3.yticks = ([1e0, 1e-5, 1e-10, 1e-15, 1e-20], ["10⁰", "10⁻⁵", "10⁻¹⁰", "10⁻¹⁵", "10⁻²⁰"])

# ylims!(ax1, 1e-13, 10^(-0.5))
ax1.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
ax2.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
# ax1.yticks = ([1e0, 1e-5, 1e-10, 1e-20], ["10⁰", "10⁻⁵", "10⁻¹⁰", "10⁻²⁰"])
# ax2.yticks = ([1e0, 1e-5, 1e-10, 1e-20], ["10⁰", "10⁻⁵", "10⁻¹⁰", "10⁻²⁰"])


ax2.xticks = ([1e-4, 1e-2, 1e0], ["10⁻⁴", "10⁻²", "10⁰"])

ax2.yticklabelsvisible=false
ax2.ylabelvisible=false

colgap!(fig.layout, 1, 10)
colgap!(fig.layout, 2, 15)
colgap!(fig.layout, 3, 15)
trim!(fig.layout)

save(joinpath(DIR, "henonheiles_workprecision.pdf"), fig, pt_per_unit=1)
