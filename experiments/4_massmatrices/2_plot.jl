using CairoMakie
using ColorSchemes
using JLD

include("../theme.jl")
DIR = @__DIR__

d = load(joinpath(DIR, "pendulum_wps.jld"))
wps_pendulum = d["wps"]
d = load(joinpath(DIR, "rober_wps.jld"))
wps_rober = d["wps"]


C1, C2, C3 = ColorSchemes.Accent_3
C1, C2, C3 = ColorSchemes.Dark2_3
# C1, C2, C3 = ColorSchemes.Set1_3
# C1, C2, C3 = ColorSchemes.Set2_3
# C1, C2, C3 = ColorSchemes.seaborn_pastel[1:3]
LINEALPHA = 0.6
STYLE = Dict(
    "R23 (DAE)" => (:dtriangle, :gray),
    "Rodas5 (DAE)" => (:utriangle, :gray),
    "Rodas4 (DAE)" => (:rtriangle, :gray),
    "Radau (DAE)" => (:rect, :gray),
    "QNDF (DAE)" => (:diamond, :gray),
    "EK1(2) (DAE)" => (:star4, C1),
    "EK1(3) (DAE)" => (:star5, C2),
    "EK1(5) (DAE)" => (:star8, C3),
    "EK1(2) (ODE + cb)" => (:diamond, C1),
    "EK1(3) (ODE + cb)" => (:diamond, C2),
    "EK1(5) (ODE + cb)" => (:diamond, C3),
)


npn_keys, npn_labels = zip([
    ("Rodas5 (DAE)", "Rodas5"),
    ("QNDF (DAE)", "QNDF"),
]...)

pn_keys, pn_labels = zip([
    ("EK1(2) (DAE)", "EK1(2)"),
    ("EK1(3) (DAE)", "EK1(3)"),
    ("EK1(5) (DAE)", "EK1(5)"),
    # ("EK1(2) (ODE + cb)", "EK1(2) ODE + cb"),
    # ("EK1(3) (ODE + cb)", "EK1(3) ODE + cb"),
    # ("EK1(5) (ODE + cb)", "EK1(5) ODE + cb"),
]...)


fig = Figure(
    resolution=(300,100),
    figure_padding=5,
)
# p0ax1 = fig[1, 1] = Axis(fig, xscale=log10, xticklabelsvisible=false)
# p0ax2 = fig[2, 1] = Axis(fig, xscale=log10, xticklabelsvisible=false)
# p0ax3 = fig[3, 1] = Axis(fig, xscale=log10, xlabel=L"t")
ax1 = fig[:, 2] = Axis(fig, yscale=log10, xscale=log10)
ax2 = fig[:, 3] = Axis(fig, yscale=log10, xscale=log10)
ax1.ylabel = "Final error"
ax1.xlabel = "Number of evaluations"
ax2.xlabel = "Runtime [s]"

wps = wps_rober

######################################################################################
npn_sclines = Dict()
pn_sclines = Dict()
for (i, x, ax) in ((1, :nf, ax1), (2,:time, ax2))
    y = :final

    for (j,k) in enumerate(npn_keys)
        (i == 1 && occursin("scipy", k)) && continue
        wp = wps[k]
        m, c = STYLE[k]
        scl = scatterlines!(
            ax,
            [r[x] for r in wp[1:end-1]], [r[y] for r in wp[1:end-1]],
            color=(:gray, 0.5),
            markersize=4,
            markercolor=:gray,
            marker=m,
        )
        npn_sclines[k] = scl
    end

    for (j,k) in enumerate(pn_keys)
        wp = wps[k]
        m, c = STYLE[k]
        scl = scatterlines!(
            ax,
            [r[x] for r in wp], [r[y] for r in wp],
            color=(c, LINEALPHA),
            markercolor=c,
            # linestyle=[:solid, :dash, :dot][j],
            marker=m,
            markersize=10,
        )
        pn_sclines[k] = scl
    end
end


######################################################################################
leg = Legend(
    fig[:,4],
    [[pn_sclines[k] for k in pn_keys]..., [npn_sclines[k] for k in npn_keys]...],
    [pn_labels..., npn_labels...],
#     # orientation=:horizontal,
#     # framevisible=true,
    rowgap=0,
)


# ax1.xticks = ([1e3, 1e4, 1e5], ["10³", "10⁴", "10⁵"])
# ax1.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
# ax2.xticks = ([1e-2, 1e0, 1e2], ["10⁻²", "10⁰", "10²"])
# ax2.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
ax1.xticks = ([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4],
              ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴"])
# ax1.xticks = ([1e2, 1e4], ["10²", "10⁴"])
# xlims!(ax1, 8e1, 3e4)
ax2.xticks = ([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1],
              ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"])
# ax2.xticks = ([1e-4, 1e-1], ["10⁻⁴", "10⁻¹"])
# xlims!(ax2, 8e-5, 1.2e-1)
ax2.ylabelvisible=false
ax2.yticklabelsvisible=false


colgap!(fig.layout, 10)
trim!(fig.layout)


save(joinpath(DIR, "dae.pdf"), fig, pt_per_unit=1)
