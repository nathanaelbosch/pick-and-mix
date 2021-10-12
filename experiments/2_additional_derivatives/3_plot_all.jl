using CairoMakie
using ColorSchemes
using JLD

include("1_simple_example.jl")

include("../theme.jl")
DIR = @__DIR__

d = load(joinpath(DIR, "vdp_wps.jld"))
wps_vdp = d["wps"]
d = load(joinpath(DIR, "lv_wps.jld"))
wps_lv = d["wps"]

# C1, C2, C3, C4 = ColorSchemes.seaborn_muted[[4,3,4,3]]
# C1, C2, C3, C4 = ColorSchemes.Accent_3[[2,3,2,3]]
# C1, C2, C3, C4 = ColorSchemes.BrBG_4[[2,3,1,4]]
# ColorSchemes.PRGn_4
# C1, C2, C3, C4 = ColorSchemes.Paired_4[[1,3,2,4]]
# C1, C2, C3, C4 = ColorSchemes.Paired_4[[1,2,3,4]]
# C1, C2, C3, C4 = ColorSchemes.BrBG_4[[2,3,1,4]]
C1, C2, C3, C4 = ColorSchemes.PuOr_4[[2,3,1,4]]
# C1, C2, C3, C4 = ColorSchemes.Dark2_4
# C1, C2, C3, C4 = ColorSchemes.Set1_4
# C1, C2, C3, C4 = ColorSchemes.Set2_4
# C1, C2, C3, C4 = ColorSchemes.Set3_7[4:7]
# C1, C2 = GRUVBOX_DARK[1], GRUVBOX_LIGHT[2]
LINEALPHA = 0.6
MARKERS = Dict(
    "Tsit5" => :dtriangle,
    "RK45-scipy" => :utriangle,
    "RadauIIA5" => :rect,
    "Radau-scipy" => :pentagon,
    "LSODA-scipy" => :hexagon,
    "3" => :diamond,
    "3,3" => :star4,
    "5" => :pentagon,
    "5,3" => :star5
)



fig = Figure(
    resolution=(600,250),
    figure_padding = 5,
)
gl1 = fig[1:4, 1] = GridLayout()
gl2 = fig[1:2, 2:3] = GridLayout()
g3 = fig[3:4, 2:3] = GridLayout()
gl = fig[:, 4] = GridLayout()


######################################################################################
# Left: simple plot
######################################################################################
p1ax1 = Axis(gl1[1,1], ylabel=L"Y^{(0)}", xticksvisible=false, xticklabelsvisible=false, xgridvisible=false, ygridvisible=false)
p1ax2 = Axis(gl1[2,1], ylabel=L"Y^{(0)} - y", xticksvisible=false, xticklabelsvisible=false, xgridvisible=false, ygridvisible=false)
p1ax3 = Axis(gl1[3,1], ylabel=L"Y^{(1)} - \dot{y}", xticksvisible=false, xticklabelsvisible=false, xgridvisible=false, ygridvisible=false)
p1ax4 = Axis(gl1[4,1], ylabel=L"Y^{(2)} - \ddot{y}", xlabel=L"t", xgridvisible=false, ygridvisible=false)

lines!(p1ax1, 0:0.1:3, true_solution.(0:0.1:3), linestyle=:dash, color=:black)
hlines!(p1ax2, [0], linestyle=:dash, color=:black)
hlines!(p1ax3, [0], linestyle=:dash, color=:black)
hlines!(p1ax4, [0], linestyle=:dash, color=:black)

plot_sol(sol1, 1;
         color=C1,
         marker=MARKERS["3"],
         axes=(p1ax1, p1ax2, p1ax3, p1ax4))
plot_sol(sol2, 2;
         color=C2,
         marker=MARKERS["5"],
         axes=(p1ax1, p1ax2, p1ax3, p1ax4))

p1ax1.xticks = [0,3]
p1ax2.xticks = [0,3]
p1ax3.xticks = [0,3]
p1ax4.xticks = [0,3]

ylims!(p1ax2, -0.0125, 0.0125)
p1ax2.yticks = [-0.01, 0.01]
ylims!(p1ax3, -0.025, 0.025)
p1ax3.yticks = [-0.02, 0.02]
ylims!(p1ax4, -0.25, 0.25)
p1ax4.yticks = [-0.2, 0.2]
ylims!(p1ax1, -0.2, 1.2)
p1ax1.yticks = [0,1]

rowgap!(gl1, 10)
# valign = :bottom, padding = (0, 0, 5, 0))

######################################################################################
# Top: Lotka-Volterra work-precision
######################################################################################
p2ax1 = Axis(gl2[1,1], xscale=log10, yscale=log10,
             xlabel = "Number of evaluations",
             ylabel = "Final error",
             )
p2ax2 = Axis(gl2[1,2], xscale=log10, yscale=log10,
             xlabel = "Runtime [s]",
             ylabel = "Final error",
             )

npn_keys = ("Tsit5", "RadauIIA5")
pn_labels = pn_keys = ("3", "3,3", "5", "5,3")

npn_sclines = Dict()
pn_sclines = Dict()
for (i, x, ax) in ((1, :nf, p2ax1), (2,:time, p2ax2))
    y = :final

    for (j,k) in enumerate(npn_keys)

        wp = wps_lv[k]

        (i == 1 && occursin("scipy", k)) && continue

        scl = scatterlines!(
            ax,
            [r[x] for r in wp[1:end-1]], [r[y] for r in wp[1:end-1]],
            color=(:gray, 0.5),
            markersize=4,
            markercolor=:gray,
            marker=MARKERS[k],
        )

        npn_sclines[k] = scl
    end

    for (j, k) in enumerate(pn_keys)

        wp = wps_lv[k]

        scl = scatterlines!(
            ax,
            [r[x] for r in wp], [r[y] for r in wp],
            color=([C1, C2, C3, C4][j], LINEALPHA),
            markercolor=[C1, C2, C3, C4][j],
            # linestyle=[:solid, :dash, :dot][j],
            marker=MARKERS[k],
            markersize=10,
        )
    end
end


p2ax2.ylabelvisible=false
p2ax2.yticklabelsvisible=false
linkyaxes!(p2ax1, p2ax2)
colgap!(gl2, 10)


p2ax1.xticks = ([1e2, 1e3], ["10²", "10³"])
p2ax2.yticks = p2ax1.yticks = ([1e-5, 1e-10], ["10⁻⁵", "10⁻¹⁰"])
p2ax2.xticks = ([1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
                ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"])


######################################################################################
# Bottom: Van-der-Pol work-precision
######################################################################################
p3ax1 = Axis(g3[1,1], xscale=log10, yscale=log10,
             xlabel = "Number of evaluations",
             ylabel = "Final error",
             )
p3ax2 = Axis(g3[1,2], xscale=log10, yscale=log10,
             xlabel = "Runtime [s]",
             ylabel = "Final error",
             )


npn_labels = npn_keys = ("RadauIIA5",)
pn_labels = pn_keys = ("3", "3,3", "5", "5,3")

for (i, x, ax) in ((1, :nf, p3ax1), (2,:time, p3ax2))
    y = :final

    for (j,k) in enumerate(npn_keys)

        wp = wps_vdp[k]

        (i == 1 && occursin("scipy", k)) && continue

        scl = scatterlines!(
            ax,
            [r[x] for r in wp[1:end-1]], [r[y] for r in wp[1:end-1]],
            color=(:gray, 0.5),
            markersize=4,
            markercolor=:gray,
            marker=MARKERS[k],
        )

        npn_sclines[k] = scl
    end

    for (j,k) in enumerate(pn_keys)

        wp = wps_vdp[k]

        scl = scatterlines!(
            ax,
            [r[x] for r in wp], [r[y] for r in wp],
            color=([C1, C2, C3, C4][j], LINEALPHA),
            markercolor=[C1, C2, C3, C4][j],
            # linestyle=[:solid, :dash, :dot][j],
            marker=MARKERS[k],
            markersize=10,
        )
        pn_sclines[k] = scl
    end
end

p3ax2.ylabelvisible=false
p3ax2.yticklabelsvisible=false
linkyaxes!(p3ax1, p3ax2)
colgap!(g3, 10)


p3ax1.xticks = ([1e4, 1e5], ["10⁴", "10⁵"])
p3ax1.yticks = p3ax2.yticks = ([1e-2, 1e-6, 1e-10], ["10⁻²", "10⁻⁶", "10⁻¹⁰"])
p3ax2.xticks = ([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1],
                ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"])


######################################################################################
# Right: Legend
######################################################################################

leg = Legend(
    gl[1,1],
    [
        pn_sclines["3"],
        pn_sclines["5"],
        pn_sclines["3,3"],
        pn_sclines["5,3"],
        npn_sclines["Tsit5"],
        # npn_sclines["RK45-scipy"],
        npn_sclines["RadauIIA5"],
        # npn_sclines["Radau-scipy"],
        # npn_sclines["LSODA-scipy"],
    ],
    [
        "EK1(3)̈",
        "EK1(5)",
        "EK1(3), augmented",
        "EK1(5), augmented",
        "Runge-Kutta (Tsit)",
        "Implicit RK (RadauIIA5)",
    ],
    # orientation=:horizontal,
    # tellwidth=false,
    # tellheight=true,
)


######################################################################################
# General tweaks
######################################################################################
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
Label(gl1[1,1,TopLeft()], "A", font=noto_sans_bold, textsize=10)
Label(gl2[1,1,TopLeft()], "B", font=noto_sans_bold, textsize=10)
Label(g3[1,1,TopLeft()], "C", font=noto_sans_bold, textsize=10)
Label(gl1[:, 1, Top()], "Logistic Equation",
      padding = (0, 0, 5, 0))
Label(gl2[:, :, Top()], "Lotka-Volterra",
      padding = (0, 0, 5, 0))
Label(g3[:, :, Top()], "Van-der-Pol",
      padding = (0, 0, 5, 0))

colgap!(fig.layout, 3, 15)
trim!(fig.layout)

save(joinpath(DIR, "fdb_plot.pdf"), fig, pt_per_unit=1)
