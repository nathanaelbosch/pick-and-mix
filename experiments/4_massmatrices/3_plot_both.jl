using CairoMakie
using ColorSchemes
using JLD

include("../theme.jl")
DIR = @__DIR__

d = load(joinpath(DIR, "pendulum_wps.jld"))
wps_pendulum = d["wps"]
d = load(joinpath(DIR, "rober_wps.jld"))
wps_rober = d["wps"]
wps_rober["Rosenbrock23"] = wps_rober["R23"]
@info "??" keys(wps_rober) keys(wps_pendulum)

C1, C2, C3 = ColorSchemes.Set1_3
LINEALPHA = 0.6
STYLE = Dict(
    "Rosenbrock23" => (:dtriangle, :gray),
    "Rodas5" => (:utriangle, :gray),
    "Rodas4" => (:rtriangle, :gray),
    "Radau" => (:rect, :gray),
    "QNDF" => (:diamond, :gray),
    "EK1(2)" => (:star4, C1),
    "EK1(3)" => (:star5, C2),
    "EK1(5)" => (:star8, C3),
)

npn_keys, npn_labels = zip([
    ("Rosenbrock23", "Rosen-\nbrock23"),
    ("Rodas5", "Rodas5"),
    ("QNDF", "QNDF"),
]...)

pn_keys, pn_labels = zip([
    ("EK1(2)", "EK1(2)"),
    ("EK1(3)", "EK1(3)"),
    ("EK1(5)", "EK1(5)"),
]...)


fig = Figure(
    resolution=(300,200),
    figure_padding=5,
)

gl1 = fig[1, 1:2] = GridLayout()
gl2 = fig[2, 1:2] = GridLayout()
gl3 = fig[:, 3] = GridLayout()


######################################################################################
npn_sclines = Dict()
pn_sclines = Dict()
function plot_wps!(wps, axes)
    ax1, ax2 = axes
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
end


############################################################################################
p1ax1 = Axis(gl1[1,1], yscale=log10, xscale=log10,
             # xticklabelsvisible=false
             )
p1ax2 = Axis(gl1[1,2], yscale=log10, xscale=log10,
             # xticklabelsvisible=false,
             yticklabelsvisible=false)
p1ax1.ylabel = "Final error"
# p1ax1.xlabel = "Number of evaluations"
# p1ax2.xlabel = "Runtime [s]"
plot_wps!(wps_rober, (p1ax1, p1ax2))

############################################################################################
p2ax1 = Axis(gl2[1,1], yscale=log10, xscale=log10)
p2ax2 = Axis(gl2[1,2], yscale=log10, xscale=log10, yticklabelsvisible=false)
p2ax1.ylabel = "Final error"
p2ax1.xlabel = "Number of evaluations"
p2ax2.xlabel = "Runtime [s]"
plot_wps!(wps_pendulum, (p2ax1, p2ax2))

######################################################################################
leg = Legend(
    gl3[1,1],
    [[pn_sclines[k] for k in pn_keys]..., [npn_sclines[k] for k in npn_keys]...],
    [pn_labels..., npn_labels...],
    # rowgap=0,
)

######################################################################################
linkyaxes!(p1ax1, p1ax2)
linkyaxes!(p2ax1, p2ax2)
# linkxaxes!(p1ax1, p2ax1)
# linkxaxes!(p1ax2, p2ax2)

colgap!(gl1, 10)
colgap!(gl2, 10)
rowgap!(fig.layout, 10)
colgap!(fig.layout, 2, 10)

Label(gl1[1,1,TopLeft()], "A", font=noto_sans_bold, textsize=10)
Label(gl2[1,1,TopLeft()], "B", font=noto_sans_bold, textsize=10)
Label(gl1[:, :, Top()], "Robertson DAE", padding = (0, 0, 5, 0))
Label(gl2[:, :, Top()], "Pendulum DAE", padding = (0, 0, 5, 0))


# p1ax2.xticks = p2ax2.xticks = ([1e-4, 1e-2, 1e0], ["10⁻⁴", "10⁻²", "10⁰"])
p1ax1.yticks = p1ax2.yticks = ([1e-6, 1e-9, 1e-12], ["10⁻⁶", "10⁻⁹", "10⁻¹²"])
p1ax1.xticks = ([1e3, 1e4], ["10³", "10⁴"])

save(joinpath(DIR, "dae.pdf"), fig, pt_per_unit=1)
