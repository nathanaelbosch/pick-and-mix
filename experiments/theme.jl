using CairoMakie
using Plots: RGB


# COLORS = parse.(
#     RGB, ["#107D79", "#FF9933", "#1F77B4", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"]
# )

# Gruvbox
# COLORS = parse.(RGB, ["#a54242", "#8c9440", "#de935f", "#5f819d", "#85678f", "#5e8d87", "#707880"])
COLORS = parse.(RGB, [
    "#5f819d",
    "#de935f",
    "#85678f",
    "#5e8d87",
    "#707880",
])

set_theme!()
theme = Theme(
    fontsize=10,
    Axis=(
        # xgridvisible=false,
        # ygridvisible=false,
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
        patchsize=(13,13),
    )
)
set_theme!(theme)
