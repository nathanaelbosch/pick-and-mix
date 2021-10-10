using CairoMakie
using ColorSchemes
using Plots: RGB


PN_COLORS = parse.(
    RGB, ["#107D79", "#FF9933", "#1F77B4", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"]
)

# Gruvbox
COLORS = parse.(RGB, [
    "#a54242", # red
    "#8c9440", # green
    "#de935f", # yellow
    "#5f819d", # blue
    "#85678f", # purple
    "#5e8d87", # aqua
    "#707880", # white
])

# Other gruvbox: https://camo.githubusercontent.com/410b3ab80570bcd5b470a08d84f93caa5b4962ccd994ebceeb3d1f78364c2120/687474703a2f2f692e696d6775722e636f6d2f776136363678672e706e67
GRUVBOX_DARK = parse.(RGB, [
    "#cc241d", # red
    "#98971a", # green
    "#d79921", # yellow
    "#458588", # blue
    "#b16286", # purple
    "#689d6a", # aqua
    "#d65d0e", # orange
])
# Lighter colors
GRUVBOX_LIGHT = parse.(RGB, [
    "#fb4934", # red
    "#b8bb26", # green
    "#fabd2f", # yellow
    "#83a598", # blue
    "#d3869b", # purple
    "#8ec07c", # aqua
    "#fe8019", # orange
])

set_theme!()
theme = Theme(
    fontsize=10,
    Figure=(
        figure_padding=5
    ),
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
        xticklabelsize=8,
        yticklabelsize=8,
    ),
    Legend=(
        labelsize=9,
        framevisible=false,
        patchsize=(10,10),
        padding=(0,0,0,0),
    ),
    Scatter = (
        strokewidth=0.5,
    ),
    ScatterLines = (
        linewidth=3,
        strokewidth=0.5,
    )
)
set_theme!(theme)
