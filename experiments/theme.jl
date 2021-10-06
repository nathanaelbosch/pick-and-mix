using CairoMakie

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
        patchsize=(13,13),
    )
)
set_theme!(theme)
