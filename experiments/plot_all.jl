#=
Run all the plotting scripts

Useful after tweaking the plotting theme.
=#

include(joinpath(@__DIR__, "0_kepler_samples/main.jl"))

include(joinpath(@__DIR__, "1_secondorder_odes/2.2_plot_solutions.jl"))

include(joinpath(@__DIR__, "2_additional_derivatives/3_plot_all.jl"))

include(joinpath(@__DIR__, "3_conserved_quantities/wps_plot.jl"))
include(joinpath(@__DIR__, "3_conserved_quantities/longtermplot.jl"))
lnclude(joinpath(@__DIR__, "3_conserved_quantities/inflated_kepler_samples.jl"))

include(joinpath(@__DIR__, "4_massmatrices/2_plot.jl"))
