# include(joinpath(@__DIR__), "0_kepler_samples/main.jl")
# include(joinpath(@__DIR__), "1_secondorder_odes/1_simple_example.jl")
include(joinpath(@__DIR__, "1_secondorder_odes/2.1_compute_solutions.jl"))
include(joinpath(@__DIR__, "1_secondorder_odes/2.2_plot_solutions.jl"))
include(joinpath(@__DIR__, "2_additional_derivatives/1_simple_example.jl"))
include(joinpath(@__DIR__, "2_additional_derivatives/2.1_lotkavolterra.jl"))
include(joinpath(@__DIR__, "2_additional_derivatives/2.2_vanderpol.jl"))
include(joinpath(@__DIR__, "2_additional_derivatives/3_plot_all.jl"))
include(joinpath(@__DIR__, "4_massmatrices/1.2_rober.jl"))
include(joinpath(@__DIR__, "4_massmatrices/2_plot.jl"))


include(joinpath(@__DIR__, "3_conserved_quantities/workprecision.jl"))
include(joinpath(@__DIR__, "3_conserved_quantities/wps_plot.jl"))
