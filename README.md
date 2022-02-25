# Pick-and-Mix Information Operators for Probabilistic ODE Solvers - Code

This repo contains the code which was used to compute the results of the paper **"Pick-and-Mix Information Operators for Probabilistic ODE Solvers"**.


---

__To solve differential equations in Julia with probabilistic numerical solvers, please use
[ProbNumDiffEq.jl](https://github.com/nathanaelbosch/ProbNumDiffEq.jl)!__<br />
The code in this repository is not meant to be used as generic ODE solvers, whereas
[ProbNumDiffEq.jl](https://github.com/nathanaelbosch/ProbNumDiffEq.jl)
is a Julia package under active development.
It is more stable and documented, its solvers are more efficent, and it contains more features.
The DE solvers it provides are compatible with the
[DifferentialEquations.jl](https://docs.sciml.ai/stable/)
ecosystem.

---

A __Python__ implementation of these solvers, as well as of additional probabilistic numerical methods, is maintained in [ProbNum](https://github.com/probabilistic-numerics/probnum).


## Usage
### Figure 1
Run
```
julia --project="." -e "include("experiments/0_kepler_samples/main.jl")
```
to obtain `./experiments/0_kepler_samples/figure1.pdf`:
<img alt="Figure 1" src="./experiments/0_kepler_samples/figure1.svg" width="50%">

### Experiment 1: Second-order ODEs
The experiment script is `./experiments/1_secondorder_odes/1_compute_solutions.jl`.
Then, with `./experiments/1_secondorder_odes/2_plot_solutions.jl` you obtain the plot
<img alt="Figure 2" src="./experiments/1_secondorder_odes/figure2_secondorder_workprecision.svg" width="100%">

### Experiment 2: Additional second-derivative information
The experimentc can be run with
- `./experiments/2_additional_derivatives/2.1_lotkavolterra.jl`
- `./experiments/2_additional_derivatives/2.2_vanderpol.jl`

Then, to plot run `./experiments/2_additional_derivatives/3_plot_all.jl`

<img alt="Figure 3" src="./experiments/2_additional_derivatives/figure3_additional_derivative_info.svg" width="100%">

## Reference
```
TBD
```
