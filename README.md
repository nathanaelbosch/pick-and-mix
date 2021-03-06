# Pick-and-Mix Information Operators for Probabilistic ODE Solvers - Code

This repo contains the code which was used to compute the results of the paper [**"Pick-and-Mix Information Operators for Probabilistic ODE Solvers"**](https://proceedings.mlr.press/v151/bosch22a.html), accepted at AISTATS 2022.


<p align="center">
<img alt="Figure 1" src="./experiments/0_kepler_samples/figure1.svg" width="80%">
</p>


---

__To solve differential equations in Julia with probabilistic numerical solvers, have a look at our package
[ProbNumDiffEq.jl](https://github.com/nathanaelbosch/ProbNumDiffEq.jl).__
Much of the new functionality presented in this paper is already available in ProbNumDiffEq.jl and will be covered in the documentation.

---

A __Python__ implementation of probabilistic numerical ODE solvers, as well as many other probabilistic numerical methods, is maintained in __[ProbNum](https://github.com/probabilistic-numerics/probnum)__ (`pip install probnum`).

---


## Running the code to recreate our experiments
First open `julia`, activate the local environment, and instantiate it to install all the packages:
```
julia> ]
(v1.7) pkg> activate .
(v1.7) pkg> instantiate
```
and you can quit the `pkg` environment by hitting backspace.

To run a julia script from the Julia REPL, do
```
julia> include("myscript.jl")
```

### Generating figure 1
Run the script `./experiments/0_kepler_samples/main.jl`
to obtain `./experiments/0_kepler_samples/figure1.pdf` (shown above).


### Case study 1: Second-order ODEs
The experiment script is `./experiments/1_secondorder_odes/1_compute_solutions.jl`.
Then, with `./experiments/1_secondorder_odes/2_plot_solutions.jl` you obtain the plot

<p align="center">
<img alt="Figure 2" src="./experiments/1_secondorder_odes/figure2_secondorder_workprecision.svg" width="100%">
</p>


### Case study 2: Additional second-derivative information
The experimentc can be run with
- `./experiments/2_additional_derivatives/2.1_lotkavolterra.jl`
- `./experiments/2_additional_derivatives/2.2_vanderpol.jl`

Then, to plot run `./experiments/2_additional_derivatives/3_plot_all.jl`

<p align="center">
<img alt="Figure 3" src="./experiments/2_additional_derivatives/figure3_additional_derivative_info.svg" width="100%">
</p>


### Case study 3: Dynamical systems with conserved quantities
For figure 4, run (in order):
- `./experiments/3_conserved_quantities/1.1_generate_workprecision_data.jl`
- `./experiments/3_conserved_quantities/1.2_workprecision_plot.jl`

<p align="center">
<img alt="Figure 4" src="./experiments/3_conserved_quantities/figure4_conservedquantity_workprecision.svg" width="100%">
</p>

`./experiments/3_conserved_quantities/2_longtermplot.jl` to get figure 5:

<p align="center">
<img alt="Figure 5" src="./experiments/3_conserved_quantities/figure5_longtermplot.svg" width="80%">
</p>

And `./experiments/3_conserved_quantities/3_inflated_kepler_samples.jl` to get figure 6:

<p align="center">
<img alt="Figure 6" src="./experiments/3_conserved_quantities/figure6_conservedquantity_and_samples.svg" width="70%">
</p>


### Case study 4: Mass-matrix differential-algebraic equations
To run the experiments and make the plot, run (in order)
- `./experiments/4_massmatrices/1.1_pendulum.jl`
- `./experiments/4_massmatrices/1.2_rober.jl`
- `./experiments/4_massmatrices/2_plot.jl`

<p align="center">
<img alt="Figure 7" src="./experiments/4_massmatrices/figure7_massmatrix_dae.svg" width="70%">
</p>


## Reference
```
@InProceedings{pmlr-v151-bosch22a,
  title = 	 { Pick-and-Mix Information Operators for Probabilistic ODE Solvers },
  author =       {Bosch, Nathanael and Tronarp, Filip and Hennig, Philipp},
  booktitle = 	 {Proceedings of The 25th International Conference on Artificial Intelligence and Statistics},
  pages = 	 {10015--10027},
  year = 	 {2022},
  editor = 	 {Camps-Valls, Gustau and Ruiz, Francisco J. R. and Valera, Isabel},
  volume = 	 {151},
  series = 	 {Proceedings of Machine Learning Research},
  month = 	 {28--30 Mar},
  publisher =    {PMLR},
  pdf = 	 {https://proceedings.mlr.press/v151/bosch22a/bosch22a.pdf},
  url = 	 {https://proceedings.mlr.press/v151/bosch22a.html},
  abstract = 	 { Probabilistic numerical solvers for ordinary differential equations compute posterior distributions over the solution of an initial value problem via Bayesian inference. In this paper, we leverage their probabilistic formulation to seamlessly include additional information as general likelihood terms. We show that second-order differential equations should be directly provided to the solver, instead of transforming the problem to first order. Additionally, by including higher-order information or physical conservation laws in the model, solutions become more accurate and more physically meaningful. Lastly, we demonstrate the utility of flexible information operators by solving differential-algebraic equations. In conclusion, the probabilistic formulation of numerical solvers offers a flexible way to incorporate various types of information, thus improving the resulting solutions. }
}
```
