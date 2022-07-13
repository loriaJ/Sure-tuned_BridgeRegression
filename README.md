# Sure-tuned_BridgeRegression
**"SURE-tuned bridge regression" Authors: Jorge Loría, Anindya Bhadra.** 
Preprint: []()
# include the link to the preprint

The functions to replicate the simulations are in the files: 
- `functions_correlated_design_matrices.R`: has the functions that generates the simulations and runs the comparisons for a specific setting (e.g., number of observations, coefficient size, number of coefficients, etc...),
- `run_correlated_design_matrices.R`: which calls the functions to run the simulations.
These are used several times, using bash calls, to perform the simulations shown in the preprint.

The main functions to run a linear regression are in: `multivariate_sure_no_svd.R`. Specifically, the main function to use is: `t_multivariate_sure_schains`. This function takes four required parameters:

- `alpha`: the penalty coefficient,
- `x`: the design matrix,
- `y`: the observed responses,
- `sigma2`: the assumed variance.

It can also receive the optionals parameters: 
- `n_sim`: the number of simulations, predefined to be 10000,
- `nu_vect`: the values over which it will attempt the optimization, by default is the vector: 10^(-20:3),
- `nchains`: the number of separate simulations it will run ($C$, in Section ? of the manuscript), by default is 2.

The function to generate polynomially tilted positive alpha stable random variables ($PS(\alpha,1/2)$) are in the file: `t_alpha_beta_simulation.R`, which follows the procedure of Devroye (2008) for simulating these variables.




