# Code for "Extreme value theory for singular subspace estimation in the matrix denoising model" by Junhyung Chang and Joshua Cape

### Introduction

The files in this repository are used to generate figures and tables in the paper. 
All files with names ending in `_mat.m` generate simulated data. They should be run first (this might take a while), then the corresponding files that do not have `_mat` at the end should run last.

### Functions

* `generate_U1.m`: inputs an orthonormal matrix U and a value t in [0, 1], and outputs an orthonormal matrix U1 that is closer to U as t is closer to 0. 
* `gumbel_quantile.m`: computes standard Gumbel quantile values.
* `procrustes_sol_F.m`: inputs two orthonormal matrices of the same size and computes the Frobenius norm optimal orthogonal Procrustes problem solution.
* `sbmmean.m`: outputs a block-structured matrix that has constant entries in each block.
* `tinorm.m`: computes the two-to-infinity norm of a matrix.

### Main convergence plots (Figure 1)

Description: create the three histogram + QQ-plot pairs in Figure 1.

* Step 1: run `main_convergence_mat.m` to create simulated data (this will create a new data file named `main_convergence.mat`).
* Step 2: run `main_convergence_plot.m` to create the plots.

### Power comparison contour plot (Figure 2)

Description: create two contour plots for empirical power of our test and Frobenius norm-based test, respectively.

* Step 1: run `pwr_comparison_mat.m` to create simulated data (this will create a new data file named `pwr_comparison.mat`).
* Step 2: run `pwr_comparison_plot.m` to create the plots.

### Power comparison table (Table 2)

Description: generate entries for the power comparison table in Table 2.

* Run `pwr_comparison_table.m`.

### Power analysis contour plot (Figure 3)

Description: create contour plot for Figure 3. Comment out line 26 for weak signal. Comment out line 25 for strong signal.

* Step 1: run `pwr_contour_mat.m` to create simulated data (this will create a new data file named `pwr_contour.mat`).
* Step 2: run `pwr_contour_plot.m` to create the plots.

### t-distributed noise (Figure 4)

Description: create the three histogram + QQ-plot pairs in Figure 4.

* Step 1: run `t_dist_mat.m`.
* Step 2: run `t_dist_plot.m`.


### Experiments with stochastic block model (Figures 2,3 and 4 in Supplement)

Description: create three histogram + QQ-plot pairs for Bernoulli, Poisson, Gaussian SBMs each. In the following steps, replace X with {bern, pois, norm} for each SBM edge distribution.

* Step 1: run `sbm_X_mat.m`.
* Step 2: run `sbm_X_plot.m`.