# multimethod_removal
Estimating animal abundance using multiple removal methods

Below are descriptions of each of the scripts/functions used to fit multimethod removal models or simulate data from them. This model is described in full by Davis et al. 2022 in Ecological Applications.

## rem_mod_script.R 
Main script for implementing all functions described here. Specify simulation parameters and/or data, specify model structures, fit model using MCMC, conduct posterior predictive checks.

## rem_sim.R 
Simulates data from multimethod removal model to test MCMC on. User can specify number of sites, periods, removals, a time-varying covariate, a spatially-varying covariate, and model structure.

## rem_mod_setup.R 
Specifies removal model and translates into nimble code, constants, and initial values for running MCMC algoirthm and fitting model to data. 

## run_mcmc.R
General nimble setup function that fits model to data, using outputs from rem_mod_setup and data or rem_sim.

## traceplot_fun.R
Generates traceplots and cumulative mean plots to assess chain convergence.

## postpred_checks.R
Generates posterior predictive plots that plot model-predicted abundance over time and 95% credible intervals.

## storing_log_probs.R
Internal nimble function used in nimble MCMC to store log probabilities of data, used in model validation. Sourced in run_mcmc and in nimble model code.

## avail_fun.R 
Internal nimble function for performing availability correction within nimble-generated MCMC. Sourced in run_mcmc and in nimble model code.
