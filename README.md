# Dose Deescalation

This R package provides code to simulate the designs considered in 
*Targeted randomization de-escalation trials enable fractional dosing of scarce drugs*
by Boonstra, Tabarrok, and Strohbehn (2022). 

Install the package using `devtools::install_github("psboonstra/DoseDeescalation", build_vignettes = TRUE)`

There are four main functions in this package:

  1. `gen_curves()` generates a list of curves from the Hill equations
  2. `get_posterior_dist_xi()` is a clever way to sample from the posterior of $\xi$ given draws from the prior and the log-likelihood of data
  3. `get_posterior_prob_mdse()` is a dlever way to calculate the posterior probability of each dose level being the MDSE given draws from the prior and the log-likelihood of data
  4. `sim_model_deesc()` is the function to actually do the simulations
    
Run `vignette("RunSimulations", package = "DoseDeescalation")`


