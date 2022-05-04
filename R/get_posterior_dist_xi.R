#' @title Posterior sampling of \eqn{\xi}
#' 
#' @description Clever way to sample from the posterior of \eqn{\xi} given draws
#' from the prior and the log-likelihood of data
#'
#'
#' @param params A tibble coming from calling `rstan::extract(x)$xi`, where
#' x is a stanfit object from either the horseshoe or horseshoe plus priors
#' @param n_draws An integer giving the desired number of draws from the posterior
#' distribution
#' @param data_grouped A tibble with three named columns: `n`, `resp`, and `dose_num`,
#' corresponding to the number of observations, number of responses among those
#' observations, and dose number label
#'
#' @return A tibble with the same number of columns as `params`
#' @export
#' 
#' @importFrom dplyr pull


get_posterior_dist_xi <- function(params, 
                                  n_draws = nrow(params) / 2, 
                                  data_grouped = NULL) {
  
  tparams = t(params)
  
  if(is.null(data_grouped)) {
    loglikelihood <- numeric(nrow(params))
  } else {
    loglikelihood <- 
      colSums(pull(data_grouped, resp) * log(tparams) + 
                (pull(data_grouped, n) - pull(data_grouped, resp)) * log(1 - tparams))
  }
  
  likelihood <- exp(loglikelihood - max(loglikelihood) + 10)
  
  foo <- params[sample(nrow(params), n_draws, replace = T, prob = likelihood),]
  
}

