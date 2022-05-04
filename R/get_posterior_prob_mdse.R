#' @title Posterior probability of being minimum dose with satisfactory efficacy (MDSE)
#' 
#' @description Clever way to calculate the posterior probability of
#' each dose level being the MDSE given draws
#' from the prior and the log-likelihood of data
#'
#'
#' @param params A tibble coming from calling `rstan::extract(x)$xi`, where
#' x is a stanfit object from either the horseshoe or horseshoe plus priors. 
#' @param relative_threshold A number between 0 and 1 indicating the smallest
#' relative efficacy as compared to the largest dose level that would still
#' be considered acceptable. 
#' @param data_grouped A tibble with three named columns: `n`, `resp`, and `dose_num`,
#' corresponding to the number of observations, number of responses among those
#' observations, and dose number label
#'
#' @return A simplex vector of probabilities
#' @export
#'
#' @importFrom dplyr pull


get_posterior_prob_mdse <- function(params,
                                   relative_threshold, 
                                   data_grouped = NULL) {
  
  tparams = t(params)
  
  
  if(is.null(data_grouped)) {
    loglikelihood <- numeric(nrow(params))
  } else {
    loglikelihood <- 
      colSums(pull(data_grouped, resp) * log(tparams) + 
                (pull(data_grouped, n) - pull(data_grouped, resp)) * log(1 - tparams))
  }
  
  foo1 = (params >= params[,ncol(params)] * relative_threshold)
  foo2 = cbind(T, (params[,-ncol(params)] < params[,ncol(params)] * relative_threshold))
  
  likelihood <- exp(loglikelihood - max(loglikelihood) + 10)
  
  colSums(foo1 * foo2 * likelihood) / sum(likelihood)
  
}

