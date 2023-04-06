#' @title Posterior probability of being minimum dose with satisfactory efficacy (MDSE)
#' 
#' @description Clever way to calculate the posterior probability of each dose
#'   level being the MDSE given draws from the prior and the log-likelihood of
#'   data
#'
#'
#' @param params A matrix that normally comes from calling
#'   `rstan::extract(x)$xi`, where `x` is a stanfit object from either
#'   `iso_horseshoe` or `iso_horseshoe_plus`. The rows of this matrix represent
#'   draws from the *prior* distribution, and the columns represent the
#'   individual components of \eqn{\xi}.
#' @param relative_threshold A number between 0 and 1 indicating the smallest
#'   relative efficacy as compared to the largest dose level that would still be
#'   considered acceptable. Its default value is set to 0.8, which is used in the
#'   manuscript.
#' @param data_grouped A tibble with two named columns: `n` and `resp`,
#'   corresponding to the number of observations and number of responses among
#'   those observations. The number of rows must be equal to the number of
#'   columns in `params`, corresponding to the number of dose levels. It is
#'   assumed that the rows are ordered by dose number, with the first row
#'   corresponding to the smallest/lowest dose and the last row corresponding to
#'   the largest/highest dose. Additional columns are allowed and will be
#'   ignored. `data_grouped` can alternatively be `NULL`, in which case you will
#'   get draws from the priors.
#'
#' @return A simplex vector of probabilities, interpreted as the posterior
#'   probability that each dose level is the MDSE.
#' @export
#'


get_posterior_prob_mdse <- function(params,
                                    relative_threshold = 0.8, 
                                    data_grouped = NULL) {
  
  tparams <- t(params)
  
  if(is.null(data_grouped)) {
    loglikelihood <- numeric(nrow(params))
  } else {
    stopifnot(nrow(data_grouped) == ncol(params))
    stopifnot(c("n", "resp") %in% colnames(data_grouped))
    stopifnot(all(data_grouped$n >= data_grouped$resp))
    
    loglikelihood <- 
      colSums(data_grouped$resp * log(tparams) + 
                (data_grouped$n - data_grouped$resp) * log(1 - tparams))
  }
  
  foo1 <- (params >= params[,ncol(params)] * relative_threshold)
  foo2 <- cbind(T, (params[,-ncol(params)] < params[,ncol(params)] * relative_threshold))
  
  likelihood <- exp(loglikelihood - max(loglikelihood) + 10)
  
  colSums(foo1 * foo2 * likelihood) / sum(likelihood)
  
}

