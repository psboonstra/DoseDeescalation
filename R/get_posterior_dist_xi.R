#' @title Posterior sampling of \eqn{\xi}
#' 
#' @description Clever way to sample from the posterior of \eqn{\xi} given draws
#'   from the prior and the log-likelihood of data
#'
#'
#' @param params A matrix that normally comes from calling
#'   `rstan::extract(x)$xi`, where `x` is a stanfit object from either
#'   `iso_horseshoe` or `iso_horseshoe_plus`. The rows of this matrix represent
#'   draws from the *prior* distribution, and the columns represent the
#'   individual components of \eqn{\xi}.
#' @param n_draws An integer giving the desired number of draws from the
#'   posterior distribution.
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
#' @return A matrix containing draws from the posterior. It will have the same
#'   number of columns as `params` and `n_draws` rows.
#' 
#' @export
#' 


get_posterior_dist_xi <- function(params, 
                                  n_draws = nrow(params) / 2, 
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
  
  likelihood <- exp(loglikelihood - max(loglikelihood) + 10)
  
  params[sample(nrow(params), n_draws, replace = TRUE, prob = likelihood),]
  
}

