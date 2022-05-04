#' @title Generate Hill equations
#' 
#' @description Generate a list of curves from the Hill equations
#' 
#' @param number_dose_levels Positive integer giving number of doses per curve
#' @param a The 'a' parameter of the Hill equation. Must be a single number
#' @param b_seq The 'b' parameter. May be a vector, which will result in one curve per value
#' @param c The 'c' parameter. Must be a single number
#' @param d_seq The 'd' parameter. May be a vector, which will result in one curve per value
#'
#' @return A list with length `length(b_seq) * length(d_seq)`, each with `number_dose_levels`
#' elements, giving the true probabilities at equal-spaced points along the x-axis
#' 
#' @export
#' 
#' @importFrom dplyr %>% 
#' @importFrom purrr map map2


gen_curves <- function(number_dose_levels = 6,
                       a = 0.2, 
                       b_seq = c(0.85, 0.7, 0.55, 0.4),
                       c = 0.7,
                       d_seq = c(1, 2, 3, 4, 5, 7, 10, 20)) {
  
  four_par_logistic = function(x, a, b, c, d) {
    b + (a - b) / (1 + (x / c) ^ d)
  }
  
  four_par_logistic_inv = function(y, a, b, c, d) {
    c * ((a - b) / (y - b) - 1) ^ (1 / d)
  }
  
  curve_list <- NULL
  
  for(curr_b in b_seq) {
    
    curve_list <- 
      c(curve_list,
        map(
          # 1. given fixed asymptote b, identify dose level x that yields 99% of asymptote (one for each value of d)
          four_par_logistic_inv(0.99 * curr_b, a = a, b = curr_b, c = c, d = d_seq),
          # 2. construct equidistant set of dose levels between 0 and x 
          seq, 
          from = 0,
          length = number_dose_levels) %>%
          # 3. map back to the probability scale 
          map2(.y = d_seq,
               function(x, y) four_par_logistic(x = x, a = a, b = curr_b, c = c, d = y)) %>%
          map(round, digits = 4)
      )
  }
  
  curve_list
}

