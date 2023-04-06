#' @title Simulate many dose deescalation trials
#' 
#' @description This function will simulate an arbitrary number
#' of dose deescalation trials for a given design configuration
#'
#'
#' @param n_sim Number of simulations to conduct
#' @param true_eff_curve Vector of non-decreasing probabilities denoting the
#' true response probabilities. The order is presumed to give the corresponding
#' dose level
#' @param relative_threshold A number between 0 and 1 indicating the smallest
#' relative efficacy as compared to the largest dose level that would still
#' be considered acceptable
#' @param cohort_sizes The numbers of subjects to enroll in each cohort before
#' re-estimating the assigned dose level. The total sample size of a
#' trial will be `sum(cohort_sizes)`. The initial cohort, which requires
#' assignments prior to having any data, will also be deterministically
#' allocated in descending order, starting from the top dose and proceeding
#' downwards, recycling as necessary. 
#' @param cohort_dose_span Non-negative integers giving the number of dose 
#' levels below the current best dose (see `be_greedy`) to randomize to. Should
#' be equal in length to `cohort_sizes`. Set to 0 to always choose the best dose. 
#' Any integers equal to or greater than `length(true_eff_curve)` will be
#' soft-truncated above at `length(true_eff_curve)-1`
#' @param be_greedy A logical. If TRUE, the best dose is defined as the
#' lowest dose level with posterior probability of being the MDSE greater
#' than one over the number of dose levels. If FALSE, the best dose is 
#' defined as the dose level with highest posterior probability of being
#' the MDSE
#' @param prior_type A character, either 'horseshoe' or 'horseshoe_plus', 
#' indicating the type of prior to use
#' @param mc_warmup Number of warmup iterations for the sampler. Note
#' the sampler will only be run once to sample from the prior; the posterior
#' is approximated using these prior draws and the likelihood. Thus, `mc_warmup`
#' and `mc_iter` (below) can and probably should be a lot larger than usual
#' @param mc_iter Number of total iterations for the sampler. See also
#' `mc_warmup` above 
#' @param mc_chains Number of separate chains to run.  See also
#' `mc_warmup` above
#' @param seed Random seed to set
#'
#' @return A named list with three objects:
#' 
#' `subject_results` is a tibble containg subject-level information: the dose
#' assignment (`dose_num`), the subject ID (`subj_id`), the recorded response
#' (`resp`), the true response probability at the subject's assigned dose level
#' (`true_xi`), the estimated posterior probability that their assigned dose
#' number was the MDSE just prior to their enrollment (`prob_mdse_just_before`),
#' the estimated posterior mean of the response probability at the subject's
#' assigned dose level **as of the end of the trial** (`est_xi`), the estimated
#' posterior probability that the subject's assigned dose level is true MDSE
#' **as of the end of the trial** (`prob_mdse`). Other information
#' that is not specific to the subject is an integer label for each simulated
#' trial (`sim_id`), the prior used (`prior_type`).
#' 
#' 
#' `mdse_results` is a tibble containing trial-level information: `sim_id`,
#' `prior_type`, a character giving the cohort sizes separated by commas
#' (`cohort_sizes`), a character giving the dose spans used for each cohort
#' separated by commas (`cohort_dose_span`), the total planned sample size of
#' the trial (`n_patients`), the true MDSE in the trial (`true_mdse`), the
#' estimated MDSE as of the end of the trial (`est_mdse`), the estimated
#' **absolute** response probability at the estimated MDSE as of the end of the
#' trial (`est_xi_at_est_mdse`), the estimated **relative** response probability
#' at the estimated MDSE as compared to the response probability at the largest
#' dose level (`est_rel_xi_at_est_mdse`), the estimated posterior probability
#' that the subject's assigned dose level is true MDSE as of the end of the
#' trial (`prob_mdse`), the actual number of subjects enrolled in the trial
#' (`enrollment`), the number of subjects enrolled to the dose number that was
#' the estimated MDSE at the end of the trial (`num_at_est_mdse`), below the
#' estimated MDSE (`num_below_est_mdse`), above the estimated MDSE
#' (`num_above_est_mdse`), at the **true** MDSE (`num_at_true_mdse`), below the
#' **true** MDSE (`num_below_true_mdse`), and above the true MDSE
#' (`num_above_true_mdse`), and the average of the true response probabilities
#' across all subjects (`avg_true_xi`).
#' 
#' @export
#'
#' @examples 
#' # See vignette for examples
#' 
#' @importFrom dplyr %>% bind_cols bind_rows mutate pull group_by summarize arrange filter select left_join n
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr drop_na
#' @importFrom tidyselect everything
#' @importFrom stats rbinom

sim_model_deesc <- 
  function(n_sim,
           true_eff_curve,
           relative_threshold,
           cohort_sizes = c(34, 33, 33),
           cohort_dose_span = c(3, 5, 5),
           be_greedy = FALSE,
           prior_type = c("horseshoe", "horseshoe_plus"),
           mc_warmup = 1e4, 
           mc_iter = 2.1e5,
           mc_chains = 3,
           seed = sample(.Machine$integer.max,1)) {
    
    
    begin = Sys.time() 
    
    prior_type = match.arg(prior_type)
    stopifnot(isTRUE(prior_type %in% c("horseshoe", "horseshoe_plus")))
    stopifnot(isTRUE(length(cohort_sizes) == length(cohort_dose_span)))
    stopifnot(all(true_eff_curve >= 0) && all(true_eff_curve <= 1))
    stopifnot(relative_threshold > 0 && relative_threshold <= 1)
    # ensure non-negative integers
    stopifnot(all(cohort_dose_span >= 0) && all(cohort_dose_span %% 1 == 0))
    # force all spans to be no more than number of dose levels
    cohort_dose_span = pmin(cohort_dose_span, length(true_eff_curve) - 1)
    
    n_subjects = sum(cohort_sizes);
    num_dose_levels = length(true_eff_curve)
    subject_results = NULL
    mdse_results = 
      cbind(sim_id = 1:n_sim,
            matrix(nrow = n_sim, ncol = 11, 
                   dimnames = list(NULL,
                                   c("est_mdse",
                                     "est_xi_at_est_mdse",
                                     "est_rel_xi_at_est_mdse",
                                     "prob_mdse",
                                     "enrollment",
                                     "num_at_est_mdse",
                                     "num_below_est_mdse", 
                                     "num_above_est_mdse",
                                     "num_at_true_mdse",
                                     "num_below_true_mdse", 
                                     "num_above_true_mdse"))))
    
    
    true_mdse = 
      min(which(true_eff_curve / max(true_eff_curve) >= relative_threshold))
    
   curr_fit <- sampling(object = stanmodels[[paste0("iso_",prior_type)]],
                         data = list(n_groups_stan = num_dose_levels,
                                     n_per_group_stan = as.array(numeric(num_dose_levels)),
                                     y_stan = as.array(numeric(num_dose_levels)),
                                     local_dof_stan = 1, 
                                     global_dof_stan = 1,
                                     alpha_scale_stan = 1e-5,
                                     only_prior_stan = 0),
                         warmup = mc_warmup,
                         iter = mc_iter,
                         chains = mc_chains,
                         verbose = F,
                         refresh = 0)
    
    params = rstan::extract(curr_fit)$xi
    
    get_new_dose <- 
      function(
    posterior_prob_mdse,
    num_pat,
    dose_span
      ) {
        
        if(!be_greedy) {
          best_dose <- which.max(posterior_prob_mdse)
        } else {    
          best_dose <- min(which(posterior_prob_mdse >= 1 / length(posterior_prob_mdse)))
        }
        
        new_dose_options <- 
          intersect(seq_along(posterior_prob_mdse), 
                    c(best_dose + c(0, seq(-1, -num_dose_levels, by = -1)))[1:(1+dose_span)]) %>%
          sort()
        
        if(length(new_dose_options) == 1) {
          rep(new_dose_options, length = num_pat)
        } else {
          sample(new_dose_options, num_pat, prob = posterior_prob_mdse[new_dose_options], replace = T)
        }
      }
    
    
    
    set.seed(seed);
    
    for(i in 1:n_sim) {
      posterior_prob_mdse <-
        get_posterior_prob_mdse(params = params, 
                                relative_threshold = relative_threshold,
                                data_grouped = NULL)
      
      
      initial_assignments = 
        rep((num_dose_levels : (num_dose_levels - cohort_dose_span[1])),
            length = cohort_sizes[1])
      initial_outcomes =
        rbinom(cohort_sizes[1], 1, true_eff_curve[initial_assignments])
      
      curr_dose = min(initial_assignments)
      
      curr_results_ungrouped = 
        tibble(dose_num = initial_assignments, 
               subj_id = 1:cohort_sizes[1],
               resp = initial_outcomes, 
               true_xi = true_eff_curve[initial_assignments],
               prob_mdse_just_before = posterior_prob_mdse[initial_assignments])
      
      if(cohort_sizes[1] < n_subjects) {
        curr_results_ungrouped <- 
          bind_rows(
            curr_results_ungrouped,
            tibble(dose_num = NA_integer_, 
                   subj_id = (cohort_sizes[1]+1):(n_subjects),
                   resp = NA_integer_,
                   true_xi = NA_real_)
          )
      }
      
      curr_results <- 
        curr_results_ungrouped %>%
        drop_na(dose_num) %>%
        mutate(dose_num = factor(dose_num, levels = 1:num_dose_levels)) %>%
        group_by(dose_num, .drop = FALSE) %>%
        summarize(n = n(), 
                  resp = sum(resp)) %>%
        mutate(dose_num = as.integer(dose_num)) %>%
        arrange(dose_num)
      
      j = 2
      curr_last_subject = cohort_sizes[1]
      while(j <= length(cohort_sizes)) {
        posterior_prob_mdse <-
          get_posterior_prob_mdse(params = params, 
                                  relative_threshold = relative_threshold,
                                  data_grouped = curr_results)
        
        curr_dose <- 
          get_new_dose(posterior_prob_mdse = posterior_prob_mdse,
                       num_pat = cohort_sizes[j], 
                       dose_span = cohort_dose_span[j])
        
        new_outcome <- rbinom(cohort_sizes[j], 1, true_eff_curve[curr_dose])
        
        #
        
        curr_results_ungrouped[(curr_last_subject+1):(curr_last_subject+cohort_sizes[j]),] <- 
          tibble(dose_num = curr_dose, 
                 subj_id = (curr_last_subject+1):(curr_last_subject+cohort_sizes[j]), 
                 resp = new_outcome, 
                 true_xi = true_eff_curve[curr_dose],
                 prob_mdse_just_before = posterior_prob_mdse[curr_dose])
        
        curr_results <- 
          curr_results_ungrouped %>%
          drop_na(dose_num) %>%
          mutate(dose_num = factor(dose_num, levels = 1:num_dose_levels)) %>%
          group_by(dose_num, .drop = FALSE) %>%
          summarize(n = n(), 
                    resp = sum(resp), 
                    .groups = "drop") %>%
          mutate(dose_num = as.integer(dose_num)) %>% 
          arrange(dose_num)
        
        curr_last_subject = curr_last_subject + cohort_sizes[j]
        j = j + 1
      }
      rm(j, curr_last_subject)
      
      posterior_prob_mdse <-
        get_posterior_prob_mdse(params = params, 
                                relative_threshold = relative_threshold,
                                data_grouped = curr_results)
      
      posterior_dist_xi <- 
        get_posterior_dist_xi(params = params, 
                              data_grouped = curr_results)
      
      posterior_mean_xi <-
        colMeans(posterior_dist_xi)
      
      posterior_mean_relative_xi <- 
        colMeans(posterior_dist_xi / posterior_dist_xi[,num_dose_levels])
      
      
      mdse_results[i,"est_mdse"] <-
        which.max(posterior_prob_mdse)
      mdse_results[i,"est_xi_at_est_mdse"] <-
        posterior_mean_xi[which.max(posterior_prob_mdse)]
      mdse_results[i,"est_rel_xi_at_est_mdse"] <-
        posterior_mean_relative_xi[which.max(posterior_prob_mdse)]
      mdse_results[i,"prob_mdse"] <-
        max(posterior_prob_mdse)
      mdse_results[i,"enrollment"] <-
        sum(curr_results %>% pull(n))
      #
      mdse_results[i,"num_at_est_mdse"] <-
        curr_results %>% 
        filter(dose_num == which.max(posterior_prob_mdse)) %>%
        pull(n)
      mdse_results[i,"num_below_est_mdse"] <-
        sum(curr_results %>% 
              filter(dose_num < which.max(posterior_prob_mdse)) %>%
              pull(n))
      mdse_results[i,"num_above_est_mdse"] <-
        sum(curr_results %>% 
              filter(dose_num > which.max(posterior_prob_mdse)) %>%
              pull(n))
      #
      mdse_results[i,"num_at_true_mdse"] <-
        curr_results %>% 
        filter(dose_num == true_mdse) %>%
        pull(n)
      mdse_results[i,"num_below_true_mdse"] <-
        sum(curr_results %>% 
              filter(dose_num < true_mdse) %>%
              pull(n))
      mdse_results[i,"num_above_true_mdse"] <-
        sum(curr_results %>% 
              filter(dose_num > true_mdse) %>%
              pull(n))
      
      
      
      subject_results <- 
        bind_rows(
          subject_results, 
          bind_cols(sim_id = i, 
                    curr_results_ungrouped %>%
                      left_join(
                        tibble(dose_num = seq_along(posterior_prob_mdse),
                               prob_mdse = posterior_prob_mdse, 
                               est_xi = posterior_mean_xi),
                        by = "dose_num"))
        )
      if(i %% max(1, round(n_sim / 10)) == 0) {cat(paste0(formatC(100 * i / n_sim), "% complete\n"));}
      rm(curr_dose, initial_outcomes, curr_results, curr_results_ungrouped)
    }
    
    mdse_results <- 
      bind_cols(prior_type = prior_type,
                cohort_sizes = paste0(cohort_sizes,collapse = ","),
                cohort_dose_span = paste0(cohort_dose_span,collapse = ","),
                n_subjects = n_subjects,
                true_mdse = true_mdse,
                as_tibble(mdse_results) %>%
                  left_join(
                    subject_results %>%
                      group_by(sim_id) %>%
                      summarize(avg_true_xi = mean(true_xi)), 
                    by = "sim_id")) %>%
      select(sim_id, everything())
    
    subject_results <- 
      bind_cols(prior_type = prior_type, 
                subject_results) %>%
      select(sim_id, everything())
    
    
    return(list(subject_results = subject_results,
                mdse_results = mdse_results,
                seed = seed, 
                running_time = Sys.time() - begin));
    
  }

