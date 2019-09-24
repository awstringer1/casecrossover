### Compute post-hoc quantities ###
# Compute densities, means and variances.



#' Approximation to hyperparameter log-posterior
#'
#' @description This function computes log pi(theta|y) for all thetas in the grid.
#'
#' @param optresults Optimization results, a tibble() output by optimize_all_thetas_parallel().
#' @param model_data ccmodeldata object output by model_setup()
#'
#' @return A tibble() in the same format as optresults, with columns added for the log-posterior
#' of theta and sigma.
#'
#' @export
#'
add_log_posterior_values <- function(optresults,model_data) {
  optresults <- dplyr::ungroup(optresults)
  hessian_structure <- hessian_log_likelihood_structure(rep(0,model_data$Wd),model_data)
  # Log posterior for theta
  logposttheta <- optresults %>%
    purrr::pmap(~log_posterior_theta(unlist(..1),unlist(..4),model_data,hessian_structure)) %>%
    purrr::reduce(c)
  optresults$theta_logposterior <- logposttheta
  optresults <- optresults %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sigma = list(exp(-.5 * unlist(.data[["theta"]]))),
           sigma_logposterior = length(unlist(.data[["sigma"]])) * log(2) - sum(log(unlist(.data[["sigma"]]))) + .data[["theta_logposterior"]])
  optresults
}

#' Normalize log-posterior via simple numerical integration
#'
#' @description Simple implementation of the trapezoid rule for integration. The key is that
#' this function accepts grid points and LOGGED function values and computes the integral in
#' a numerically stable manner, using matrixStats::logSumExp().
#'
#' @param pp logged function values
#' @param tt grid points at which these function values were computed
#'
#' @return A number giving the LOG of the integral of the function used to compute pp.
#' This is the log-normalizing constant.
#'
#' @export
#'
normalize_log_posterior <- function(pp,tt) {
  # pp: log posterior
  # tt: theta values at which pp is evaluated
  # function returns the LOG of the normalizing constant

  # Make sure tt is sorted
  df <- dplyr::tibble(pp = pp,tt = tt) %>% dplyr::arrange(tt)
  tt <- df$tt
  pp <- df$pp

  lp <- length(pp)
  matrixStats::logSumExp(c(
    matrixStats::logSumExp(pp[-1] + log(diff(tt)) + log(1/2)),
    matrixStats::logSumExp(pp[-lp] + log(diff(tt)) + log(1/2))
  ))
}



