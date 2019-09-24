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

#' Helper function to return the correct indices for latent variables
#'
#' @description The organization of the latent field is confusing. There
#' are random effects and fixed effects, plus the actual linear predictor values
#' (which you usually don't want).
#'
#' This function accepts a ccmodeldata object and returns a list telling you
#' which indices in the vector of latent variables correspond to which model elements.
#'
#' Details: the latent Gaussian variables have dimension Wd = Nd + M + p. Nd
#' is the total number of control days, sum(model_data$control_days), of all subjects.
#' M is the total number of random effects, i.e. unique values of smooth terms, in the
#' order in which they appear in the model. p is the number of regression coefficients.
#' To get the means/variances of the random effects, the correct call is
#' i = (model_data$Nd+1):(model_data$Nd+model_data$M). To get the regression coefficients,
#' the correct call is i = (model_data$Nd+model_data$M+1):(model_data$Nd+model_data$M+model_data$p).
#'
#' @param model_data ccmodeldata object as returned by model_setup().
#'
#' @export
#'
get_indices <- function(model_data) {
  model_elements <- model_data$model_elements
  out <- structure(list(), class = "ccindex")
  if (model_data$M == 0) {
    # No random effects.
    out$linear <- (model_data$Nd+1):(model_data$Nd+model_data$p)
    names(out$linear) <- model_elements$linear
    out$smooth <- c()
  } else if (model_data$p == 0) {
    # No linear terms
    out$smooth <- (model_data$Nd+1):(model_data$Nd+model_data$M)
    # Figure out the number of smooth effects for each covariate
    numterms <- model_data$A %>%
      purrr::map("u") %>%
      purrr::map(unique) %>%
      purrr::map(length)
    numtermsvec <- purrr::reduce(numterms,c)
    names(numtermsvec) <- names(numterms)
    # Account for the fact that we've potentially set one of the first
    # smooth covariate's effects to zero, through a linear constraint.
    if (!is.null(model_data$vectorofcolumnstoremove)) {
      if (model_data$vectorofcolumnstoremove != 0) {
        numtermsvec[1] <- numtermsvec[1] - 1
      }
    }

    names(out$smooth) <- rep(names(numtermsvec),numtermsvec)
    out$linear <- c()
  } else {
    # Both linear and smooth terms
    out$smooth <- (model_data$Nd+1):(model_data$Nd+model_data$M)
    numterms <- model_data$A %>%
      purrr::map("u") %>%
      purrr::map(unique) %>%
      purrr::map(length)
    numtermsvec <- purrr::reduce(numterms,c)
    names(numtermsvec) <- names(numterms)
    if (!is.null(model_data$vectorofcolumnstoremove)) {
      if (model_data$vectorofcolumnstoremove != 0) {
        numtermsvec[1] <- numtermsvec[1] - 1
      }
    }

    names(out$smooth) <- rep(names(numtermsvec),numtermsvec)
    out$linear <- (model_data$Nd+model_data$M+1):(model_data$Nd+model_data$M+model_data$p)
    names(out$linear) <- model_elements$linear
  }
  out
}


#' Compute marginal means and variances
#'
#' @description Given optimization results with log posteriors computed, compute the
#' marginal means and variances. These are the final outputs of the whole estimation
#' procedure.
#'
#' Linear constraints are corrected for at this point. Marginal variances for linear
#' combinations of latent variables are also available.
#'
#' @param i Either a 1) vector giving the indices of the latent variables for which you want marginal
#' means and variances or 2) an object of class ccindex output by get_indices() which prescribes which
#' terms you want means/variances for.
#' @param model_results Output of optimization; can be before or after you call add_log_posterior_values().
#' @param model_data ccmodeldata object output by model_setup()
#' @param constrA Either a sparse matrix whose columns contain linear constraints under which you would
#' like to compute means/variances, or NULL. If NULL, any linear constraints will be pulled from model_data.
#' @param lincomb Either a sparse matrix whose columns contain linear combinations of the latent variables
#' whose means/variances you would like to compute, or an object of class cclincomb output by make_model_lincombs()
#'
#' @export

compute_marginal_means_and_variances <- function(i,model_results,model_data,constrA = NULL,lincomb = NULL) {
  if (nrow(model_results) > 1) {
    # Add log posterior values for theta if not present
    if (!("theta_logposterior" %in% names(model_results))) {
      model_results <- add_log_posterior_values(model_results,model_data)
    }
    # Normalize
    thetanormconst <- normalize_log_posterior(model_results$theta_logposterior,model_results$theta)
    model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
    # Get the integration weights
    dx1 <- diff(model_results$theta) # dx1[i] = x[i+1] - x[i]
    ld <- length(dx1)
    intweights <- c(
      dx1[1]/2,
      (dx1[1:(ld-1)] + dx1[2:ld])/2,
      dx1[ld]/2
    )
  }

  # Compute the precision matrices for each theta
  precision_matrices <- model_results %>%
    purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                      theta = ..1)
    )

  # Compute the hessians for each theta
  hessian_structure <- hessian_log_likelihood_structure(W = model_results$solution[[1]],model_data = model_data)
  hessians <- model_results %>%
    purrr::pmap(~list(
      C = hessian_log_likelihood(W = ..4,model_data = model_data,structure = hessian_structure),
      theta = ..1)
    )
  # If linear combinations required, set up the relevant functions
  if (!is.null(lincomb)) {
    compute_var_one_lincomb <- function(a,Q) {
      # a <- cbind(a)
      ZZ <- solve(Q,a)
      as.numeric(crossprod(a,ZZ))
    }
    compute_var_all_lincombs <- function(A,Q) {
      # Coerce to list of sparse vectors
      AA <- list()
      for (j in 1:ncol(lincomb)) AA[[j]] <- as(lincomb[ ,j],"sparseVector")
      AA %>% purrr::map(~compute_var_one_lincomb(.x,Q)) %>% purrr::reduce(c)
    }
    compute_one_lincomb_correction <- function(a,WW,VV) {
      as.numeric(crossprod(crossprod(VV,a),crossprod(WW,a)))
    }
    compute_all_lincomb_correction <- function(A,WW,VV) {
      AA <- list()
      for (j in 1:ncol(lincomb)) AA[[j]] <- as(lincomb[ ,j],"sparseVector")
      AA %>% purrr::map(~compute_one_lincomb_correction(.x,WW,VV)) %>% purrr::reduce(c)
    }
  }

  # If no linear constraints, compute the marginal means and variances as normal
  if (is.null(constrA)) {
    margmeans <- model_results %>%
      purrr::pmap(~..4) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)

    # Marginal variances: add the precision and the hessian and get diagOfInv
    margvars <- purrr::map2(precision_matrices,hessians,~.x[["Q"]] + .y[["C"]]) %>%
      purrr::map(~diagOfInv(x = .x,constrA = NULL,i = i)) %>%
      purrr::reduce(rbind)

    # If there are linear combinations, compute their variances separately from diagOfInv
    if (!is.null(lincomb)) {
      # lincomb is a column matrix. Change to list and map over the columns
      lincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]])) %>%
        purrr::map("lincombvar") %>%
        purrr::reduce(rbind)
    }
  } else {
    # If there are linear constraints, compute the corrected mean and variance
    # First compute the uncorrected mean
    uncorrectedmean <- purrr::pmap(model_results,~list(theta = ..1,mode = ..4))

    # Get the precision matrix of the GMRF- Q + C
    QpC <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]]))

    # Compute the correction term
    WW <- purrr::map(QpC,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))

    YY <- purrr::map2(WW,uncorrectedmean,~list(
      YY = solve(t(constrA) %*% .x[["WW"]],t(constrA) %*% .y[["mode"]]),
      theta = .y[["theta"]]
    ))
    correction_mean <- purrr::map2(WW,YY,~list(correction = .x[["WW"]] %*% .y[["YY"]],theta = .x[["theta"]]))

    # Now correct
    margmeans <- purrr::map2(uncorrectedmean,correction_mean,
                             ~.x[["mode"]] - .y[["correction"]]) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)

    # Now compute the variances
    # Add the corrected mean to the model_results
    model_results$corrected_mean <- vector(mode = "list",length = nrow(model_results))
    for (k in 1:nrow(model_results)) model_results$corrected_mean[[k]] <- margmeans[k, ]

    # Re-compute the hessians
    corrected_hessians <- model_results %>%
      purrr::pmap(~list(
        C = hessian_log_likelihood(W = ..12,model_data = model_data,structure = hessian_structure),
        theta = ..1)
      )
    # Get the corrected precision matrix of the GMRF- Q + C_correct
    QpC_corrected <- purrr::map2(precision_matrices,corrected_hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]]))

    # uncorrectedvariances <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = NULL,i = i))
    margvars <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = constrA,i = i)) %>%
      purrr::reduce(rbind)
    if (!is.matrix(margvars)) margvars <- matrix(margvars,nrow = 1)

    # If we require marginal variances for linear combinations, compute them separately from diagOfInv
    if (!is.null(lincomb)) {
      uncorrectedlincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]]))

      # Compute the corrections
      WW <- purrr::map(QpC_corrected,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
      VV <- purrr::map(WW,~list(VV = solve(t(.x[["WW"]]) %*% constrA,t(.x[["WW"]])),theta = .x[["theta"]])) %>%
        purrr::map(~list(VV = t(.x[["VV"]]),theta = .x[["theta"]]))

      lincombvarcorrections <- purrr::map2(WW,VV,~list(lincombvarcorrection = compute_all_lincomb_correction(lincomb,.x[["WW"]],.y[["VV"]]),
                                                       theta = .x[["theta"]]))

      lincombvars <- purrr::map2(uncorrectedlincombvars,lincombvarcorrections,
                                 ~list(lincombvar = .x[["lincombvar"]] - .y[["lincombvarcorrection"]],
                                       theta = .x[["theta"]])) %>%
        purrr::map("lincombvar") %>%
        purrr::reduce(rbind)
    }


  }
  if (nrow(margmeans) == 1) {
    finalmeans <- as.numeric(margmeans)[i]
    finalvars <- as.numeric(margvars)
    finallincombvars <- NULL
    if (!is.null(lincomb)) finallincombvars <- as.numeric(lincombvars)

  } else {
    postvals <- exp(model_results$theta_logposterior + log(intweights))
    finalmeans <- sweep(margmeans,1,postvals,"*") %>% apply(2,sum)
    finalvars <- sweep(margvars,1,postvals,"*") %>% apply(2,sum)
    finallincombvars <- NULL
    if (!is.null(lincomb)) finallincombvars <- sweep(lincombvars,1,postvals,"*") %>% apply(2,sum)
    finalmeans <- finalmeans[i]
  }

  list(mean = finalmeans,
       variance = finalvars,
       lincombvars = finallincombvars)

}


