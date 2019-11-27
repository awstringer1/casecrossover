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

  out <- optresults %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sigma = list(exp(-.5 * unlist(.data[["theta"]]))),
           sigma_logposterior = length(unlist(.data[["sigma"]])) * log(2) - sum(log(unlist(.data[["sigma"]]))) + .data[["theta_logposterior"]])

  attr(out,"thetagrid") <- attributes(optresults)$thetagrid
  out
}

#' Normalize log-posterior via simple numerical integration
#'
#' @description Simple implementation of the trapezoid rule for integration. The key is that
#' this function accepts grid points and LOGGED function values and computes the integral in
#' a numerically stable manner, using matrixStats::logSumExp().
#'
#' For a single hyperparameter, the function works on uniform and non-uniform grids. For
#' a multidimensional hyperparameter, the function only works on a uniform grid (in all dimensions).
#'
#' @param pp logged function values
#' @param tt grid points at which these function values were computed, as a NIGrid object
#'
#' @return A number giving the LOG of the integral of the function used to compute pp.
#' This is the log-normalizing constant.
#'
#' @export
#'
normalize_log_posterior <- function(pp,tt) {
  # tt: grid returned by mvQuad::createNIGrid
  # pp: log posterior evaluated at these points
  ww <- mvQuad::getWeights(tt)
  matrixStats::logSumExp(log(ww) + pp)
}

#' Normalize the log posterior returned by the optimization
#'
#' @description Wrapper around normalize_log_posterior() that takes in
#' the dataframe of optimization results and returns the same dataframe
#' but with the theta_logposterior column normalized.
#'
#' @param optresults Optimization results, a tibble() output by optimize_all_thetas_parallel().
#'
#' @return Data frame of optimization results with the theta_logposterior column normalized,
#' i.e. satisfying logSumExp(theta_logposterior) = 0
#'
#' @export
#'
normalize_optresults_logpost <- function(optresults) {
  thetanormconst <- normalize_log_posterior(optresults$theta_logposterior,attributes(optresults)$thetagrid)
  optresults$theta_logposterior <- optresults$theta_logposterior - thetanormconst

  out <- optresults %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sigma = list(exp(-.5 * unlist(.data[["theta"]]))),
                  sigma_logposterior = length(unlist(.data[["sigma"]])) * log(2) - sum(log(unlist(.data[["sigma"]]))) + .data[["theta_logposterior"]])

  attr(out,"thetagrid") <- attributes(optresults)$thetagrid

  out
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
#' All this is way too complicated for the user, and even for the developer, so here is a
#' function to parse the model data and return an object which enumerates which elements
#' of the (internal) latent field correspond to which quantities. This is mostly
#' used internally.
#'
#' @param model_data ccmodeldata object as returned by model_setup().
#' @param removezeros Optional. Should the indices be returned with hard-zeroes included or removed from the
#' random effect vector? Depends on where this function is being called.
#'
#' @return A list of class ccindex containing vectors of indices for the linear and smooth terms.
#'
#' @export
#'
get_indices <- function(model_data,removezeros = TRUE) {
  if (is.null(model_data$vectorofcolumnstoremove)) {
    removezeros <- FALSE
  } else if (length(model_data$vectorofcolumnstoremove) == 0) {
    removezeros <- FALSE
  } else if (all(model_data$vectorofcolumnstoremove == 0)) {
    removezeros <- FALSE
  }
  model_elements <- model_data$model_elements
  out <- structure(list(), class = "ccindex")
  if (model_data$M == 0) {
    # No random effects.
    out$linear <- (model_data$Nd+1):(model_data$Nd+model_data$p)
    names(out$linear) <- model_elements$linear
    out$smooth <- c()
  } else if (model_data$p == 0) {
    # No linear terms
    if (removezeros) {
      numsmooth <- model_data$M
    } else {
      numsmooth <- model_data$M + length(model_data$vectorofcolumnstoremove)
    }
    out$smooth <- (model_data$Nd+1):(model_data$Nd+numsmooth)
    covvalues <- model_data$A %>%
      purrr::map("u") %>%
      purrr::map(unique) %>%
      purrr::map(sort)
    numterms <- covvalues %>%
      purrr::map(length)
    numtermsvec <- purrr::reduce(numterms,c)
    names(numtermsvec) <- names(numterms)
    if (removezeros) {
      if (!is.null(model_data$vectorofcolumnstoremove)) {
        if (all(model_data$vectorofcolumnstoremove == 0)) {
          # Do nothing
        } else {
          for (vr in names(model_data$vectorofcolumnstoremove)) {
            numtermsvec[vr] <- numtermsvec[vr] - 1
            toremove <- which(covvalues[[vr]] == model_data$control$linear_constraints[[vr]]$whichzero[1])
            covvalues[[vr]] <- covvalues[[vr]][-toremove]
          }
        }
      }
    }

    names(out$smooth) <- rep(names(numtermsvec),numtermsvec)
    out$linear <- c()
    out$covvalues <- covvalues
  } else {
    # Both linear and smooth terms
    if (removezeros) {
      numsmooth <- model_data$M
    } else {
      numsmooth <- model_data$M + length(model_data$vectorofcolumnstoremove)
    }
    out$smooth <- (model_data$Nd+1):(model_data$Nd+numsmooth)
    covvalues <- model_data$A %>%
      purrr::map("u") %>%
      purrr::map(unique) %>%
      purrr::map(sort)
    numterms <- covvalues %>%
      purrr::map(length)
    numtermsvec <- purrr::reduce(numterms,c)
    names(numtermsvec) <- names(numterms)
    if (removezeros) {
      if (!is.null(model_data$vectorofcolumnstoremove)) {
        if (all(model_data$vectorofcolumnstoremove == 0)) {
          # Do nothing
        } else {
          for (vr in names(model_data$vectorofcolumnstoremove)) {
            numtermsvec[vr] <- numtermsvec[vr] - 1
            toremove <- which(covvalues[[vr]] == model_data$control$linear_constraints[[vr]]$whichzero[1])
            covvalues[[vr]] <- covvalues[[vr]][-toremove]
          }
        }
      }
    }

    names(out$smooth) <- rep(names(numtermsvec),numtermsvec)
    out$linear <- (model_data$Nd+numsmooth+1):(model_data$Nd + numsmooth +model_data$p)
    degrees <- get_polynomial_degree(model_elements$linear_formula)
    names(out$linear) <- rep(names(degrees),degrees)
    out$covvalues <- covvalues
  }
  out
}



#' Marginal means/variances for linear combinations of latent variables
#'
#' @description The compute_marginal_means_and_variances() function lets you specify a sparse matrix of
#' linear combinations of latent variables to compute the marginal means and variances of. The most
#' common use case for this is when you have included an effect as both linear and smooth, and you
#' need the marginal mean and variance of the linear predictor. You need to account for the correlation
#' between the posterior regression coefficient and posterior random effect.
#'
#' Because this is a common use case, this function takes in your model_data and returns a correctly
#' formatted matrix which specifies that you want means/variances for the linear predictors from all
#' terms which appear in the model both as linear and smooth effects.
#'
#' @param model_data ccmodeldata object output by model_setup()
#'
#' @return A sparse matrix with one column per necessary linear combination.
#'
#' @export
#'
make_model_lincombs <- function(model_data) {
  # Check to make sure both linear and smooth terms in the model
  if (length(model_data$model_elements$linear) == 0 | length(model_data$model_elements$smooth) == 0) {
    stop("You should only be looking at linear combinations if there are both linear and smooth terms in your model.")
  }
  # Check to see if the SAME term is included as both linear and smooth.
  if (length(intersect(model_data$model_elements$linear,model_data$model_elements$smooth)) == 0) {
    stop("You have linear and smooth terms in your model, but I don't see any overlap. You should only use this function to create linear combinations for terms in your model that are included as both linear and smooth")
  }

  # Determine the degree of polynomial for each covariate
  polydegrees <- get_polynomial_degree(model_data$model_elements$linear_formula)

  # Get the indices for all terms
  indices <- get_indices(model_data,removezeros = TRUE)
  indices_zeroes <- get_indices(model_data,removezeros = FALSE)
  # Wd <- model_data$Nd + length(indices$smooth) + length(indices$linear)
  Wd <- model_data$Wd
  smooth_terms <- unique(names(indices$smooth))
  linear_terms <- unique(names(indices$linear))
  # The terms we use are the terms that appear both in smooth and linear
  terms_to_use <- intersect(smooth_terms,linear_terms)

  # Create the linear combinations
  lincomblist <- list()
  outmatlist <- list()
  for (nm in terms_to_use) {
    degree <- polydegrees[names(polydegrees) == nm]
    idx <- indices$smooth[names(indices$smooth) == nm]
    linear_idx <- indices$linear[names(indices$linear) == nm]
    u <- indices$covvalues[[nm]]

    for (j in 1:length(u)) {
      betavec <- u[j]^(1:degree)
      ll <- sparseVector(
        x = c(1,betavec),
        i = c(idx[j],linear_idx),
        length = Wd
      )
      lincomblist <- c(lincomblist,ll)
    }
    outmat <- lincomblist %>%
      purrr::map(~as(.,"sparseMatrix")) %>%
      purrr::reduce(cbind)

    # Correct for removed zeros
    if (!is.null(model_data$vectorofcolumnstoremove)) {
      if (length(model_data$vectorofcolumnstoremove) > 0 & !all(model_data$vectorofcolumnstoremove == 0)) {
        constr <- model_data$control$linear_constraints[[nm]]
        wcol <- which(constr$u == constr$whichzero[1])
        nr <- nrow(outmat)
        bindvec <- as(sparseVector(x = constr$u[wcol]^(1:degree),i = linear_idx,length = Wd),"sparseMatrix")
        if (wcol == 1) {
          outmat <- cbind(bindvec,outmat)
        } else if (wcol == ncol(outmat) + 1) {
          outmat <- cbind(outmat,bindvec)
        } else {
          outmat <- cbind(outmat[ ,1:(wcol-1)],bindvec,outmat[ ,wcol:ncol(outmat)])
        }
      }
    }

    outmatlist <- c(outmatlist,outmat)
    lincomblist <- list()
  }

  # Note: the following is not the fastest way to do this. See stackoverflow:
  # https://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors#8844057
  # It's about 2x - 3x faster in benchmarking; not worth introducing new code.
  outmatlist %>% purrr::reduce(cbind)
}


#' Make a matrix of linear constraints
#'
#' @description Take the linear constraints from the control argument of a ccmodeldata object
#' and create a sparse matrix of linear constraints.
#'
#' @param model_data ccmodeldata object created by model_setup. Must have control$linear_constraints
#'
#' @return A sparse matrix where each column represents one linear constraint. If the return value is A
#' then AW = 0.
#'
#' @export
#'
make_linear_constraints <- function(model_data) {
  if (is.null(model_data$control$linear_constraints)) stop("No linear constraints provided in model_data$control")
  if (length(model_data$control$linear_constraints) == 0) stop("No linear constraints provided in model_data$control")

  # Pull the constraints from the model data, retaining only those that
  # have more than one value. Single-value constraints are automatically done in the precision matrix
  constr <- model_data$control$linear_constraints %>% purrr::map("whichzero") %>%
    purrr::keep(~length(.x) > 1)
  if (length(constr) == 0) return(0) # No additional constraints beyond what was manually set to zero.
  # Get the indices of all unique covariate values that haven't already been set to zero:
  idx <- get_indices(model_data)
  # For the constraints that appear in covvalues (because they haven't already been
  # manually removed), figure out which element of the covvalues it is, then take
  # that element from the idx$smooth for that covariate (!)
  get_constraint_index <- function(nm) {
    cvec <- constr[[nm]]
    ivec <- idx$covvalues[[nm]]
    idxvec <- idx$smooth[names(idx$smooth) == nm]
    cvec <- cvec[cvec %in% ivec]
    out <- numeric(length(cvec))
    for (j in 1:length(cvec)) out[j] <- idxvec[which(cvec[j] == ivec)]
    unname(out)
  }
  # Create a sparse matrix of constraints
  names(constr) %>%
    purrr::map(get_constraint_index) %>%
    purrr::map(~Matrix::sparseVector(x = 1,i = .x,length = model_data$Wd)) %>%
    purrr::map(as,Class = "sparseMatrix") %>%
    purrr::reduce(cbind) %>%
    as(Class = "dgTMatrix")

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
#' whose means/variances you would like to compute, or an object of class cclincomb output by make_model_lincombs().
#' If NULL, will be computed automatically. Set lincomb = FALSE in order to prevent this.
#'
#'
#' @export
#'

# i <- index11
# model_results <- opt_11
# model_data <- model_data11

compute_marginal_means_and_variances <- function(model_results,model_data,i = NULL,constrA = NULL,lincomb = NULL) {
  # If a ccindex object provided, change to numeric
  if (is.null(i)) i <- get_indices(model_data)
  if (class(i) == "ccindex") {
    idx <- c(i$smooth,i$linear)
  } else if (is.numeric(i)) {
    idx <- i
  } else {
    stop(stringr::str_c("i must be an object of class ccindex, or a numeric vector. You provided an object of class ",class(i)))
  }
  idx <- unname(idx)

  # Pull constraints from model_data
  if (is.null(constrA)) {
    # Check if model has any linear constraints
    if (!is.null(model_data$control$linear_constraints)) {
      if (length(model_data$control$linear_constraints) > 0) {
        constrA <- make_linear_constraints(model_data)
      }
    }
  } else if (!constrA) {
    # If constrA is FALSE, nullify it now
    constrA <- NULL
  } else {
    if (!inherits(constrA,"sparseMatrix")) {
      stop("constrA must either be NULL, logical, or an object inheriting from sparseMatrix")
    }
  }

  # Pull linear combinations from model data
  if (is.null(lincomb)) {
    # Check if the model has both linear and smooth terms for the SAME covariate
    if (length(intersect(model_data$model_elements$smooth,model_data$model_elements$linear)) > 0) {
      lincomb <- make_model_lincombs(model_data)
    }
  } else if (is.logical(lincomb)) {
    # If lincomb set to FALSE, nullify it now
    if (!lincomb) lincomb <- NULL
  } else {
    if (!inherits(lincomb,"sparseMatrix")) {
      stop("lincomb must either be NULL, logical, or an object inheriting from sparseMatrix")
    }
  }

  # Function can take a tibble with only a single row, representing the "eb" setting from INLA.
  # This is mostly to accommodate models with only linear terms.
  if (nrow(model_results) > 1) {
    # Add log posterior values for theta if not present
    if (!("theta_logposterior" %in% names(model_results))) {
      model_results <- add_log_posterior_values(model_results,model_data)
    }
    # Normalize
    thetanormconst <- normalize_log_posterior(model_results$theta_logposterior,attributes(model_results)$thetagrid)
    model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
    # Note that doing this still leaves sigma_logposterior unnormalized! But it's not used in this function.

    # Get the integration weights
    # dx1 <- diff(model_results$theta) # dx1[i] = x[i+1] - x[i]
    # ld <- length(dx1)
    # intweights <- c(
    #   dx1[1]/2,
    #   (dx1[1:(ld-1)] + dx1[2:ld])/2,
    #   dx1[ld]/2
    # )
    intweights <- mvQuad::getWeights(attributes(model_results)$thetagrid)[ ,1]
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
  if (is.null(constrA) | all(constrA == 0)) {
    margmeans <- model_results %>%
      purrr::pmap(~..4) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)

    # Marginal variances: add the precision and the hessian and get diagOfInv
    margvars <- purrr::map2(precision_matrices,hessians,~.x[["Q"]] + .y[["C"]]) %>%
      purrr::map(~diagOfInv(x = .x,constrA = NULL,i = idx)) %>%
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
        C = hessian_log_likelihood(W = ..10,model_data = model_data,structure = hessian_structure),
        theta = ..1)
      )
    # Get the corrected precision matrix of the GMRF- Q + C_correct
    QpC_corrected <- purrr::map2(precision_matrices,corrected_hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]]))

    # uncorrectedvariances <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = NULL,i = i))
    margvars <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = constrA,i = idx)) %>%
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
    finalmeans <- as.numeric(margmeans)[idx]
    finalvars <- as.numeric(margvars)
    finallincombvars <- NULL
    if (!is.null(lincomb)) finallincombvars <- as.numeric(lincombvars)

  } else {
    postvals <- exp(model_results$theta_logposterior + log(intweights))
    # postvals <- exp(model_results$theta_logposterior)
    finalmeans <- sweep(margmeans,1,postvals,"*") %>% apply(2,sum)
    finalvars <- sweep(margvars,1,postvals,"*") %>% apply(2,sum)
    finallincombvars <- NULL
    if (!is.null(lincomb)) finallincombvars <- sweep(lincombvars,1,postvals,"*") %>% apply(2,sum)
    finalmeans <- finalmeans[idx]
  }

  # If elements of the latent field were manually removed, now is the time to add them
  # back in!
  if (!is.null(model_data$vectorofcolumnstoremove)) {
    if (!all(model_data$vectorofcolumnstoremove == 0)) {
      finalmeans <- stitch_zero_vector(finalmeans,model_data$vectorofcolumnstoremove)
      finalvars <- stitch_zero_vector(finalvars,model_data$vectorofcolumnstoremove)
    }
  }

  list(mean = finalmeans,
       variance = finalvars,
       lincombvars = finallincombvars,
       lincombmatrix = lincomb)
}

#' Marginal hyperparameter posterior
#'
#' @description Get the marginal density of a single coordinate of the hyperparameter vector theta using
#' quadrature.
#'
#' @param j Integer index indicating which marginal to evaluate. If theta is K-dimensional then this should be from 1 to K.
#' @param ccfit Object of class ccfit returned by casecrossover()
#' @param quantiles Vector of numbers between 0 and 1 representing quantiles of the marginal posterior that are to be estimated
#'
#' @return A list containing elements
#'         - margpost: dataframe of logged theta_j posterior values, evaluated on a grid of the same accuracy as those used
#' to construct the full joint.
#'         - quantiles: dataframe containing the quantiles of the marginal posterior. Easier to do these here since they
#'         require a grid.
#'
#' @export
#'
marginal_hyperparameter_posterior <- function(j,ccfit,quantiles = c(2.5,97.5)/100) {
  thetagridfull <- attr(ccfit$optimization,"thetagrid")
  S <- thetagridfull$dim
  # If it's already one-dimensional, don't need to do anything new, but do compute quantiles
  if (S == 1) {
    outmat <- dplyr::tibble(
      theta = purrr::reduce(ccfit$optimization$theta,rbind),
      thetalogmargpost = ccfit$optimization$theta_logposterior,
      sigma = purrr::reduce(ccfit$optimization$sigma,rbind),
      sigmalogmargpost = ccfit$optimization$sigma_logposterior
    )

    thetacumsum <- cumsum(mvQuad::getWeights(thetagridfull) * exp(outmat$thetalogmargpost))
    thetaquantiles <- purrr::map(quantiles,~outmat$theta[which(thetacumsum == min(thetacumsum[thetacumsum > .x]))]) %>% purrr::reduce(c)
    sigmaquantiles <- rev(exp(-.5 * thetaquantiles))

    return(
      list(
        margpost = outmat,
        quantiles = dplyr::tibble(whichmarginal = rep(1,length(quantiles)),q = quantiles,theta = thetaquantiles,sigma = sigmaquantiles)
      )
    )
  }
  # Get the reduced grid
  thetagridreduced <- mvQuad::createNIGrid(
    dim = thetagridfull$dim - 1,
    type = thetagridfull$type[-j],
    level = as.numeric(thetagridfull$level[ ,-j]),
    ndConstruction = thetagridfull$ndConstruction,
    level.trans = thetagridfull$level.trans
  )
  mvQuad::rescale(thetagridreduced,domain = thetagridfull$features$domain[-j, ])

  # Get a 1-d grid, for computing quantiles at the end.
  thetagrid1d <- mvQuad::createNIGrid(
    dim = 1,
    type = thetagridfull$type[j],
    level = as.numeric(thetagridfull$level[ ,j]),
    ndConstruction = thetagridfull$ndConstruction,
    level.trans = thetagridfull$level.trans
  )
  mvQuad::rescale(thetagrid1d,domain = thetagridfull$features$domain[j, ])

  # Return the marginal posterior and its evaluation points
  # In the optimization results, we have a matrix of theta values which matches the full grid,
  # and the log posterior evaluated at these values.
  nodesfull <- purrr::reduce(ccfit$optimization$theta,rbind)
  nodesfull <- cbind(nodesfull,ccfit$optimization$theta_logposterior) # Theta logposterior is now the last column
  nodesfull <- nodesfull[order(nodesfull[ ,j]), ]
  colnames(nodesfull) <- c(paste0("theta",1:S),"thetalogpost")
  nodesfull <- as.data.frame(nodesfull)

  # Now we have a matrix of thetavalues, nodes, and function values-- add on the weights
  nodesmulti <- mvQuad::getNodes(thetagridreduced)
  nodesmulti <- cbind(nodesmulti,mvQuad::getWeights(thetagridreduced))
  colnames(nodesmulti) <- c(paste0("theta",(1:S)[-j]),"weights")
  nodesmulti <- as.data.frame(nodesmulti)

  suppressMessages({# It prints what it's joining by, which is all columns, and I don't want to see this printed
    thetamargposts <- dplyr::left_join(nodesfull,nodesmulti) %>%
      dplyr::group_by(.data[[stringr::str_c("theta",j)]]) %>%
      dplyr::summarize(thetalogmargpost = matrixStats::logSumExp(.data[["thetalogpost"]] + log(.data[["weights"]])))
  })

  thetamargposts$whichmarginal <- rep(j,nrow(thetamargposts))
  # Now add on the sigmas
  outmat <- thetamargposts %>%
    dplyr::mutate(sigma = exp(-.5 * .data[[paste0("theta",j)]]),
                  sigmalogmargpost = log(2/.data[["sigma"]]) + .data[["thetalogmargpost"]]
    ) %>%
      dplyr::rename(theta = .data[[paste0("theta",j)]])

  # Quantiles
  thetacumsum <- cumsum(mvQuad::getWeights(thetagrid1d) * exp(outmat$thetalogmargpost))
  thetaquantiles <- purrr::map(quantiles,~outmat$theta[which(thetacumsum == min(thetacumsum[thetacumsum > .x]))]) %>% purrr::reduce(c)
  sigmaquantiles <- rev(exp(-.5 * thetaquantiles))

  list(
    margpost = outmat,
    quantiles = dplyr::tibble(whichmarginal = rep(j,length(quantiles)),q = quantiles,theta = thetaquantiles,sigma = sigmaquantiles)
  )
}


