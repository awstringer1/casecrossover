### casecrossover function ###
# Full, user-facing function to fit case crossover models.
# Stitches together all the previous steps.

#' Fit a case crossover model with linear and non-parametric terms.
#'
#' @description Fit a case crossover model as described in Stringer et. al. (2019).
#' A case crossover model is used to quantify the association between mortality risk
#' and short-term exposure to environmental risk factors such as extreme temperatures
#' or air pollution. Multiple exposure records are available for a collection of
#' subjects who have died, including exposure on the death day and one or more
#' previous days. Higher exposure on date of death indicates a positive association
#' between exposure and mortality risk.
#'
#' This function fits a case crossover model with linear and/or non-parametric terms.
#' It has a typical R formula interface, y ~ x, with the following special terms allowed:
#'
#' - s(x) fits a smooth (non-parametric) risk curve to covariate x, using a RW2 model,
#' - strata(id) indicates that the variable id groups together records from the same
#' subject. Your formula must have an id.
#'
#' See the package vignette for detailed examples.
#'
#' @param formula An R formula, which can contain s() terms and must contain a strata() term.
#' @param data An object inheriting from class data.frame, typically a tibble or a data.table.
#'             Must contain multiple rows for each subject, with an id variable which states
#'             which records are from the same subject. Must have one case day and one or more
#'             control days per subject. There should be a case variable which equals 0 for control
#'             days and equals a positive integer for case days. See the vignette.
#' @param control An object created using cc_control(). This is used to specify priors and linear
#'                constraints for the s() terms. You can also set internal control parameters for the
#'                optimization.
#' @param verbose Logical, print progress updates and debugging information? Default FALSE.
#'
#' @return An object of class "cc_fit", with methods for summary() and plot().
#'
#' @export
#'

# # formula <- case1 ~ s(x) + s(x2) + poly(x,2,raw = TRUE) + poly(x2,3,raw=TRUE) + strata(id)
# formula <- case1 ~ x + strata(id)
# data <- sampledata
# # control = controlsmooth2
# control = cc_default_control()

casecrossover <- function(formula,data,control = cc_default_control(),verbose = FALSE) {
  # Set up the model
  if (verbose) cat("Setting up model...\n")
  model_data <- model_setup(formula,data,control)

  # Sort the data by id, case
  if (verbose) cat("Sorting data...\n")
  data <- data %>% dplyr::arrange(.data[[model_data$model_elements$strata]],.data[[model_data$model_elements$response]])

  # Choose the theta grid
  # TODO: add this to control
  if (verbose) cat("Creating the grid for numerical integration of theta posterior...\n")
  K <- length(model_data$model_elements$smooth)
  if (K == 0) {
    # Create a blank grid as a placeholder
    thetagrid <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 1)
  } else {
    thetagrid <- mvQuad::createNIGrid(dim = K,type = "GLe",level = model_data$control$thetaaccuracy)
    if (is.matrix(model_data$control$thetarange)) {
      mvQuad::rescale(thetagrid,domain = model_data$control$thetarange)
    } else if (is.numeric(model_data$control$thetarange)) {
      mvQuad::rescale(thetagrid,domain = matrix(rep(model_data$control$thetarange,K),ncol = K,byrow = TRUE))
    } else {
      stop("Specify thetarange either as a matrix with 2 columns, or as a vector of length dim(theta).")
    }
  }

  # Optimization
  if (verbose) cat("Performing optimization...\n")
  opt <- optimize_all_thetas_parallel(thetagrid,
                                      model_data,
                                      optcontrol = model_data$control$opt_control,
                                      doparallel = model_data$control$doparallel)

  # Post-hoc quantities
  if (verbose) cat("Computing post-hoc quantities...\n")
  # Add log posterior values before computing means/variances, so they can be returned.
  optwithlogpost <- add_log_posterior_values(opt,model_data) %>% normalize_optresults_logpost()
  posthoc <- compute_marginal_means_and_variances(optwithlogpost,model_data)

  # If we removed columns from the precision matrix, we added the appropriate zeroes back in inside
  # compute_marginal_means_and_variances(). So now, update the dimensions
  # if (!is.null(model_data$vectorofcolumnstoremove)) {
  #   ll <- length(model_data$vectorofcolumnstoremove)
  #   model_data$M <- model_data$M + ll
  #   model_data$Wd <- model_data$Nd + model_data$M + model_data$p
  #   # model_data$vectorofcolumnstoremove <- NULL # No longer needed
  # }

  # Build final output object, for plotting and printing etc.
  if (verbose) cat("Building output object...\n")
  out <- list(
    posthoc = posthoc,
    optimization = optwithlogpost,
    thetagrid = thetagrid,
    modeldata = model_data
  )
  structure(out,class = "cc_fit")
}

#' Summarize results of a casecrossover fit.
#'
#' @method summary cc_fit
#'
#' @description Generic summary method for a cc_fit object returned by casecrossover(). Prints similar
#' information to the familiar lm() and glm() summary functions, or inla().
#'
#' @param object Object of class cc_fit returned by casecrossover().
#' @param ... Not used
#'
#' @return A list of summary information of class cc_summary, with a print method.
#'
#' @export
#'
summary.cc_fit <- function(object,...) {
  idx <- get_indices(object$modeldata,removezeros = FALSE)

  call <- object$modeldata$model_elements$call

  summarytablefixed <- summarytablerandom <- summarytablelincomb <- summarytablehyper <- margpostsummarylist <- NULL

  if ("linear" %in% names(idx)) {
    linearidx <- idx$linear - object$modeldata$Nd
    mn <- object$posthoc$mean[linearidx]
    sd = sqrt(object$posthoc$variance[linearidx])
    summarytablefixed <- data.frame(
      mean = mn,
      sd = sd,
      q2.5 = stats::qnorm(.025,mean = mn,sd = sd),
      q97.5 = stats::qnorm(.975,mean = mn,sd = sd)
    )
    summarytablefixed$covariate <- colnames(object$modeldata$X)
  }

  if ("smooth" %in% names(idx)) {
    # Random effect summary
    smoothidx <- idx$smooth - object$modeldata$Nd
    mn <- object$posthoc$mean[smoothidx]
    sd <- sqrt(object$posthoc$variance[smoothidx])

    summarytablerandom <- data.frame(
      mean = mn,
      sd = sd,
      q2.5 = stats::qnorm(.025,mean = mn,sd = sd),
      q97.5 = stats::qnorm(.975,mean = mn,sd = sd)
    )
    # Create the proper rownames
    covvalues <- purrr::reduce(idx$covvalues,c)
    # covvalues <- stitch_zero_vector(covvalues,object$modeldata$vectorofcolumnstoremove)
    covnames <- names(idx$smooth)
    summarytablerandom$covariate <- covnames
    summarytablerandom$covariate_value <- covvalues

    summarytablerandom <- summarytablerandom[ ,c("covariate","covariate_value","mean","sd","q2.5","q97.5")]

    # Hyperparameter summary
    ww <- mvQuad::getWeights(attr(object$optimization,"thetagrid"))[ ,1]
    thetamean <- apply(ww * purrr::reduce(object$optimization$theta,rbind) * exp(object$optimization$theta_logposterior),2,sum)
    sigmamean <- apply(ww * purrr::reduce(object$optimization$sigma,rbind) * exp(object$optimization$sigma_logposterior),2,sum)
    thetasd <- sqrt(apply(ww * (purrr::reduce(object$optimization$theta,rbind) - thetamean)^2 * exp(object$optimization$theta_logposterior),2,sum))
    sigmasd <- sqrt(apply(ww * (purrr::reduce(object$optimization$sigma,rbind) - sigmamean)^2 * exp(object$optimization$sigma_logposterior),2,sum))

    # Quantiles based on the marginal posterior
    S <- length(thetamean)
    margpostsummarylist <- purrr::map(1:S,~marginal_hyperparameter_posterior(.x,object,c(2.5,97.5)/100))
    margquants <- margpostsummarylist %>%
      purrr::map("quantiles") %>%
      purrr::reduce(dplyr::bind_rows)

    q2.5theta <- dplyr::filter(margquants,.data[["q"]] == 0.025) %>% dplyr::pull(.data[["theta"]])
    q97.5theta <- dplyr::filter(margquants,.data[["q"]] == 0.975) %>% dplyr::pull(.data[["theta"]])
    q2.5sigma <- dplyr::filter(margquants,.data[["q"]] == 0.025) %>% dplyr::pull(.data[["sigma"]])
    q97.5sigma <- dplyr::filter(margquants,.data[["q"]] == 0.975) %>% dplyr::pull(.data[["sigma"]])

    summarytablehyper <- data.frame(
      covariate = rep(unique(names(idx$smooth)),2),
      variable = rep(c("theta","sigma"),each = length(unique(names(idx$smooth)))),
      mean = c(thetamean,sigmamean),
      sd = c(thetasd,sigmasd),
      q2.5 = c(q2.5theta,q2.5sigma),
      q97.5 = c(q97.5theta,q97.5sigma)
    ) %>%
      dplyr::arrange(.data[["covariate"]],.data[["variable"]])
  }

  if (!is.null(object$posthoc$lincombvars)) {
    # Linear combinations mean and variance
    lincombmatrix <- make_model_lincombs(object$modeldata) # This now correctly reflects the presence of the added-back zeroes
    summarytablelincomb <- cbind(
      as.numeric(t(lincombmatrix[(nrow(lincombmatrix) - length(object$posthoc$mean) + 1):nrow(lincombmatrix), ]) %*% cbind(object$posthoc$mean)),
      sqrt(object$posthoc$lincombvars)
    ) %>% as.data.frame()
    colnames(summarytablelincomb) <- c("mean","sd")
    summarytablelincomb$covariate <- covnames
    summarytablelincomb$covariate_value <- covvalues
  }

  out <- list(
    call = call,
    summarytablefixed = summarytablefixed,
    summarytablerandom = summarytablerandom,
    summarytablehyper = summarytablehyper,
    summarytablelincomb = summarytablelincomb,
    hypermarginals = margpostsummarylist,
    idx = idx
  )
  structure(out,class = "cc_summary")
}

# summary(pollutioncc)

#' Print method for summary.cc_fit.
#'
#' @method print cc_summary
#'
#' @description Generic print method for a cc_fit summary object. Not exported, user will call this function
#' by simply typing summary(...) in the console, like normal.
#'
#' @param x An object of class cc_summary, output by applying summary() to an object of class cc_fit returned by casecrossover().
#'
#' @return Nothing.
#'
print.cc_summary <- function(x) {
  cat("Call:\n")
  print(x$call)

  if ("linear" %in% names(x$idx)) {
    cat("\n\nFixed effects:\n")
    print(x$summarytablefixed)
  }
  if ("smooth" %in% names(x$idx)) {
    cat("\n\nRandom effects:\n")
    print(x$summarytablerandom)
    cat("\n\nHyperparameter(s):\n")
    print(x$summarytablehyper)
  }

  if (!is.null(x$summarytablelincomb)) {
    if (nrow(x$summarytablelincomb) > 10) {
      cat("\n\nSuppressing printing of linear combinations due to size. Access them with summary(...)$summarytablelincomb.\n")
    } else {
      cat("\n\nLinear combinations:\n")
      print(x$summarytablelincomb)
    }
  }
}

#' Plot results of fitting a case crossover model.
#'
#' @method plot cc_fit
#'
#' @description Plots corresponding to a case crossover model. This function computes histograms of
#' fixed-effect posteriors and smooth curves with pointwise error bars (+- 2sd) based on the marginal
#' posteriors of the smooth effects.
#'
#' @param x Object of class cc_fit returned by casecrossover().
#' @param ... Not used.
#'
#' @return A named list of lists of ggplots, one for the fixed effects and one for the smooth.
#' Only mimimal annotations are made; you can put your own titles and axis labels and further change the plot
#' features using the ggplot2 API. The plot lists themselves can be passed directly to cowplot::plot_grid,
#' or further annotated by you.
#'
#' @export
#'
plot.cc_fit <- function(x,...) {
  idx <- get_indices(x$modeldata,removezeros = FALSE)

  plotlist <- list()

  summ <- summary(x)

  if ("linear" %in% names(idx)) {
    linearidx <- idx$linear - x$modeldata$Nd
    plotlist$linear <- list()
    for (nm in names(idx$linear)) {
      mn <- x$posthoc$mean[linearidx[nm]]
      sd <- sqrt(x$posthoc$variance[linearidx[nm]])
      plt <- dplyr::tibble(x = c(mn - 3.5*sd,mn + 3.5*sd)) %>%
        ggplot2::ggplot(ggplot2::aes(x = .data[["x"]])) +
        ggplot2::theme_classic() +
        ggplot2::stat_function(fun = stats::dnorm,args = list(mean = mn,sd = sd)) +
        ggplot2::labs(title = "",x = nm,y = "density") +
        ggplot2::theme(text = ggplot2::element_text(size = 12))

      plotlist$linear[[nm]] <- plt
    }
  }

  if ("smooth" %in% names(idx)) {
    # Plot the posterior covariate effects
    plotlist$smooth <- list()
    covs <- summ$summarytablerandom$covariate %>% unique()
    plotlist$hyperparameters <- list()
    for (nm in covs) {
      vals <- sort(unique(x$modeldata$control$linear_constraints[[nm]]$u))
      plt <- summ$summarytablerandom %>%
        dplyr::filter(.data[["covariate"]] == nm) %>%
        dplyr::mutate(x = vals) %>%
        ggplot2::ggplot(ggplot2::aes(x = .data[["x"]])) +
        ggplot2::theme_classic() +
        ggplot2::geom_line(ggplot2::aes(y = .data[["mean"]])) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["q2.5"]],ymax = .data[["q97.5"]]),colour = "lightgrey",alpha = .1) +
        ggplot2::labs(title = "",x = nm,y = "Posterior mean and 95% CI") +
        ggplot2::theme(text = ggplot2::element_text(size = 12))

      plotlist$smooth[[nm]] <- plt
    }
    # Theta and sigma posteriors
    plotlist$hyperparameters$theta[[nm]] <- summ$hypermarginals[[which(covs == nm)]]$margpost %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["theta"]],y = exp(.data[["thetalogmargpost"]]))) +
      ggplot2::theme_classic() +
      ggplot2::geom_line() +
      ggplot2::labs(title = "",x = "theta",y = "density") +
      ggplot2::theme(text = ggplot2::element_text(size = 12))

    plotlist$hyperparameters$sigma[[nm]] <- summ$hypermarginals[[which(covs == nm)]]$margpost %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["sigma"]],y = exp(.data[["sigmalogmargpost"]]))) +
      ggplot2::theme_classic() +
      ggplot2::geom_line() +
      ggplot2::labs(title = "",x = "sigma",y = "density") +
      ggplot2::theme(text = ggplot2::element_text(size = 12))
  }

  # Linear combinations
  if (!is.null(x$posthoc$lincombmatrix)) {
    plotlist$linearcombinations <- list()
    covs <- summ$summarytablelincomb$covariate %>% unique()
    for (nm in covs) {
      vals <- sort(unique(x$modeldata$control$linear_constraints[[nm]]$u))
      plt <- summ$summarytablelincomb %>%
        dplyr::filter(.data[["covariate"]] == nm) %>%
        dplyr::mutate(x = vals) %>%
        ggplot2::ggplot(ggplot2::aes(x = .data[["x"]])) +
        ggplot2::theme_classic() +
        ggplot2::geom_line(ggplot2::aes(y = .data[["mean"]])) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["mean"]] - 2 * .data[["sd"]],ymax = .data[["mean"]] + 2 * .data[["sd"]]),colour = "lightgrey",alpha = .1) +
        ggplot2::labs(title = "",x = nm,y = "Posterior mean and 95% CI") +
        ggplot2::theme(text = ggplot2::element_text(size = 12))

      plotlist$linearcombinations[[nm]] <- plt
    }
  }

  structure(plotlist,class = c("cc_plot","list"))
}

