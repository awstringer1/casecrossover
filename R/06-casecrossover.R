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
    mvQuad::rescale(thetagrid,domain = matrix(c(rep(-10,K),rep(10,K)),ncol=K))
  }

  # Optimization
  if (verbose) cat("Performing optimization...\n")
  opt <- optimize_all_thetas_parallel(thetagrid,
                                      model_data,
                                      optcontrol = model_data$control$opt_control,
                                      doparallel = model_data$control$doparallel)

  # Post-hoc quantities
  if (verbose) cat("Computing post-hoc quantities...\n")
  posthoc <- compute_marginal_means_and_variances(opt,model_data)

  # Build final output object, for plotting and printing etc.
  if (verbose) cat("Building output object...\n")
  out <- list(
    posthoc = posthoc,
    optimization = opt,
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
  idx <- get_indices(object$modeldata)

  call <- object$modeldata$model_elements$call

  summarytablefixed <- data.frame(
    mean = object$posthoc$mean,
    sd = sqrt(object$posthoc$variance),
    q2.5 = stats::qnorm(.025,mean = object$posthoc$mean,sd = sqrt(object$posthoc$variance)),
    q97.5 = stats::qnorm(.975,mean = object$posthoc$mean,sd = sqrt(object$posthoc$variance))
  )
  rownames(summarytablefixed) <- colnames(object$modeldata$X)
  out <- list(
    call = call,
    summarytablefixed = summarytablefixed,
    idx = idx
  )
  structure(out,class = "cc_summary")
}

# summary(pollutioncc)

#' Print method for summary.cc_fit.
#'
#'
#' @description Generic print method for a cc_fit summary object. Not exported, user will call this function
#' by simply typing summary(...) in the console, like normal.
#'
#' @param ccsummary An object of class cc_summary, output by summary().
#'
#' @return Nothing.
#'
print.cc_summary <- function(ccsummary) {
  cat("Call:\n")
  print(ccsummary$call)

  if ("linear" %in% names(ccsummary$idx)) {
    cat("\n\nFixed effects:\n")
    print(ccsummary$summarytablefixed)
  }
  if ("smooth" %in% names(ccsummary$idx)) {
    cat("\n\nRandom effects:\n")
    print(ccsummary$summarytablefixed)
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
  idx <- get_indices(x$modeldata)

  plotlist <- list()

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
  structure(plotlist,class = c("cc_plot","list"))
}

