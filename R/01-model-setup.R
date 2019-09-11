### MODEL SETUP ###
# This file contains functions which take in data and
# model specifications and return objects usable
# by the internal likelihood/optimization/summary functions
# in the casecrossover package.

## Main model setup
# This subsection contains the main function for model setup. It is a
# "level 2" function, meaning it will be exported to users, but will only
# need to be used by intermediate/advanced users. It will be called by the
# "level 1" casecrossover::casecrossover() function, which is the main function
# users will use.
# Included in this section are several "level 3" functions, which are completely
# internal and not exported. They are only documented in an ad-hoc way.
#
# This section also includes some special imports that I couldn't figure out
# how to put anywhere else :/

#' Set control parameters for case-crossover models
#'
#' @description This function returns a list of control parameters
#' suitable for input into case-crossover model setup and fitting functions.
#' You can call it directly to see what parameters should be set and what their
#' default values are. You can pass any of these parameters to the function
#' directly to set their values. Arguments include:
#'
#' smooth_prior: Specify the prior on the smooth terms in your model. This should
#' be a named list with elements corresponding to each smooth term in the model. Each list
#' element should further be a named list naming the prior and its parameters. For an
#' explanation of available priors and their parameters, see the "priors" vignette.
#' Note that the package will throw an error if priors aren't specified for smooth terms.
#' There are no default priors.
#'
#' @param ... Arguments used to override the defaults output by cc_control().

#'
#' @examples
#' cc_control(smooth_prior = list(z = list(prior = "pc_prec",parameters = c(u = 3,a = 1))))
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @import stats
#' @import utils
#' @export
cc_control <- function(...) {
  outargs <- cc_default_control()
  if (length(list(...)) == 0) return(outargs)
  for (nm in names(list(...))) {
    outargs[[nm]] <- list(...)[[nm]]
  }
  outargs
}

#' Set up a casecrossover model for fitting.
#'
#' @description Set up a casecrossover model for fitting: creates all internal objects and stores them in a list with class "cc_model_data".
#' May contain standard linear terms "x" parsable
#' by model.matrix() and smooth terms "s(x)" to which RW2 models are fit
#' (see documentation for s()).
#'
#' The data will be sorted using
#' dplyr::arrange(id,case) where id and case are the names of the variables identifying
#' the subject and the case day respectively.
#'
#' The case column in your dataframe should
#' contain a 0 for control days, and a positive integer for case days. If the case day is
#' a 1, then the model is fit exactly as described in the paper; if it is a positive
#' integer greater than 1, then this is equivalent to having a dataset where that subject
#' is replicated that many times (i.e. the full log-likelihood for that subject is
#' multiplied by case). This is a convenience to save dataset storage space/memory.
#'
#' @param formula A standard R formula; see description.
#' Must contain a strata() term which defines the variable in "data" used to
#' distringuish individual subjects.
#' @param data A data.frame or tibble containing at minimum one column which groups
#' subjects together, one column indicating the case day, and one or more columns
#' containing covariates.
#' @examples
#' model_setup(y~x,data.frame(y = c(0,1),x = c(1,2)))
#'
#' @export
# model_setup <- function(formula,data,control = cc_control()) {
#   # TODO: implement this.
#   # TODO: document this.
#   model_data <- structure(list(), class = "ccmodeldata")
#
#   # Parse the formula
#   model_elements <- parse_formula(formula)
#   # Check that the smooth and strata terms exist in the data
#   # The linear terms will be passed to model.matrix, which has its own
#   # error checking.
#   extra_model_vars <- model_elements[c("smooth","strata")] %>% purrr::reduce(c)
#   if (!all(extra_model_vars %in% colnames(data))) {
#     missing_vars <- extra_model_vars[!(extra_model_vars %in% colnames(data))]
#     stop(paste0("The following variables were provided in the model formula but not in the data: ",
#                 stringr::str_c(missing_vars,collapse = ", ")))
#   }
#
#   # Create the smooth terms- design matrix
#   Alist <- list()
#   if (length(model_elements$smooth) > 0) {
#     for (nm in names(model_elements$smooth)) {
#       Alist[[nm]] <- create_alist_element(data[[nm]])
#     }
#   }
#   model_data$A <- Alist
#   model_data$M <- Alist %>% purrr::map("A") %>% purrr::map(ncol) %>% purrr::reduce(sum)
#   # TODO: when implementing the linear constraint part, where one element of the
#   # precision matrix is set to 0, make sure to make the appropriate reduction in M
#
#   # Number of subjects
#   n <- length(unique(data[model_elements$strata]))
#
#   # Linear terms
#   model_data$X <- Matrix::sparse.model.matrix(model_elements$linear_formula,data = data)
#   model_data$p <- ncol(model_data$X)
#
#   # Create the vector of control days
#   # A named vector where the names are the subject ids and the values are the number
#   # of control days that each has in the data
#
#   control_days <- data %>%
#     filter(.data[[model_elements$response]] == 0) %>%
#     group_by(stratum) %>%
#     summarize(control_days = n())
#
#   model_data$control_days <- control_days
#   names(model_data$control_days) <- sim1data$subject %>% unique() %>% sort()
#   model_data$Nd <- sum(model_data$control_days)
#   model_data$Ne <- model_data$Nd + model_data$n
#   model_data$Wd <- model_data$M + model_data$p + model_data$Nd
#   model_data$Wdf <- model_data$M + model_data$p + model_data$Ne
#   # Priors. Prior u is P(sigma > u) = alpha. It's the prior median if alpha = .5
#   # So set u = log(1.25), 50% chance that a unit increase in exposure yields a
#   # 25% increase in risk.
#   model_data$theta_logprior <- function(theta,prior_alpha = .75,prior_u = log(20)) {
#     # In this model, theta is the LOG PRECISION of the rw2 smoothing variance
#     # Implement the PC prior directly.
#     # P(sigma > u) = alpha.
#     # See inla.doc("pc.prec")
#     lambda <- -log(prior_alpha)/prior_u
#     log(lambda/2) - lambda * exp(-theta/2) - theta/2
#   }
#
#   # log(precision) for prior on beta
#   model_data$beta_logprec <- log(.05)
#
#
#   # Differenced matrices...
#   model_data$diffmat <- create_diff_matrix(model_data$control_days)
#   model_data$lambdainv <- create_full_dtcp_matrix(model_data$control_days)
#   model_data$A$exposure$Ad <- model_data$diffmat %*% model_data$A$exposure$A
#   model_data$Xd <- model_data$diffmat %*% model_data$X
#   # Random effect model specification data
#   model_data$modelspec <- model_data$A %>%
#     purrr::map("model") %>%
#     purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
#     purrr::reduce(bind_rows)
#
#   model_data$vectorofcolumnstoremove <- round(RW2BINS/2)
# }
