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
#' linear_constraints: Specify linear constraints on your smooth model. You should
#' input a named list of vectors. The name of the list is the name of the covariate that you
#' want the constraint on. The length of the vector should be equal to the number of unique
#' values of this covariate. The vector should satisfy a^T u = 0, where a is the vector and
#' u is a SORTED (ascending) vector of unique covariate values. For help creating this list,
#' see create_linear_combinations().
#'
#' @param ... Arguments used to override the defaults output by cc_control().
#'
#' @examples
#' cc_control(smooth_prior = list(z = list(prior = "pc_prec",params = c(u = 3,a = 1))))
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stats formula terms
#' @import Matrix
#' @export
cc_control <- function(...) {
  outargs <- cc_default_control()
  if (length(list(...)) == 0) return(outargs)
  for (nm in names(list(...))) {
    outargs[[nm]] <- list(...)[[nm]]
  }
  outargs
}

#' Supported prior distributions.
#'
#' @description This function returns a list of all supported prior distributions,
#' including their names, parameters, and actual function call.
#'
#' @param nm Optional; if you know the name of the
#' @examples
#' supported_prior_distributions()
#' supported_prior_distributions("pc.prec")
#'
#' @export
supported_prior_distributions <- function(nm = "") {
  supported_priors <- list(
    pc.prec = list(
      name = "pc.prec",
      params = c("u","alpha"),
      call = prior_calls[["pc.prec"]]
    )
  )

  if (nm != "") return(supported_priors[[nm]])
  supported_priors
}

#' Validate a specified prior distribution.
#'
#' @description Pass in your control parameters, and this function will tell you whether you specified
#' your priors correctly. It just checks supported_prior_distributions() for you; so also, see that function.
#' You won't get a formal error if your prior is determined to be incorrectly specified; you'll get a
#' message that tries to help you fix it.
#'
#' @param prior The smooth_prior element of a control list as returned by cc_control().
#'
#' @examples
#' mycontrol <- cc_control(smooth_prior = list(name = "pc.prec",params = c(u = 3,alpha = .5)))
#' mycontrol_missinginfo <- cc_control(smooth_prior = list(name = "pc.prec"))
#' mycontrol_badprior <- cc_control(smooth_prior = list(name = "foo"))
#'
#' validate_prior_distribution(mycontrol)
#' validate_prior_distribution(mycontrol_missinginfo)
#' validate_prior_distribution(mycontrol_badprior)
#'
#' @export
validate_prior_distribution <- function(prior) {
  # Check it's a list with the correct names
  islistwithcorrectnames <- is.list(prior) & length(names(prior)) == 3 & "name" %in% names(prior) & "params" %in% names(prior) & "call" %in% names(prior)
  if (!islistwithcorrectnames) {
    print("The correct format for a prior is a named list with elements 'name', 'params', and 'call'.")
    return(FALSE)
  }
  PRIOR_VALID <- TRUE
  supportedpriors <- supported_prior_distributions()
  # Check prior is supported
  if (prior$name %in% names(supportedpriors)) {
    print(stringr::str_c("Prior ",prior$name," is supported!"))
  } else {
    print(stringr::str_c("Prior ",prior$name," is not supported. Try running supported_prior_distributions() to see supported priors."))
    return(FALSE)
  }

  # Check parameters specified correctly
  in_user_not_in_supported <- dplyr::setdiff(names(prior$params),supportedpriors[[prior$name]]$params)
  in_supported_not_in_user <- dplyr::setdiff(supportedpriors[[prior$name]]$params,names(prior$params))
  if (length(in_user_not_in_supported) != 0) {
    print(stringr::str_c("The following parameters were specified by you but not supported: ",in_user_not_in_supported))
    PRIOR_VALID <- FALSE
  }
  if (length(in_supported_not_in_user) != 0) {
    print(stringr::str_c("The following parameters were not specified by you but are required: ",in_supported_not_in_user))
    PRIOR_VALID <- FALSE
  }

  if (PRIOR_VALID) {
    print(paste0("Your prior specification is valid!"))
    return(TRUE)
  }

  print(paste0("Your prior specification is not valid. See above messages for how to fix it."))
  return(FALSE)

}

#' The PC Prior for the smoothing log precision.
#'
#' @description The PC prior for the smoothing log precision is described in the package vignette.
#' This is a convenience function for specifying this prior in the casecrossover software.
#'
#' The prior is specified as P(sigma > u) = 1 - alpha.
#'
#' @param u The "u" in P(sigma > u) = 1 - alpha.
#' @param alpha The "alpha" in P(sigma > u) = 1 - alpha.
#'
#' @examples
#' pc_prior(u = 3,alpha = .5)
#'
#' @export
pc_prior <- function(u,alpha) {
  list(
    name = "pc.prec",
    params = c(u = u,alpha = alpha),
    call = supported_prior_distributions("pc.prec")$call
  )
}

#' Create a linear combination vector for a single-element-zero constraint on a smooth term
#'
#' @description The most common and easiest to interpret constraint on a random walk model for
#' a smooth term in a case-crossover model is simply to set a single element of the random effect
#' to zero. This means that effects are interpreted as relative effects compared to this "reference" value.
#' The reason this is the default constraint, as opposed to the usual "sum-to-zero" constraint used in
#' random walk models, is because the main advantage of the "sum-to-zero" constraint is that it is
#' orthogonal to the intercept in the model. In a case-crossover model, the intercept is not estimable,
#' and while it's still totally possible to use a sum-to-zero constraint, it becomes less clear how
#' to interpret it.
#'
#' The function returns a named list of sparseVectors suitable for input into cc_control(). Specifically, the list items
#' are themselves lists, containing the sorted unique values of your covariate, and a sparseVector
#' implementing the constraint.
#'
#' @param u Covariate vector. You can pass it in raw (like data$u) or as a sorted vector of unique values.
#' @param whichzero Actual values of u for which you want the random effect to be zero.
#'
#' @examples
#' temperature <- c(10,15,20,25,30,35,40)
#' create_linear_constraints(temperature,30)
#' create_linear_constraints(temperature,c(30,35))
#'
#' @export
create_linear_constraints <- function(u,whichzero) {
  u <- sort(unique(u))
  if (!all(whichzero %in% u)) {
    notfound <- stringr::str_c(whichzero[which(!(whichzero %in% u))],collapse = ", ")
    stop(stringr::str_c("Value(s) ",notfound," not found in the covariate you supplied. Try rounding your data?"))
  }

  lu <- length(u)

  purrr::map(whichzero,~list(u = u,constraint = sparseVector(x = 1,i = which(u == .x),length = lu)))
}

#' Validate your created linear constraints.
#'
#' @description This function checks the format of what you're inputting into cc_control() for
#' linear constraints. This is to help you format your model correctly and prevent downstream
#' errors, which can be harder to debug.
#'
#' @param constraints What you're planning on inputting to cc_control(linear_constraints = ...)
#'
#' @examples
#' validate_linear_constraints(create_linear_constraints(u = c(1,2,3),whichzero = 1)) # Valid
#' validate_linear_constraints(list(x = c(1,2,3))) # Not valid
#'
#' @export
validate_linear_constraints <- function(constraints) {
  VALID_CONSTRAINT <- TRUE
  # Check if input is a list
  if (!is.list(constraints)) {
    print(stringr::str_c("Expected a list but you provided an object of class ",class(constraints),
                        ". Even if you have only one constraint, pass a list of sparseVectors."))
  return(FALSE)
  }

  # Check if the elements of the list each contain the correct items
  itemcheck <- constraints %>% purrr::map_lgl(~length(.x) == 2 & "u" %in% names(.x) & "constraint" %in% names(.x))
  if (any(!itemcheck)) {
    print("The correct format for linear constraints is a list of lists, each sub-list containing exactly elements u and constraint.")
    return(FALSE)
  }

  # Check if all constraint vectors in the list are sparseVectors
  classcheck <- constraints %>%
    purrr::map("constraint") %>%
    purrr::map_lgl(~inherits(.x,"sparseVector"))

  if (any(!classcheck)) {
    print("One or more vectors from the list you provided does not inherit from class sparseVector. Check your vectors!")
    return(FALSE)
  }

  # Check if each vector is the same length as the covariate
  lengthcheck <- constraints %>%
    purrr::map_lgl(~length(.x$u) == length(.x$constraint))

  if (any(!lengthcheck)) {
    print("One or more of your constraints isn't the same length as the covariate")
    return(FALSE)
  }

  print("Your linear constraint appears valid!")
  return(TRUE)
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
#' @param control A list containing control parameters. See cc_control().
#' @examples
#' model_setup(y~x + strata(id),data.frame(y = c(0,1),x = c(1,2),id = c(1,1)))
#'
#' @export
# model_setup <- function(formula,data,control = cc_control(smooth_prior = pc_prior(3,.5))) {}
model_setup <- function(formula,data,control = cc_control()) {
  # TODO: implement this.
  # TODO: document this.
  model_data <- structure(list(), class = "ccmodeldata")

  # Parse the formula
  model_elements <- parse_formula(formula)
  # Check that the smooth and strata terms exist in the data
  # The linear terms will be passed to model.matrix, which has its own
  # error checking.
  extra_model_vars <- model_elements[c("smooth","strata")] %>% purrr::reduce(c)
  if (!all(extra_model_vars %in% colnames(data))) {
    missing_vars <- extra_model_vars[!(extra_model_vars %in% colnames(data))]
    stop(paste0("The following variables were provided in the model formula but not in the data: ",
                stringr::str_c(missing_vars,collapse = ", ")))
  }

  # Create the smooth terms- design matrix
  Alist <- list()
  if (length(model_elements$smooth) > 0) {
    for (nm in names(model_elements$smooth)) {
      Alist[[nm]] <- create_alist_element(data[[nm]])
    }
  }
  model_data$A <- Alist
  if (length(Alist) == 0) {
    model_data$A <- NULL
    model_data$M <- 0 # No smooth terms
  } else {
    model_data$M <- Alist %>% purrr::map("A") %>% purrr::map(ncol) %>% purrr::reduce(sum)
  }
  # TODO: when implementing the linear constraint part, where one element of the
  # precision matrix is set to 0, make sure to make the appropriate reduction in M

  # Number of subjects
  n <- length(unique(data[model_elements$strata]))

  # Linear terms
  if (length(model_elements$linear) == 0) {
    model_data$X <- NULL
    model_data$p <- 0 # No linear terms
  } else {
    model_data$X <- Matrix::sparse.model.matrix(model_elements$linear_formula,data = data)
    model_data$p <- ncol(model_data$X)
  }
  # Safety check: ncol(X) > 0.
  if (ncol(model_data$X) == 0) {
    model_data$X <- NULL
    model_data$p <- 0 # No linear terms
  }


  # Create the vector of control days
  # A named vector where the names are the subject ids and the values are the number
  # of control days that each has in the data

  control_days <- data %>%
    dplyr::arrange(.data[[model_elements$strata]],.data[[model_elements$response]]) %>%
    dplyr::filter(.data[[model_elements$response]] == 0) %>%
    dplyr::group_by(.data[[model_elements$strata]]) %>%
    dplyr::summarize(control_days = n())

  # Create the vector of case days
  # A named vector where the names are the subject ids and the values are the number
  # of control days that each has in the data

  case_days <- data %>%
    dplyr::arrange(.data[[model_elements$strata]],.data[[model_elements$response]]) %>%
    dplyr::filter(.data[[model_elements$response]] != 0) %>%
    dplyr::group_by(.data[[model_elements$strata]]) %>%
    dplyr::summarize(case_days = n())


  model_data$control_days <- control_days$control_days
  model_data$case_days <- case_days$case_days
  names(model_data$control_days) <- data[[model_elements$strata]] %>% unique() %>% sort()
  names(model_data$case_days) <- data[[model_elements$strata]] %>% unique() %>% sort()
  model_data$Nd <- sum(model_data$control_days)
  model_data$Ne <- model_data$Nd + model_data$n
  model_data$Wd <- model_data$M + model_data$p + model_data$Nd
  model_data$Wdf <- model_data$M + model_data$p + model_data$Ne

  # Priors.
  # Currently only pc prec prior is implemented.
  # Prior u is P(sigma > u) = alpha.
  #
  # First, validate the user's prior is correctly specified
  priorvalid <- TRUE
  if (model_data$M > 0) {
    validate_prior_distribution(control$smooth_prior)
    if (!priorvalid) stop("Please specify a valid prior for your smooth terms.")
    # If it's valid, grab the actual function call and set the parameters
    model_data$theta_logprior <- function(theta) {
      callargs <- as.list(c(theta = theta,control$smooth_prior$params))
      do.call(supported_prior_distributions(control$smooth_prior$name)$call,callargs)
    }
  } else {
    model_data$theta_logprior <- function() {return(0)} # Placeholder
  }

  # log(precision) for prior on beta. Specified in control
  model_data$beta_logprec <- control$beta_prior_logprec


  # Differenced matrices:
  model_data$diffmat <- create_diff_matrix(model_data$control_days)
  model_data$lambdainv <- create_full_dtcp_matrix(model_data$control_days)
  if (model_data$M > 0) {
    for (nm in names(model_data$A)) {
      model_data$A[[nm]]$Ad <- model_data$diffmat %*% model_data$A[[nm]]$A
    }
  }
  if (model_data$p > 0) {
    model_data$Xd <- model_data$diffmat %*% model_data$X
  }

  # Random effect model specification data
  # model_data$modelspec <- NULL
  # if (model_data$M > 0) {
  #   model_data$modelspec <- model_data$A %>%
  #     purrr::map("model") %>%
  #     purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
  #     purrr::reduce(dplyr::bind_rows)
  # }

  # TODO: linear constraints. The first one should be taken out explicitly then put back
  # in. The others are corrected after. This behaviour does not need to be exposed to the user.
  # model_data$vectorofcolumnstoremove <- round(RW2BINS/2)

  # Check for linear constraints. Validate and then add them to the model data.
  # If there are linear constraints, take the first one for each variable and create the
  # model_data$vectorofcolumnstoremove element.


}
