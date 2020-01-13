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
#' see create_linear_constraints().
#'
#' thetaaccuracy: the accuracy of the grid used for numerical integration. If there are K hyperparameters
#' (dim(theta) = K), then thetaaccuracy^K regular gridpoints will be used. Default is 3, which is
#' probably too low.
#'
#' sparsetheta: logical. Use a sparse grid for the numerical integration of theta? This can dramatically
#' speed up computation time. Default FALSE.
#'
#' thetarange: numeric vector or matrix specifying the range of numerical integration for theta. If a vector,
#' give a lower and upper endpoint; if dim(theta) > 1 then these are recycled, creating a cube. If dim(theta) > 1
#' you can also specify a rectangle by giving a matrix with one column for the lower end points and one for the
#' upper.
#'
#' optcontrol: control parameters passed to trustOptim::trust.optim(). See their vignette, which is excellent.
#' However, the defaults here are probably fine for most problems.
#'
#' @param ... Arguments used to override the defaults output by cc_control().
#'
#' @examples
#' cc_control(smooth_prior = list(z = list(prior = "pc_prec",params = c(u = 3,a = 1))))
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stats formula terms
#' @importFrom methods as new
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
#' @param nm Optional; if you know the name of the prior you want, specify it.
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
#' @param verbose Logical. Should the function actually print error messages? TRUE by default; turned
#' off when used programatically (inside other functions).
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
validate_prior_distribution <- function(prior,verbose = TRUE) {
  # Check it's a list with the correct names
  islistwithcorrectnames <- is.list(prior) &
    length(names(prior)) == 3 &
    "name" %in% names(prior) &
    "params" %in% names(prior) &
    "call" %in% names(prior)
  if (!islistwithcorrectnames) {
    if (verbose) print("The correct format for a prior is a named list with elements 'name', 'params', and 'call'.")
    return(FALSE)
  }
  PRIOR_VALID <- TRUE
  supportedpriors <- supported_prior_distributions()
  # Check prior is supported
  if (prior$name %in% names(supportedpriors)) {
    if (verbose) print(stringr::str_c("Prior ",prior$name," is supported!"))
  } else {
    if (verbose) print(stringr::str_c("Prior ",prior$name," is not supported. Try running supported_prior_distributions() to see supported priors."))
    return(FALSE)
  }

  # Check parameters specified correctly
  in_user_not_in_supported <- dplyr::setdiff(names(prior$params),supportedpriors[[prior$name]]$params)
  in_supported_not_in_user <- dplyr::setdiff(supportedpriors[[prior$name]]$params,names(prior$params))
  if (length(in_user_not_in_supported) != 0) {
    if (verbose) print(stringr::str_c("The following parameters were specified by you but not supported: ",in_user_not_in_supported))
    PRIOR_VALID <- FALSE
  }
  if (length(in_supported_not_in_user) != 0) {
    if (verbose) print(stringr::str_c("The following parameters were not specified by you but are required: ",in_supported_not_in_user))
    PRIOR_VALID <- FALSE
  }

  if (PRIOR_VALID) {
    if (verbose) print(paste0("Your prior specification is valid!"))
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

#' Fortify a list of priors
#'
#' @description This function takes a list of priors, validates them using validate_prior_distribution(),
#' and stacks them together into one function that takes a vector argument, representing the joint prior
#' on all hyperparameters. This enforces prior independence between hyperparameters.
#'
#' The user shouldn't have to
#' use this function; it is exported in case advanced users want to define their own priors.
#'
#' @return A function which computes the joint log-prior of all hyperparameters, assuming prior
#' independence.
#'
#' @param priorlist A list of priors, each as created by e.g. pc_prior() or in the same format. Its
#' elements will be passed to validate_prior_distribution().
#'
#' @export
#'
fortify_priors <- function(priorlist) {
  # Validate them
  priorsvalid <- priorlist %>%
    purrr::map_lgl(~validate_prior_distribution(.x,verbose = FALSE)) %>%
    all()

  if (!priorsvalid) {
    stop("One or more priors failed validation. Check them each using validate_prior_distribution().")
  }

  # Fortify each of them
  flist <- list()
  for (prior in priorlist) {
    ff <- function(theta) {
      do.call(supported_prior_distributions(prior$name)$call,as.list(c(theta = theta,prior$params)))
    }
    environment(ff) <- new.env()
    environment(ff)$prior <- prior
    flist <- c(flist,ff)
  }
  flist <- flist %>% purrr::map(rlang::as_function) # Error checking

  # Return one function that computes the joint log-prior
  f <- function(theta) {
    # theta is a VECTOR now.
    # Check the number of priors provided is the same as the number of thetas
    if (length(theta) != length(flist)) stop(stringr::str_c(length(flist)," priors provided but ",length(theta)," hyperparameters detected."))

    flist %>%
      purrr::map2(theta,~.x(.y)) %>%
        purrr::reduce(sum)
  }
  rlang::as_function(f)
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
#' contain the sorted unique values of your covariate, the index of the zero value(s), the name of the variable,
#' and a list of sparseVector(s) implementing the constraint(s).
#'
#' @param u Covariate vector. You can pass it in raw (like data$u) or as a sorted vector of unique values.
#' @param whichzero Actual values of u for which you want the random effect to be zero.
#' @param nm The name of the covariate in your dataframe.
#'
#' @examples
#' temperature <- c(10,15,20,25,30,35,40)
#' create_linear_constraints(temperature,30,"temperature")
#' create_linear_constraints(temperature,c(30,35),"temperature")
#'
#' @export
create_linear_constraints <- function(u,whichzero,nm = "") {
  u <- sort(unique(u))
  if (!all(whichzero %in% u)) {
    notfound <- stringr::str_c(whichzero[which(!(whichzero %in% u))],collapse = ", ")
    stop(stringr::str_c("Value(s) ",notfound," not found in the covariate you supplied. Try rounding your data?"))
  }

  lu <- length(u)

  out <- list(list(
    u = u,
    whichzero = whichzero,
    constraint = purrr::map(whichzero,~sparseVector(x = 1,i = which(u == .x),length = lu))
  ))
  names(out) <- nm
  out
}

#' Validate your created linear constraints.
#'
#' @description This function checks the format of what you're inputting into cc_control() for
#' linear constraints. This is to help you format your model correctly and prevent downstream
#' errors, which can be harder to debug.
#'
#' @param constraints What you're planning on inputting to cc_control(linear_constraints = ...)
#' @param verbose Logical. Should debugging information be printed?
#'
#' @examples
#' validate_linear_constraints(create_linear_constraints(u = c(1,2,3),whichzero = 1)) # Valid
#' validate_linear_constraints(list(x = c(1,2,3))) # Not valid
#'
#' @export
validate_linear_constraints <- function(constraints,verbose = TRUE) {
  VALID_CONSTRAINT <- TRUE
  # Check if input is a list
  if (!is.list(constraints)) {
    if (verbose) print(stringr::str_c("Expected a list but you provided an object of class ",class(constraints),
                        ". Even if you have only one constraint, pass a list of constraints."))
  return(FALSE)
  }
  # If the input is only a single linear constraint, remind the user to use a named list.
  if (length(constraints) == 3 &
      "u" %in% names(constraints) &
      "constraint" %in% names(constraints) &
      "whichzero" %in% names(constraints)) {
    if (verbose) print("It appears you input only a single linear constraint object, not in a list. Remember to wrap your linear constraint in a named list")
    VALID_CONSTRAINT <- FALSE
  }

  # Check if the elements of the list each contain the correct items
  itemcheck <- constraints %>% purrr::map_lgl(~length(.x) == 3 &
                                                "u" %in% names(.x) &
                                                "constraint" %in% names(.x) &
                                                "whichzero" %in% names(.x))
  if (any(!itemcheck)) {
    if (verbose) print("The correct format for linear constraints is a list containing exactly elements 'u', 'whichzero', and 'constraint'.")
    return(FALSE)
  }

  # Check if all constraint vectors in the list are sparseVectors
  classcheck <- constraints %>%
    purrr::map("constraint") %>%
    purrr::map(1) %>%
    purrr::map_lgl(~inherits(.x,"sparseVector"))

  if (any(!classcheck)) {
    if (verbose) print("One or more vectors from the list you provided does not inherit from class sparseVector. Check your vectors!")
    return(FALSE)
  }

  # Check if each vector is the same length as the covariate
  lengthtester <- function(cnst) {
    lu <- length(cnst$u)
    cnstrlengths <- cnst$constraint %>% purrr::map_dbl(~length(.x))
    all(lu == cnstrlengths)
  }
  lengthcheck <- constraints %>% purrr::map_lgl(lengthtester)

  if (any(!lengthcheck)) {
    if (verbose) print("One or more of your constraints isn't the same length as the covariate")
    return(FALSE)
  }

  if (verbose) print("Your linear constraint appears valid!")
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
#' @param verbose Logical. Print progress and diagnostic information? Useful for debugging or
#' keeping up with what the function is doing.
#'
#' @return An object of class cc_modeldata, for passing into future functions.
#'
#' @examples
#' model_setup(y~x + strata(id),data.frame(y = c(0,1),x = c(1,2),id = c(1,1)))
#'
#' @export
# model_setup <- function(formula,data,control = cc_control(smooth_prior = pc_prior(3,.5))) {}
model_setup <- function(formula,data,control = cc_control(),verbose = FALSE) {
  model_data <- structure(list(), class = "cc_modeldata")

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
    for (nm in model_elements$smooth) {
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

  # Number of subjects
  n <- length(unique(data[model_elements$strata]))

  # Linear terms
  if (length(model_elements$linear) == 0) {
    model_data$X <- NULL
    model_data$p <- 0 # No linear terms
  } else {
    model_data$X <- Matrix::sparse.model.matrix(model_elements$linear_formula,data = data)
    model_data$p <- ncol(model_data$X)
    # Safety check: ncol(X) > 0.
    if (ncol(model_data$X) == 0) {
      model_data$X <- NULL
      model_data$p <- 0 # No linear terms
    }
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
    dplyr::rename(case_days = .data[[model_elements$response]])


  model_data$control_days <- control_days$control_days
  model_data$case_days <- case_days$case_days
  names(model_data$control_days) <- data[[model_elements$strata]] %>% unique() %>% sort()
  names(model_data$case_days) <- data[[model_elements$strata]] %>% unique() %>% sort()
  model_data$n <- length(model_data$case_days)
  model_data$Nd <- sum(model_data$control_days)
  model_data$Ne <- model_data$Nd + model_data$n

  # Priors.
  # Currently only pc.prec prior is implemented.
  # Prior u is P(sigma > u) = alpha.

  if (length(control$smooth_prior) != length(model_elements$smooth)) stop(stringr::str_c(length(control$smooth_prior)," priors provided for ",length(model_elements$smooth)," hyperparameters."))

  if (model_data$M > 0) {
    model_data$theta_logprior <- fortify_priors(control$smooth_prior)
  } else {
    model_data$theta_logprior <- function(theta) {force(theta); return(0)} # Placeholder
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

  # Check for linear constraints. Validate and then add them to the model data.
  # If there are linear constraints, take the first one for each variable and create the
  # model_data$vectorofcolumnstoremove element.
  if (length(model_elements$smooth) > 0) {
    if (length(control$linear_constraints) != length(model_elements$smooth)) {
      warning("Smooth terms, but no linear constraints, specified. You should add one or more constraints. See create_linear_constraints().")
    } else {
      # For each smooth term with a constraint, Take the first one and use it to set one to zero
      k <- 0
      s <- 1
      model_data$vectorofcolumnstoremove <- numeric()
      names(model_data$vectorofcolumnstoremove) <- character()
      for (nm in model_elements$smooth) {
        whichzero <- which(control$linear_constraints[[nm]]$whichzero[1] == control$linear_constraints[[nm]]$u)
        model_data$vectorofcolumnstoremove <- c(model_data$vectorofcolumnstoremove,whichzero + k)
        k <- k + length(control$linear_constraints[[nm]]$u)
        model_data$M <- model_data$M - 1
        names(model_data$vectorofcolumnstoremove)[s] <- nm
        s <- s + 1
      }
      # nm <- model_elements$smooth[1]
      # model_data$vectorofcolumnstoremove <- control$linear_constraints[[nm]]$whichzero[1]
      # Adjust M
      # model_data$M <- model_data$M - 1
      # Remove this constraint from the list of constraints
      # ACTUALLY: do I have to...? If I add it back in later I don't think it matters. We'll see.
    }
  }

  # Set final dimensions
  model_data$Wd <- model_data$M + model_data$p + model_data$Nd
  model_data$Wdf <- model_data$M + model_data$p + model_data$Ne

  # Add back the control list and model elements
  model_data$control <- control
  model_data$model_elements <- model_elements

  # Add the structure of the hessian
  model_data$hessian_structure <- hessian_log_likelihood_structure(W = rep(0,model_data$Wd),model_data = model_data)


  model_data
}
