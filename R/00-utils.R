# This file contains miscillaneous top-level non-exported functions
# for the casecrossover package


# INTERNAL: parse formula
parse_formula <- function(ff) {
  # Parse the formula ff into linear and smooth terms and a strata
  # Linear terms will be passed to model.matrix()
  # Smooth terms will be handled in a proprietary manner
  # The strata argument must be provided.

  # Grab the RHS elements of the formula
  ff_elements <- attributes(terms(ff))$term.labels
  # Get the names of the variables. This strips off the s() wrappers
  # and includes the respose.
  ff_variables <- all.vars(ff)

  # Split the names into response, linear, and smooth terms.
  # Throw an error if any terms are not supported.
  response <- ff_variables[1] # Response is always the first
  # Match smooth terms using regex
  # Will only match terms that look like s(variable). Extra stuff in the s()
  # will cause an error. Keep formulas clean!
  smooth <- stringr::str_extract(ff_elements,"^s\\(\\w+\\)$")
  smooth <- smooth[!is.na(smooth)]
  # Remove it from ff
  if (length(smooth) > 0) {
    for (tm in smooth) ff <- update(ff,formula(stringr::str_c(".~.-",tm)))
  }
  # Strip off the s() part
  smooth <- stringr::str_remove(smooth,"^s\\(") %>% stringr::str_remove("\\)$")
  # Match strata terms in same way as smooth
  strata <- stringr::str_extract(ff_elements,"^strata\\(\\w+\\)$")
  strata <- strata[!is.na(strata)]
  # Remove it from ff
  if (length(strata) > 0) {
    for (tm in strata) ff <- update(ff,formula(stringr::str_c(".~.-",tm)))
  }
  # Strip off the s() part
  strata <- stringr::str_remove(strata,"^strata\\(") %>% stringr::str_remove("\\)$")

  # All the terms that are left are treated as linear and passed to model.matrix
  # within the model_setup function. Any errors will come from that call.

  # If there is not exactly one strata term, throw an error
  if (length(strata) == 0) stop("No strata variable provided.")
  if (length(strata) > 1) {
    stop(paste0("You must provide exactly one strata variable. The following strata variables were parsed: ",
                stringr::str_c(strata,collapse=", ")))
  }

  # Return a named list of vectors with the names of each type of term
  list(
    linear = attributes(terms(ff))$term.labels,
    linear_formula = ff,
    smooth = smooth,
    strata = strata,
    response = response
  )
}

ff <- y ~ x + s(z) + strata(id) + a*b + log(d) + I(xx^2) + poly(zz,2)
parse_formula(ff)


# INTERNAL: prescribe default control arguments
cc_default_control <- function() {
  list(
    smooth_prior = list(),
    opt_control = list(
      prec = 1e-06,
      stop.trust.radius = 1e-03,
      report.freq = 1,
      report.level = 4,
      start.trust.radius = 10,
      contract.threshold = .25,
      contract.factor = .1,
      expand.factor = 5,
      preconditioner = 1,
      trust.iter = 200000,
      cg.tol = 1e-06,
      maxit = 1000
    )
  )
}

# INTERNAL: create a single element of the random effects structure list.
# This includes the design matrix, the actual values of the covariate,
# and information about the model. Currently only rw2 is supported.
create_alist_element <- function(u,constraint = NULL) {
  # u: covariate. NOT sorted and MAY contain ties/repeated values, in general.
  # constraint: vector containing values of u for which random effect U should be
  # constrained to be zero.
  lu <- length(u)
  A <- Diagonal(n = lu)[match(u,unique(u)),order(unique(u))]
  model <- "rw2"
  constrzero <- NULL
  if (!is.null(constraint)) {
    constrzero <- match(constraint,sort(unique(u)))
    if (any(is.na(constrzero))) warning(paste0("no match found for constraint: ",constraint[which(is.na(constrzero))]))
  }

  list(u = u,A = A,model = model,constrzero = constrzero)
}
