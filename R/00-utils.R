# This file contains miscillaneous top-level non-exported functions
# for the casecrossover package


# INTERNAL: parse formula
parse_formula <- function(ff) {
  # Parse the formula ff into linear and smooth terms and a strata
  # Linear terms will be passed to model.matrix()
  # Smooth terms will be handled in a proprietary manner
  # The strata argument must be provided.
  fullff <- ff # Save the input formula for later printing.

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

  # Finally, remove the intercept. If a -1 was already in the formula, this operation does nothing,
  # so it's safe.
  ff <- update(ff,.~.-1)

  # Return a named list of vectors with the names of each type of term
  list(
    # linear = attributes(terms(ff))$term.labels,
    linear = setdiff(all.vars(ff),response),
    linear_formula = ff,
    smooth = smooth,
    strata = strata,
    response = response,
    call = fullff
  )
}

# INTERNAL: get the degree of a polynomial from a formula
get_polynomial_degree <- function(ff) {
  ffvars <- all.vars(ff)[-1]
  ffattr <- attributes(terms(ff))$term.labels

  varnameregex <- "A-Za-z0-9_."

  degree_1 <- stringr::str_extract(ffattr,stringr::str_c("^[",varnameregex,"]+$"))
  degree_1 <- degree_1[!is.na(degree_1)]

  degree_more_than_1 <- stringr::str_extract(ffattr,stringr::str_c("^poly\\([",varnameregex,"=\\s\\,\\)]*"))
  degree_more_than_1 <- degree_more_than_1[!is.na(degree_more_than_1)]

  # Get the names
  deg_mt1_names <- stringr::str_extract(degree_more_than_1,stringr::str_c("^poly\\([",varnameregex,"]+")) %>%
    stringr::str_remove("^poly\\(")

  deg_mt1_degrees <- stringr::str_extract(degree_more_than_1,stringr::str_c("^poly\\([",varnameregex,"]+\\,\\s?[A-Za-z\\s=]*[0-9]")) %>%
    stringr::str_remove(stringr::str_c("^poly\\([",varnameregex,"]+\\,\\s?[A-Za-z\\s=]*")) %>%
      as.numeric()

  out <- c(rep(1,length(degree_1)),deg_mt1_degrees)
  names(out) <- c(degree_1,deg_mt1_names)
  out
}


# INTERNAL: prescribe default control arguments
cc_default_control <- function() {
  list(
    smooth_prior = list(),
    linear_constraints = list(),
    doparallel = TRUE,
    thetaaccuracy = 3,
    sparsetheta = FALSE,
    thetarange = c(-1,1),
    beta_prior_logprec = log(1/10),
    tau = exp(12),
    opt_control = list(
      prec = 1e-08,
      stop.trust.radius = sqrt(.Machine$double.eps),
      report.freq = 1,
      report.level = 4,
      start.trust.radius = 1000,
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
  A <- Matrix::Diagonal(n = lu)[match(u,unique(u)),order(unique(u))]
  model <- "rw2"
  constrzero <- NULL
  if (!is.null(constraint)) {
    constrzero <- match(constraint,sort(unique(u)))
    if (any(is.na(constrzero))) warning(paste0("no match found for constraint: ",constraint[which(is.na(constrzero))]))
  }

  list(u = u,A = A,model = model,constrzero = constrzero)
}

# INTERNAL: Functions to create the differencing matrix
# First create a differencing matrix of supplied dimension
create_single_diff_matrix <- function(d) {
  # d: ROW dimension. Result will be d x (d+1). d is the number of control days in our application
  cbind(Matrix::Diagonal(d,-1),Matrix::Matrix(1,d,1))
}
# Now create a function to create the whole big differencing matrix, given a vector of
# number of control days
create_diff_matrix <- function(control_days) {
  purrr::map(control_days,create_single_diff_matrix) %>% bdiag()
}

# dd <- create_diff_matrix(model_data$control_days) # dim = 6,835; size = 220 kB

# Now the inverse of the diff mat transpose crossproduct (dtcp)
# Function to create a single diffmat inverse
create_single_dtcp_inverse <- function(d) {
  Diagonal(d,1) - Matrix(1/(d+1),d,d)
}
# Function to create the whole block diagonal thing
# Will be slow because of bdiag(), but only have to do it once
create_full_dtcp_matrix <- function(control_days) {
  purrr::map(control_days,create_single_dtcp_inverse) %>% bdiag()
}



# INTERNAL: list of priors and their calls.
# The actual functions are defined in 00-prior-distributions.R
prior_calls <- list(
  pc.prec = pcprec
)

# INTERNAL: functions for numerical integration of log-posterior
normalize_log_posterior_single <- function(pp,tt) {
  df <- dplyr::tibble(pp = pp,tt = tt) %>% dplyr::arrange(tt)
  tt <- df$tt
  pp <- df$pp

  lp <- length(pp)
  matrixStats::logSumExp(c(
    matrixStats::logSumExp(pp[-1] + log(diff(tt)) + log(1/2)),
    matrixStats::logSumExp(pp[-lp] + log(diff(tt)) + log(1/2))
  ))
}

normalize_log_posterior_multiple <- function(pp,tt) {
  # Multidimensional quadrature
  # Only works for evenly-spaced grids
  M <- length(tt) # Number of grid points
  K <- length(tt[[1]]) # Dimension of integration

  # This only works for uniform grids:
  df <- tt[[2]][1] - tt[[1]][1]
  ww <- log(df) + c(log(1/2),rep(0,M-2),log(1/2))
  matrixStats::logSumExp(K*ww + pp)
}

# f <- function(x,y,z) mvtnorm::dmvnorm(x = c(x,y,z),log = TRUE)
#
# ttg <- expand.grid(seq(-5,5,by=.5),seq(-5,5,by=.5),seq(-5,5,by=.5))
# pp <- ttg %>% rowwise() %>% mutate(fx = f(Var1,Var2,Var3)) %>% pull(fx)
# tt <- list()
# for (i in 1:nrow(ttg)) {
#   tt[[i]] <- as.numeric(ttg[i, ])
# }
#
# normalize_log_posterior_multiple(pp,tt) # Should be zero, i.e. log(1)

# INTERNAL: take a vector of values, a vector of indices, and return a vector
# containing the values of the first with zeroes spliced in at the indices
# specified by the second
stitch_zero_vector <- function(x,z) {
  if (any(z > length(x) + 1)) stop("You are asking for zeroes at positions not covered by x")
  y <- as(x,"sparseVector")
  # This drops zero values in x automatically- add them back in
  y@i <- 1:length(x)
  y@x <- x
  for (j in z) {
    y@i[y@i >= j] <- y@i[y@i >= j] + 1
  }
  as.numeric(y)
}

# INTERNAL: compute quantiles given a vector of grid points and logged density evaluations
compute_quantiles <- function(tt,pp,origgrid) {
  pp <- exp(pp)

}





