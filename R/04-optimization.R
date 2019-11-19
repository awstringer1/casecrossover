### Optimization
# This script contains the functions for computing the conditional mode W-hat
# for a user-specified grid of theta values.

#' Convert the output of one of my optimization functions into a tibble for easy analysis.
#'
#' Not exported; for development use only.
#'
#' @param lst A list formatted according to the output of optimize_latentfield_trustoptim().
#'
optlist_to_tibble <- function(lst) {
  # If lst is a single element of opt, wrap it in a list so the below works
  if ("optimizer" %in% names(lst)) lst <- list(lst)
  dplyr::tibble(
    theta = purrr::map(lst,"theta"),
    optimizer = purrr::map(lst,"optimizer"),
    starting = purrr::map(lst,"starting"),
    solution = purrr::map(lst,"solution"),
    function_value = purrr::map_dbl(lst,"function_value"),
    iterations = purrr::map_dbl(lst,"iterations")
  )
}

#' Optimize the conditional posterior of the latent Gaussian variables for one hyperparameter configuration
#'
#' @description Compute the conditional mode W-hat for one value of the hyperparameters theta. Uses
#' trust region methods implemented in trustOptim::trust.optim. Makes use of the sparsity of the Hessian.
#'
#' @param theta Vector, representing one configuration of hyperparameters (optimize_latentfield_trustoptim). List of vectors, each
#' representing one configuration of hyperparameters (optimize_all_thetas_parallel).
#' @param model_data ccmodeldata object as output by model_setup().
#' @param hessian_structure Optional, provide the sparsity structure of the Hessian. This doesn't change with different thetas, so when
#' doing all hyperparameter configurations, only need to compute this once. It will be computed if not provided.
#' @param Q Optional, provide pre-computed Q matrix for this theta. Will be computed if not provided.
#' @param optcontrol Optional. Override the casecrossover default control parameters for trust.optim(). See cc_control()$opt_control.
#'
#' @export
#'
optimize_latentfield_trustoptim <- function(theta,model_data,hessian_structure = NULL,Q = NULL,optcontrol = NULL) {

  # Zero is the prior mean. It's a reasonable starting value.
  # My experience has been this particular optimization problem is pretty robust to this.
  startingvals <- rep(0,model_data$Wd)

  # Set up the functions. Negated, functions of W only with other arguments fixed.
  optfunc <- function(W) -1 * log_posterior_W(W,theta,model_data,Q)
  optfuncgrad <- function(W) -1 * grad_log_posterior_W(W,theta,model_data,Q)

  # Get the Q matrix if not provided
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  # Get the hessian structure, if not provided
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(startingvals,model_data)
  }

  optfunchess <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = hessian_structure)


  # Set default control arguments if not provided.
  # Note: this is really for custom development use only. In the main casecrossover() function, these are set using
  # the cc_control() function.
  if (is.null(optcontrol)) {
    optcontrol <- cc_control()$opt_control
  }

  # Perform the optimization

  opt <- trustOptim::trust.optim(
    x = startingvals,
    fn = optfunc,
    gr = optfuncgrad,
    hs = optfunchess,
    method = "Sparse",
    control = optcontrol
  )

  # Return a custom-formatted list with optimization results.
  # This is done to enable easy stacking of one of these per theta value in a dataframe.

  list(
    optimizer = list(
      name = "trust_optim",
      hessian = opt$hessian,
      trust_radius = opt$trust.radius,
      norm_grad = sqrt(sum(opt$gradient^2)),
      scaled_norm_grad = sqrt(sum(opt$gradient^2)) / sqrt(model_data$Wd),
      status = opt$status,
      nnz = opt$nnz
    ),
    theta = theta,
    starting = startingvals,
    solution = opt$solution,
    function_value = opt$fval,
    iterations = opt$iterations
  )
}

#' Find the conditional model for a user-defined list of theta values
#'
#' @description This is a wrapper for optimize_latentfield_trustoptim() that computes the
#' conditional mode in parallel over a user-specified grid of theta values.
#'
#' Parallelization only works on linux and mac; the "parallel" package needs to work. If you're on windows, you'll just
#' have to wait a little longer (long enough to go to the store and buy a serious computer maybe?).
#'
#' @rdname optimize_latentfield_trustoptim
#' @param thetagrid A NIGrid object which contains the grid of theta values to be optimized for.
#' @param doparallel Logical. Use parallel::mclapply (TRUE) or lapply (FALSE) for execution? Former executes in parallel, latter in series. Default TRUE.
#'
#' @export
#'
optimize_all_thetas_parallel <- function(thetagrid,model_data,hessian_structure = NULL,optcontrol = NULL,doparallel = TRUE) {

  # Check thetagrid is formatted correctly
  if (!inherits(thetagrid,"NIGrid")) stop("theta should be a NIGrid object returned by mvQuad::createNIgrid()")
  # Create the theta list
  theta <- split(mvQuad::getNodes(thetagrid),rep(1:nrow(mvQuad::getNodes(thetagrid)),ncol(mvQuad::getNodes(thetagrid))))

  if (!all(purrr::map_lgl(theta,is.numeric))) stop("theta should be a list of numeric vectors")
  thetalengths <- purrr::map_dbl(theta,length)
  if (length(unique(thetalengths)) != 1) stop("Make sure all thetas are the same length.")

  if (is.null(optcontrol)) {
    optcontrol <- cc_control()$opt_control
  }

  do_opt <- function(theta) {
    cat("Optimizing with theta =",theta,"\n")
    optimize_latentfield_trustoptim(theta = theta,
                                    model_data = model_data,
                                    hessian_structure = hessian_structure,
                                    optcontrol = optcontrol)
  }

  cat("Performing optimization, start time: ", format(Sys.time(),"%H:%M:%S"),"\n")

  tm <- proc.time()
  if (doparallel) {
    opt <- parallel::mclapply(theta,do_opt)
  } else {
    opt <- lapply(theta,do_opt)
  }
  opt_time <- proc.time() - tm
  opt_time <- unname(opt_time["elapsed"])
  cat("Time taken for optimization with ",length(theta)," values of theta:",opt_time,"seconds.\n")

  out <- optlist_to_tibble(opt)
  attr(out,"thetagrid") <- thetagrid
  out
}
