### LATENT FIELD ###

# Functions to implement calculations relating to the latent variables, priors and posteriors.
# Most of these are not exported.

### Q matrix functions ###

#' Calculate the Q matrix for a model with only linear terms.
#'
#' @description Calculate the Q-matrix (joint prior precision matrix) of the latent Gaussian variables.
#' This is a function of theta. Functions are provided for Q-matrices with linear terms only,
#' RW2 terms only, and both.
#'
#' Note that the linear terms' prior covariance matrix is fixed (a true prior) and its log-precision
#' is specified directly in the model setup. Hence the Q_matrix_linear function does not take an argument
#' for this parameter.
#'
#' @return A sparse matrix inheriting from class CsparseMatrix containing the joint precision matrix of the
#' latent Gaussian variables.
#'
#' @param model_data A cc_modeldata object returned by model_setup()
#' @param tau The precision of the auxillary epsilon variables. Set it to something huge. Default exp(12).
#'
Q_matrix_linear <- function(model_data,tau=exp(12)) {
  Sbinv <- Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))
  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% model_data$Xd),
    cbind(-tau*t(model_data$Xd) %*% model_data$lambdainv,Sbinv + tau*crossprod(crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd))
  )
}

#' Calculate the precision matrix for ONE component of an RW2 model. These are used to construct the centre block
#' of the full Q matrix.
#'
#' @rdname Q_matrix_linear
#' @param theta Hyperparameter, scalar, the log-precision of the RW2 model.
#' @param covariate String naming the covariate from the data for which the precision matrix should
#' be calculated.
#'
Q_matrix_rw2_one_component <- function(theta,model_data,covariate) {
  u <- sort(unique(model_data$A[[covariate]]$u))

  ul <- length(u)
  du <- diff(u)

  H <- bandSparse(n = ul,
                  diagonals = list(
                    c(1/du[-(ul-1)],0),
                    c(0,-(1/du[-(ul-1)] + 1/du[-1]),0),
                    c(0,1/du[-1])
                  ),
                  k = c(-1,0,1))

  AA <- Diagonal(n = ul,x = c(2/du[1],2/(du[-(ul-1)] + du[-1]),2/du[(ul-1)]))

  exp(theta) * forceSymmetric(crossprod(H,crossprod(AA,H)))
}

#' Calculate the Q matrix for a model with only RW2 terms
#'
#' @rdname Q_matrix_linear
#'

Q_matrix_rw2 <- function(theta,model_data,tau = exp(12)) {
  # Figure out how many rw2 components there are
  if (is.null(model_data$A)) stop("no rw2 components in model")

  # whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% purrr::pull(covariate)
  whichrw2 <- model_data$model_elements$smooth

  howmanyrw2 <- length(whichrw2)

  if (length(theta) != howmanyrw2) stop(stringr::str_c("Detected ",howmanyrw2," RW2 components in model, but length(theta) = ",length(theta),". Make sure your model is correctly specified."))


  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta[.y],model_data,covariate = .x)) %>%
    bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% purrr::map("Ad") %>% purrr::reduce(cbind)

  # APPLY ZERO CONSTRAINTS
  # Directly apply constraint that U_t = 0 by deleting the t^th column of A,
  # and the t^th row and column of Suinv

  if (!(0 %in% model_data$vectorofcolumnstoremove)) {
    Ad <- Ad[ ,-model_data$vectorofcolumnstoremove]
    Suinv <- Suinv[-model_data$vectorofcolumnstoremove,-model_data$vectorofcolumnstoremove]
  }
  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% Ad),
    cbind(-tau*t(Ad)%*% model_data$lambdainv,Suinv + tau * crossprod(Ad,crossprod(model_data$lambdainv,Ad)))
  )
}


#' Q matrix for model containing both linear and rw2 terms
#'
#' @rdname Q_matrix_linear
#'
Q_matrix_both <- function(theta,model_data,tau = exp(12)) {
  if (is.null(model_data$A)) stop("no rw2 components in model")

  # whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% purrr::pull(covariate)
  whichrw2 <- model_data$model_elements$smooth

  howmanyrw2 <- length(whichrw2)

  if (length(theta) != howmanyrw2) stop(stringr::str_c("Detected ",howmanyrw2," RW2 components in model, but length(theta) = ",length(theta),". Make sure your model is correctly specified."))

  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta[.y],model_data,covariate = .x)) %>%
    bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% purrr::map("Ad") %>% purrr::reduce(cbind)

  # APPLY ZERO CONSTRAINTS
  # Directly apply constraint that U_t = 0 by deleting the t^th column of A,
  # and the t^th row and column of Suinv

  if (!(0 %in% model_data$vectorofcolumnstoremove)) {
    Ad <- Ad[ ,-model_data$vectorofcolumnstoremove]
    Suinv <- Suinv[-model_data$vectorofcolumnstoremove,-model_data$vectorofcolumnstoremove]
  }

  # Linear
  if (model_data$p == 0) stop("no linear terms in model")
  Sbinv <- Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))

  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Ad
  Q13 <- -tau*model_data$lambdainv %*% model_data$Xd
  Q22 <- Suinv + tau * crossprod(Ad,crossprod(model_data$lambdainv,Ad))
  Q23 <- tau * crossprod(Ad,crossprod(model_data$lambdainv,model_data$Xd))
  Q33 <- Sbinv + tau*crossprod(crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd)
  rbind(
    cbind(Q11,Q12,Q13),
    cbind(t(Q12),Q22,Q23),
    cbind(t(Q13),t(Q23),Q33)
  )
}

#' Compute the Q matrix for any model.
#'
#' @rdname Q_matrix_linear
#'
Q_matrix <- function(theta,model_data,tau = exp(12)) {
  if (is.null(model_data$A)) {
    if (model_data$p == 0) stop("both X and A are null...")
    mat <- Q_matrix_linear(model_data,tau)
  }
  else {
    if (model_data$p > 0) {
      mat <- Q_matrix_both(theta,model_data,tau)
    } else {
      mat <- Q_matrix_rw2(theta,model_data,tau)
    }
  }

  forceSymmetric(mat)
}

### Priors and Posteriors ###

#' log-prior of W
#'
#' @description These functions calculate priors and posteriors of the latent Gaussian variables.
#'
#' @param W Value of W = (delta,gamma,beta) to calculate the prior/posterior at. The order is as
#' it appears in the appendix of the paper- (delta_1_1,...,delta_n_Jn,gamma_1,...,gamma_M,beta_1,...,beta_p).
#' The order of gamma and beta is the same as they are listed in model_data$model_elements.
#' @param model_data A ccmodeldata object returned by model_setup().
#' @param theta Value of the hyperparameter vector theta. The Q matrix depends on this. Can leave as NULL
#' if you're passing in the Q matrix
#' @param Q The Q-matrix as returned by Q_matrix(model). If not provided, will be calculated.
#'
logprior_W <- function(W,model_data,theta = NULL,Q = NULL) {
  if (is.null(Q)) {
    if (is.null(theta)) stop("If Q not provided, theta must be provided.")
    Q <- Q_matrix(theta,model_data,tau = model_data$control$tau)
  }
  # Check dimensions and throw an informative error
  if (!(all(length(W) == dim(Q)))) stop(stringr::str_c("Length of W is ",length(W)," but dimension of Q is ",dim(Q),". Check your model data."))
  -as.numeric((1/2)*crossprod(W,crossprod(Q,W)))
}

#' Log-posterior of W, given theta and y
#'
#' @rdname logprior_W
#'
log_posterior_W <- function(W,theta,model_data,Q = NULL) {
  logprior_W(W,model_data,theta,Q) + log_likelihood(W,model_data)
}



#' Gradient of log-posterior of W, given theta and y
#'
#' @rdname logprior_W
#'
grad_log_posterior_W <- function(W,theta,model_data,Q = NULL) {
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data,tau = model_data$control$tau)
  }
  -as.numeric(crossprod(Q,W)) + grad_log_likelihood(W,model_data)
}

#' Hessian of log-posterior of W, given theta and y
#'
#' @rdname logprior_W
#' @param structure Optional. Pass in the sparse structure of the hessian of the log-likelihood
#' to save computing time. See hessian_log_likelihood_structure().
#'
hessian_log_posterior_W <- function(W,theta = NULL,Q = NULL,model_data,structure = NULL) {
  if (is.null(theta) & is.null(Q)) stop("One of Q or theta must be provided")
  if (is.null(Q)) Q <- Q_matrix(theta,model_data,tau = model_data$control$tau)
  -(Q + hessian_log_likelihood(W,model_data,structure))
}

# Q <- Q_matrix(.1,model_data)
# microbenchmark::microbenchmark(
#   hessian_log_posterior_W(W,Q = Q,model_data = model_data),
#   hessian_log_posterior_W(W,theta = .1,model_data = model_data)
# )
# Unit: milliseconds
# expr      min       lq     mean median       uq      max neval
# hessian_log_posterior_W(W, Q = Q, model_data = model_data) 863.1047 1033.178 1200.637 1198.622 1317.136 2038.524   100
# hessian_log_posterior_W(W, theta = 0.1, model_data = model_data) 881.3422 1022.531 1196.749 1239.993 1326.682 1895.039   100
# So providing Q or not makes basically no difference, it saves like 4 milliseconds on average.
# But providing the structure...
# tmpstructure <- hessian_log_likelihood_structure(W,model_data)
# tmpstructure <- list(i = tmpstructure@i,p = tmpstructure@p)
# microbenchmark::microbenchmark(
#   hessian_log_posterior_W(W,Q = Q,model_data = model_data),
#   hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = tmpstructure)
# )
# Unit: milliseconds
# expr      min lq      mean    median        uq       max neval
# hessian_log_posterior_W(W, Q = Q, model_data = model_data) 865.1010 1064.3718 1213.0569 1223.1004 1344.2056 2047.0038   100
# hessian_log_posterior_W(W, Q = Q, model_data = model_data, structure = tmpstructure) 120.0503 148.8404  167.4932  166.9949  182.7649  280.0576   100
# Order of magnitude speedup, woohoo!





### HYPERPARAMETERS ###

#' Log-posterior approximation for theta and sigma
#'
#' @description Compute the log-posterior for theta, the log-precision, and sigma, the standard deviation, of the RW2 model components.
#'
#' @param theta Vector of log-precisions at which to evaluate the log-posterior. theta = -2*log(sigma).
#' @param W The mode of the log-posterior of W|theta,y for this theta. Computed outside and passed in.
#' @param model_data a ccmodeldata object created using model_setup()
#' @param hessian_structure Optional, sparse structure of the Hessian created using hessian_log_likelihood_structure()
#' @param Q Optional, Q matrix evaluated at this theta. Will be created if not supplied.
#'
log_posterior_theta <- function(theta,W,model_data,hessian_structure = NULL,Q = NULL) {
  # W is the mode of log_posterior_W(theta)
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data,tau = model_data$control$tau)
  }
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(W,model_data)
  }
  Q_p_C <- -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = hessian_structure)

  term1 <- log_likelihood(W,model_data)
  dt <- determinant(Q,logarithm = TRUE) # For this, we DO need the determinant
  term2_det <- (1/2) * as.numeric(dt$modulus)
  term2 <- logprior_W(W,model_data,theta) # Doesn't contain the determinant
  term3 <- model_data$theta_logprior(theta)
  qcdet <- determinant(Q_p_C,logarithm = TRUE)
  term4 <- -(1/2) * as.numeric(qcdet$modulus)
  as.numeric(term1 + term2_det + term2 + term3 + term4)
}

#' Log-posterior approximation for theta and sigma
#'
#' @rdname log_posterior_theta
#' @param sigma Vector of standard deviations at which to evaluate the log-posterior. sigma = exp(-.5 * theta)
#'
log_posterior_sigma <- function(sigma,W,model_data,hessian_structure = NULL,Q = NULL) {
  length(sigma)* log(2) - sum(log(sigma)) +
    log_posterior_theta(theta = -2 * log(sigma),
                         W = W,
                         model_data = model_data,
                         hessian_structure = hessian_structure,
                         Q = Q)
}
