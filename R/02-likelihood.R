### Likelihood ###
# This script contains functions for computing the likelihood for the casecrossover model
# This includes the gradient and hessian and any helper functions
# Most (all?) of these functions are not exported and are hidden from users, though
# they are still documented.

#' Prepare the data for computing the log-likelihood
#'
#' @description Take in the parameter vector and the model_data and prepare the data
#' for computing the log-likelihood. The log-likelihood is not permutation-invariant,
#' which is a fancy math way of saying "it matters what order the parameters are in".
#' Inside model_data, you have a named vector of control_days containing each subject's
#' id and the number of control days they have in the data. The data is sorted in
#' ascending order of id, in blocks with all control days followed by the case day.
#' This function splits the parameter vector W into a list, with one item per id,
#' containing the ordered elements of W corresponding to each subjects' control and
#' case days.
#'
#' @param W Parameter vector. First n elements are eta, then Gamma and beta.
#' @param model_data A list of class "cc_modeldata" as returned by model_setup().
#'
#' @return A list with n items (where n is the number of subjects) containing vectors
#' of the parameters corresponding to each subject's control and case days.
#'
prep_data_for_log_lik <- function(W,model_data) {
  # Get delta
  delta <- W[1:model_data$Nd]

  # Split delta into a list containing the values for each control day for each subject
  split(delta,as.numeric(rep(names(model_data$control_days),model_data$control_days)))
}

#' Get a vector of case day weights from the model_data
#'
#' @description Get a vector of case day weights, repeated once for every parameter
#' ascribed to each subject. Used inside gradient and hessian calculations as weights.
#'
#' @return A numeric vector of case day weights
#'
#' @param model_data A list of class "cc_modeldata" as returned by model_setup().
#'
get_casedays <- function(model_data) unname(rep(model_data$case_days,model_data$control_days))

#' Compute the log-likelihood for the model at a given parameter configuration
#'
#' @description Compute the model log-likelihood for a given set of parameters.
#'
#' @param W Parameter vector. First n elements are eta, then Gamma and beta.
#' @param model_data A list of class "cc_modeldata" as returned by model_setup().
#'
#' @return A number, the log-likelihood.
#'
#'
log_likelihood <- function(W,model_data) {
  # Split the parameter vector
  deltasplit <- prep_data_for_log_lik(W,model_data)

  # Helper function to compute the log(1 + sum(exp(-delta))) for each person.
  # Use logsumexp to avoid underflow in the individual delta terms inside the sum
  # Then exponentiate the result- if the WHOLE SUM is so close to -Inf that exponentiating
  # it gives nearly zero, it's not a problem, since 1 + x basically equals 1 in that case
  # anyways. Underflow is only a problem when adding the exp(-delta_it) terms.
  compute_loglik_term <- function(deltavec) {
    # deltavec is a vector of deltas for each person
    # Note we are adding both minus signs in here.
    -log(1 + exp(matrixStats::logSumExp(-deltavec)))
  }

  # Now apply that function to the whole list and sum the result to get the answer
  llterms <- purrr::map(deltasplit,compute_loglik_term) %>% purrr::reduce(c)
  sum(llterms * model_data$case_days)
}


#' Gadient of the case-crossover log-likelihood.
#'
#' @description Compute the gradient of the log likelihood with respect to W = (delta,gamma,beta).
#' The gamma and beta parts are zero.
#'
#' @return A numeric vector containing the gradient. It's not stored as a sparseVector,
#' because it's (mathematically) not sparse. I guess the end is all zeroes though, so maybe
#' we should change this...?
#'
#' @inheritParams log_likelihood
#'
grad_log_likelihood <- function(W,model_data) {
  # Split the parameter vector
  deltasplit <- prep_data_for_log_lik(W,model_data)
  casedays <- get_casedays(model_data)
  # Helper to compute the gradient term
  compute_gradient_term <- function(deltavec) {
    # Compute the denominator
    denom <- 1 + exp(matrixStats::logSumExp(-deltavec))
    # Return exp(-deltavec)/denom
    exp(-deltavec)/denom
  }
  # The gradient is the concatenation of all these (vector) terms,
  # plus zeroes on the end to make it as long as W
  gradient_front <- purrr::map(deltasplit,compute_gradient_term) %>% purrr::reduce(c)
  gradient_back <- rep(0,length(W) - length(gradient_front))
  c(gradient_front * casedays,gradient_back)
}

#' Compute the negated Hessian of the log-likelihood
#'
#' @description Function to implement the (sparse) NEGATED hessian of the case-crossover log-likelihood
#' This is MINUS the second derivative. Because that's what's used in the paper.
#' So remember when optimizing: provide the NEGATIVE of log_likelihood and gradient_...,
#' but provide the POSITIVE of this function.
#'
#' Two functions are written. hessian_log_likelihood_x() provides the ELEMENTS of the hessian only
#' This is very fast (134 milliseconds on average for the air pollution example)
#' hessian_log_likelihood_structure() brute-force computes the hessian using bdiag,
#' which is an order of magnitude slower (around a second on average). But the structure
#' never changes, so we only need to compute this once.
#' hessian_log_likelihood() computes the hessian, but you can supply it the structure and it
#' will basically just wrap hessian_log_likelihood_x(). Much faster.
#'
#' @inheritParams log_likelihood
#'
#' @return A sparse matrix inheriting from class CsparseMatrix containing the negated
#' Hessian of the log-likelihood.
#'

hessian_log_likelihood_structure <- function(W,model_data) {
  # Arguments: see log_likelihood
  # Returns: a list with i and p for the sparse hessian. NOTE: currently does it
  # for dgCMatrix- i.e. not symmetric. I couldn't figure out how to do it for symmetric.
  # Still fast and memory efficient.
  # Returns the structure as (i,p), in column-compressed format. To get (i,j) triplet
  # format do triplet = TRUE

  # Split the parameter vector
  deltasplit <- prep_data_for_log_lik(W,model_data)

  # Helper function to create the blocks
  # Takes in a vector of deltasplit and returns a block
  # These blocks are dense.
  compute_hessian_block <- function(deltavec) {
    denom <- 1 + exp(matrixStats::logSumExp(-deltavec))
    expdeltavecscaled <- exp(-deltavec)/denom

    if (length(expdeltavecscaled) == 1) return(expdeltavecscaled * (1 - expdeltavecscaled))

    diag(expdeltavecscaled) - expdeltavecscaled %o% expdeltavecscaled
  }

  # Create a list of blocks
  blocklist <- purrr::map(deltasplit,compute_hessian_block)

  # Then create a big zero matrix of dimension equal to length(W) - sum(length(deltasplit))
  blocklist <- c(blocklist,Matrix(0,model_data$Wd - model_data$Nd,model_data$Wd - model_data$Nd))

  # Final result is a big block diagonal matrix consisting of the created blocks
  # concatenated with the zero block. The Matrix::bdiag function is REALLY slow for
  # the case of many small dense matrices, as we have here (see docs).

  # UPDATE: changed to symmetric structure.
  structure <- bdiag(blocklist)
  # structure <- forceSymmetric(bdiag(blocklist))
  return(list(i = structure@i,p = structure@p))
}

#' @rdname hessian_log_likelihood_structure
#'

hessian_log_likelihood_x <- function(W,model_data) {
  deltasplit <- prep_data_for_log_lik(W,model_data)
  casedays <- model_data$case_days # Note: this is intentionally different than for the gradient.

  # Helper function to create the blocks
  # Takes in a vector of deltasplit and returns a block
  # These blocks are dense.
  compute_hessian_block <- function(deltavec) {
    denom <- 1 + exp(matrixStats::logSumExp(-deltavec))
    dveclen <- length(deltavec)
    expdeltavecscaled <- exp(-deltavec)/denom

    if (length(expdeltavecscaled) == 1) return(expdeltavecscaled * (1 - expdeltavecscaled))

    outmat <- diag(expdeltavecscaled) - expdeltavecscaled %o% expdeltavecscaled
    outmat[upper.tri(outmat,diag = TRUE)]
  }

  purrr::map2(deltasplit,casedays,
              ~compute_hessian_block(.x) * .y) %>% purrr::reduce(c)

}

#' @rdname hessian_log_likelihood_structure
#' @param structure Optional. A list returned by hessian_log_likelihood_structure()
#' which contains the sparse structure of the hessian. Computing this is the slow
#' part, but only needs to be done once.

hessian_log_likelihood <- function(W,model_data,structure=NULL) {
  # structure is a list containing elements i and p
  # If not provided, will call hessian_log_likelihood_structure()
  if (is.null(structure)) {
    # Column-compressed format
    structure <- hessian_log_likelihood_structure(W,model_data)
  }

  out <- new("dgCMatrix")
  out@i <- structure$i
  out@p <- structure$p
  out@Dim <- as.integer(rep(length(structure$p) - 1,2))
  out@x <- numeric(length(out@i))
  out <- as(out,'symmetricMatrix')
  out@x <- hessian_log_likelihood_x(W,model_data)
  return(as(out,'dgCMatrix'))
}
