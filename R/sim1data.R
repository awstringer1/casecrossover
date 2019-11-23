#' A simulated dataset for fitting case-crossover models.
#'
#' @description These data were simulated from a multinomial distribution of the form described in the
#' case crossover paper, with true additive predictor eta = 3 *(exposure^2 - .5^2). It is good for
#' testing the fitting of case crossover models and learning the API for prior specification and linear
#' constraints.
#'
#' @format A tibble with 3,596 rows and 6 columns:
#'
#' \describe{
#'   \item{exposure}{The covariate used for the simulations}
#'   \item{eta}{The true value of the additive predictor, for comparison}
#'   \item{prob}{The corresponding true probability of case == 1}
#'   \item{subject}{Subject id, to be used as the strata variable}
#'   \item{exposure_binned}{The exposure covariate binned into 50 bins, for fitting RW2 models}
#' }
#'
"sim1data"
