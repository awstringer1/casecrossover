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
model_setup <- function(formula,data) {
  # TODO: implement this.
  # TODO: document this.
  return(0)
}
