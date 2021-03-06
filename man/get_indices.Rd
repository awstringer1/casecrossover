% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05-posthoc-quantities.R
\name{get_indices}
\alias{get_indices}
\title{Helper function to return the correct indices for latent variables}
\usage{
get_indices(model_data, removezeros = TRUE)
}
\arguments{
\item{model_data}{ccmodeldata object as returned by model_setup().}

\item{removezeros}{Optional. Should the indices be returned with hard-zeroes included or removed from the
random effect vector? Depends on where this function is being called.}
}
\value{
A list of class ccindex containing vectors of indices for the linear and smooth terms.
}
\description{
The organization of the latent field is confusing. There
are random effects and fixed effects, plus the actual linear predictor values
(which you usually don't want).

This function accepts a ccmodeldata object and returns a list telling you
which indices in the vector of latent variables correspond to which model elements.

Details: the latent Gaussian variables have dimension Wd = Nd + M + p. Nd
is the total number of control days, sum(model_data$control_days), of all subjects.
M is the total number of random effects, i.e. unique values of smooth terms, in the
order in which they appear in the model. p is the number of regression coefficients.
To get the means/variances of the random effects, the correct call is
i = (model_data$Nd+1):(model_data$Nd+model_data$M). To get the regression coefficients,
the correct call is i = (model_data$Nd+model_data$M+1):(model_data$Nd+model_data$M+model_data$p).

All this is way too complicated for the user, and even for the developer, so here is a
function to parse the model data and return an object which enumerates which elements
of the (internal) latent field correspond to which quantities. This is mostly
used internally.
}
