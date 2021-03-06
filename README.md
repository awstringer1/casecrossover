<!-- badges: start -->
[![Travis build status](https://travis-ci.org/awstringer1/casecrossover.svg?branch=master)](https://travis-ci.org/awstringer1/casecrossover)
<!-- badges: end -->

# casecrossover
R package supporting case crossover models as described in Approximate Bayesian Inference for Case Crossover Models (2020), accepted to Biometrics.

This package is no longer under development. I am focussing on new methods and software, notably the `aghq` package and methods ([github](https://github.com/awstringer1/aghq), [vignette](https://arxiv.org/abs/2101.04468), [methods/theory preprint](https://arxiv.org/abs/2102.06801)). Expect a `casecrossover2` package which uses this underlying framework to appear in the future!

# Description

This R package implements approximate Bayesian inference for case crossover models. The methodology and interface supports linear
and non-linear covariate effects and can handle datasets with hundreds of thousands of observations and multiple control days 
for each subject. 

# Installation

The package is not yet on CRAN. Install using `devtools::install_github("awstringer1/casecrossover")`.

The first version of this package is nearly completed. For details on what can and can't be done currently, check out the
vignette, either by installing the package and typing `vignette("casecrossover")`, or from the `vignettes/` folder in this
repository.
