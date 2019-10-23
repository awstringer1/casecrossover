<!-- badges: start -->
[![Travis build status](https://travis-ci.org/awstringer1/casecrossover.svg?branch=master)](https://travis-ci.org/awstringer1/casecrossover)
<!-- badges: end -->

# casecrossover
R package supporting case crossover models as described in ...

# Description

This R package implements approximate Bayesian inference for case crossover models. The methodology and interface supports linear
and smooth covariate effects and can handle datasets with hundreds of thousands of observations and multiple control days 
for each subject. You can compute marginal posterior distributions for effects, as well as 
global confidence envelopes for the smooth functions, and you can get samples from the marginal and joint posteriors.

This package is in development. A vignette will be made available once a stable version is released, with real-data examples and code that reproduces the simulations in the paper (the real data in the paper is private and cannot be distributed).
