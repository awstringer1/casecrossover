context("Latent Variables - Priors and Posteriors")
library(casecrossover)
library(dplyr) # For tests
source("prep-sample-data.R")

### Q matrix functions ###

# Linear Q matrix
test_that("Linear-only Q matrix is created as expected",{
  # Function runs as expected
  expect_s4_class(Q_matrix_linear(model_data1),"CsparseMatrix")
  expect_true(isSymmetric(Q_matrix_linear(model_data1)))
  expect_s4_class(Q_matrix_linear(model_data2),"CsparseMatrix")
  expect_true(isSymmetric(Q_matrix_linear(model_data2)))

  # Function returns a matrix of the correct dimension
  expect_equal(dim(Q_matrix_linear(model_data1)),rep(model_data1$Wd,2))
  expect_equal(dim(Q_matrix_linear(model_data2)),rep(model_data2$Wd,2))
})

# RW2 Q matrix - centre block for one covariate
test_that("RW2 single component creation works as expected",{
  # Function runs as expected
  expect_s4_class(Q_matrix_rw2_one_component(0,model_data7,"x"),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2_one_component(0,model_data8,"x"),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2_one_component(0,model_data7,"x2"),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2_one_component(0,model_data8,"x2"),"CsparseMatrix")

  expect_true(isSymmetric(Q_matrix_rw2_one_component(0,model_data7,"x")))
  expect_true(isSymmetric(Q_matrix_rw2_one_component(0,model_data8,"x")))
  expect_true(isSymmetric(Q_matrix_rw2_one_component(0,model_data7,"x2")))
  expect_true(isSymmetric(Q_matrix_rw2_one_component(0,model_data8,"x2")))

  # Function returns a matrix of the correct dimension
  expect_equal(dim(Q_matrix_rw2_one_component(0,model_data7,"x")),rep(length(unique(sampledata$x)),2))
  expect_equal(dim(Q_matrix_rw2_one_component(0,model_data8,"x")),rep(length(unique(sampledata$x)),2))
  expect_equal(dim(Q_matrix_rw2_one_component(0,model_data7,"x2")),rep(length(unique(sampledata$x2)),2))
  expect_equal(dim(Q_matrix_rw2_one_component(0,model_data8,"x2")),rep(length(unique(sampledata$x2)),2))
})

# RW2 full Q matrix
test_that("RW2 full Q matrix computed as expected",{
  # Function runs as expected
  expect_error(Q_matrix_rw2(0,model_data1)) # Linear terms only
  expect_error(Q_matrix_rw2(0,model_data2))
  expect_s4_class(Q_matrix_rw2(0,model_data3),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2(0,model_data4),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2(0,model_data5),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2(0,model_data6),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2(c(0,0),model_data7),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2(c(0,0),model_data8),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2(100,model_data3),"CsparseMatrix")
  expect_s4_class(Q_matrix_rw2(-100,model_data3),"CsparseMatrix")

  expect_true(isSymmetric(Q_matrix_rw2(0,model_data3)))
  expect_true(isSymmetric(Q_matrix_rw2(0,model_data4)))
  expect_true(isSymmetric(Q_matrix_rw2(0,model_data5)))
  expect_true(isSymmetric(Q_matrix_rw2(0,model_data6)))
  expect_true(isSymmetric(Q_matrix_rw2(c(0,0),model_data7)))
  expect_true(isSymmetric(Q_matrix_rw2(c(0,0),model_data8)))

  # Function returns a matrix of the correct dimension
  expect_equal(dim(Q_matrix_rw2(0,model_data3)),rep(model_data3$Wd,2))
  expect_equal(dim(Q_matrix_rw2(0,model_data4)),rep(model_data4$Wd,2))
  expect_equal(dim(Q_matrix_rw2(c(0,0),model_data7)),rep(model_data7$Wd,2))
  expect_equal(dim(Q_matrix_rw2(c(0,0),model_data8)),rep(model_data8$Wd,2))
})

# RW2 + Linear, full Q matrix
test_that("RW2 + Linear full Q matrix computed as expected",{
  # Function runs as expected
  expect_error(Q_matrix_both(0,model_data1)) # Linear terms only
  expect_error(Q_matrix_both(0,model_data2))
  expect_error(Q_matrix_both(0,model_data3)) # Smooth terms only
  expect_error(Q_matrix_both(0,model_data4))

  expect_s4_class(Q_matrix_both(0,model_data5),"CsparseMatrix")
  expect_s4_class(Q_matrix_both(0,model_data6),"CsparseMatrix")

  expect_true(isSymmetric(Q_matrix_both(0,model_data5)))
  expect_true(isSymmetric(Q_matrix_both(0,model_data6)))

  # Function returns a matrix of the correct dimension
  expect_equal(dim(Q_matrix_both(0,model_data5)),rep(model_data5$Wd,2))
  expect_equal(dim(Q_matrix_both(0,model_data6)),rep(model_data6$Wd,2))
})

# Q matrix construction, overall
test_that("Overall Q matrix construction works",{
  # Function runs as expected
  expect_s4_class(Q_matrix(0,model_data1),"CsparseMatrix")
  expect_s4_class(Q_matrix(0,model_data2),"CsparseMatrix")
  expect_s4_class(Q_matrix(0,model_data3),"CsparseMatrix")
  expect_s4_class(Q_matrix(0,model_data4),"CsparseMatrix")
  expect_s4_class(Q_matrix(0,model_data5),"CsparseMatrix")
  expect_s4_class(Q_matrix(0,model_data6),"CsparseMatrix")
  expect_s4_class(Q_matrix(c(0,0),model_data7),"CsparseMatrix")
  expect_s4_class(Q_matrix(c(0,0),model_data8),"CsparseMatrix")
  expect_s4_class(Q_matrix(100,model_data3),"CsparseMatrix")
  expect_s4_class(Q_matrix(-100,model_data3),"CsparseMatrix")

  expect_true(isSymmetric(Q_matrix(0,model_data1)))
  expect_true(isSymmetric(Q_matrix(0,model_data2)))
  expect_true(isSymmetric(Q_matrix(0,model_data3)))
  expect_true(isSymmetric(Q_matrix(0,model_data4)))
  expect_true(isSymmetric(Q_matrix(0,model_data5)))
  expect_true(isSymmetric(Q_matrix(0,model_data6)))
  expect_true(isSymmetric(Q_matrix(c(0,0),model_data7)))
  expect_true(isSymmetric(Q_matrix(c(0,0),model_data8)))

  # Function returns a matrix of the correct dimension
  expect_equal(dim(Q_matrix(0,model_data1)),rep(model_data1$Wd,2))
  expect_equal(dim(Q_matrix(0,model_data2)),rep(model_data2$Wd,2))
  expect_equal(dim(Q_matrix(0,model_data3)),rep(model_data3$Wd,2))
  expect_equal(dim(Q_matrix(0,model_data4)),rep(model_data4$Wd,2))
  expect_equal(dim(Q_matrix(0,model_data5)),rep(model_data5$Wd,2))
  expect_equal(dim(Q_matrix(0,model_data6)),rep(model_data6$Wd,2))
  expect_equal(dim(Q_matrix(c(0,0),model_data7)),rep(model_data7$Wd,2))
  expect_equal(dim(Q_matrix(c(0,0),model_data8)),rep(model_data8$Wd,2))
})


### Priors and posteriors ###

W11 <- rep(0,model_data1$Wd)
W12 <- c(rep(0,model_data1$Wd - 1),100000) # The last value, for beta, shouldn't affect the likelihood.
W13 <- 1:model_data1$Wd # Nonzero parameter values
W14 <- rnorm(model_data1$Wd) # Random values- ensure this is robust!

W21 <- rep(0,model_data2$Wd)
W22 <- c(rep(0,model_data2$Wd - 1),100000) # The last value, for beta, shouldn't affect the likelihood.
W23 <- 1:model_data2$Wd # Nonzero parameter values
W24 <- rnorm(model_data2$Wd) # Random values- ensure this is robust!

W31 <- rep(0,model_data3$Wd)
W32 <- c(rep(0,model_data3$Wd - 1),100000) # The last value, for beta, shouldn't affect the likelihood.
W33 <- 1:model_data3$Wd # Nonzero parameter values
W34 <- rnorm(model_data3$Wd) # Random values- ensure this is robust!

W51 <- rep(0,model_data5$Wd)
W52 <- c(rep(0,model_data5$Wd - 1),100000) # The last value, for beta, shouldn't affect the likelihood.
W53 <- 1:model_data5$Wd # Nonzero parameter values
W54 <- rnorm(model_data5$Wd) # Random values- ensure this is robust!

W71 <- rep(0,model_data7$Wd)
W72 <- c(rep(0,model_data7$Wd - 1),100000) # The last value, for beta, shouldn't affect the likelihood.
W73 <- 1:model_data7$Wd # Nonzero parameter values
W74 <- rnorm(model_data7$Wd) # Random values- ensure this is robust!


Q11 <- Q_matrix(0,model_data1)
Q12 <- Q_matrix(-5,model_data1)

Q21 <- Q_matrix(0,model_data2)
Q22 <- Q_matrix(-5,model_data2)

Q31 <- Q_matrix(0,model_data3)
Q32 <- Q_matrix(-5,model_data3)

Q51 <- Q_matrix(0,model_data5)
Q52 <- Q_matrix(-5,model_data5)

Q71 <- Q_matrix(c(0,0),model_data7)
Q72 <- Q_matrix(c(-5,-5),model_data7)



# Logprior for W
test_that("Log-prior for W calculated correctly",{
  # Simplest model
  expect_equal(logprior_W(W11,model_data1,0),0)
  expect_equal(logprior_W(W12,model_data1,0),as.numeric(-.5*t(W12)%*%Q11%*%W12))
  expect_equal(logprior_W(W13,model_data1,0),as.numeric(-.5*t(W13)%*%Q11%*%W13))
  expect_equal(logprior_W(W14,model_data1,0),as.numeric(-.5*t(W14)%*%Q11%*%W14))

  expect_equal(logprior_W(W11,model_data1,-5),0)
  expect_equal(logprior_W(W12,model_data1,-5),as.numeric(-.5*t(W12)%*%Q12%*%W12))
  expect_equal(logprior_W(W13,model_data1,-5),as.numeric(-.5*t(W13)%*%Q12%*%W13))
  expect_equal(logprior_W(W14,model_data1,-5),as.numeric(-.5*t(W14)%*%Q12%*%W14))

  # Make sure the different case days don't mess up the prior (prior doesn't involve the data)
  expect_equal(logprior_W(W11,model_data2,0),0)
  expect_equal(logprior_W(W12,model_data2,0),as.numeric(-.5*t(W12)%*%Q21%*%W12))
  expect_equal(logprior_W(W13,model_data2,0),as.numeric(-.5*t(W13)%*%Q21%*%W13))
  expect_equal(logprior_W(W14,model_data2,0),as.numeric(-.5*t(W14)%*%Q21%*%W14))

  expect_equal(logprior_W(W11,model_data2,-5),0)
  expect_equal(logprior_W(W12,model_data2,-5),as.numeric(-.5*t(W12)%*%Q22%*%W12))
  expect_equal(logprior_W(W13,model_data2,-5),as.numeric(-.5*t(W13)%*%Q22%*%W13))
  expect_equal(logprior_W(W14,model_data2,-5),as.numeric(-.5*t(W14)%*%Q22%*%W14))

  # Model with smooth terms
  expect_equal(logprior_W(W31,model_data3,0),0)
  expect_equal(logprior_W(W32,model_data3,0),as.numeric(-.5*t(W32)%*%Q31%*%W32))
  expect_equal(logprior_W(W33,model_data3,0),as.numeric(-.5*t(W33)%*%Q31%*%W33))
  expect_equal(logprior_W(W34,model_data3,0),as.numeric(-.5*t(W34)%*%Q31%*%W34))

  expect_equal(logprior_W(W31,model_data3,-5),0)
  expect_equal(logprior_W(W32,model_data3,-5),as.numeric(-.5*t(W32)%*%Q32%*%W32))
  expect_equal(logprior_W(W33,model_data3,-5),as.numeric(-.5*t(W33)%*%Q32%*%W33))
  expect_equal(logprior_W(W34,model_data3,-5),as.numeric(-.5*t(W34)%*%Q32%*%W34))

  # Model with linear and smooth terms
  expect_equal(logprior_W(W51,model_data5,0),0)
  expect_equal(logprior_W(W52,model_data5,0),as.numeric(-.5*t(W52)%*%Q51%*%W52))
  expect_equal(logprior_W(W53,model_data5,0),as.numeric(-.5*t(W53)%*%Q51%*%W53))
  expect_equal(logprior_W(W54,model_data5,0),as.numeric(-.5*t(W54)%*%Q51%*%W54))

  expect_equal(logprior_W(W51,model_data5,-5),0)
  expect_equal(logprior_W(W52,model_data5,-5),as.numeric(-.5*t(W52)%*%Q52%*%W52))
  expect_equal(logprior_W(W53,model_data5,-5),as.numeric(-.5*t(W53)%*%Q52%*%W53))
  expect_equal(logprior_W(W54,model_data5,-5),as.numeric(-.5*t(W54)%*%Q52%*%W54))

  # Model with two smooth terms
  # Make sure to specify two thetas!
  expect_error(logprior_W(W71,model_data7,0)) # Only one theta
  expect_equal(logprior_W(W71,model_data7,c(0,0)),0)
  expect_equal(logprior_W(W72,model_data7,c(0,0)),as.numeric(-.5*t(W72)%*%Q71%*%W72))
  expect_equal(logprior_W(W73,model_data7,c(0,0)),as.numeric(-.5*t(W73)%*%Q71%*%W73))
  expect_equal(logprior_W(W74,model_data7,c(0,0)),as.numeric(-.5*t(W74)%*%Q71%*%W74))

  expect_equal(logprior_W(W71,model_data7,c(-5,-5)),0)
  expect_equal(logprior_W(W72,model_data7,c(-5,-5)),as.numeric(-.5*t(W72)%*%Q72%*%W72))
  expect_equal(logprior_W(W73,model_data7,c(-5,-5)),as.numeric(-.5*t(W73)%*%Q72%*%W73))
  expect_equal(logprior_W(W74,model_data7,c(-5,-5)),as.numeric(-.5*t(W74)%*%Q72%*%W74))
})

# Log posterior
# These functions just return the sum of two functions that have already been unit tested
# So only need a few tests.
# Since at W = 0 the log prior is 0, I just test equality of the posterior with the likelihood
# (and gradient/hessian etc) at W = 0. This basically just ensures the functions return results
# without errors. The likelihood functions are all extensively tested.
# For the Hessian I have to actually add the Q matrix since it's not zero.

test_that("Log posterior and gradient and hessian do not throw errros",{
  expect_equal(log_posterior_W(W11,0,model_data1),log_likelihood(W11,model_data1))
  expect_equal(log_posterior_W(W21,0,model_data2),log_likelihood(W21,model_data2))
  expect_equal(log_posterior_W(W31,0,model_data3),log_likelihood(W31,model_data3))
  expect_equal(log_posterior_W(W51,0,model_data5),log_likelihood(W51,model_data5))
  expect_equal(log_posterior_W(W71,c(0,0),model_data7),log_likelihood(W71,model_data7))

  expect_equal(grad_log_posterior_W(W11,0,model_data1),grad_log_likelihood(W11,model_data1))
  expect_equal(grad_log_posterior_W(W21,0,model_data2),grad_log_likelihood(W21,model_data2))
  expect_equal(grad_log_posterior_W(W31,0,model_data3),grad_log_likelihood(W31,model_data3))
  expect_equal(grad_log_posterior_W(W51,0,model_data5),grad_log_likelihood(W51,model_data5))
  expect_equal(grad_log_posterior_W(W71,c(0,0),model_data7),grad_log_likelihood(W71,model_data7))

  expect_equal(hessian_log_posterior_W(W11,0,model_data = model_data1),-(Q_matrix(0,model_data1) + hessian_log_likelihood(W11,model_data1)))
  expect_equal(hessian_log_posterior_W(W21,0,model_data = model_data2),-(Q_matrix(0,model_data2) + hessian_log_likelihood(W21,model_data2)))
  expect_equal(hessian_log_posterior_W(W31,0,model_data = model_data3),-(Q_matrix(0,model_data3) + hessian_log_likelihood(W31,model_data3)))
  expect_equal(hessian_log_posterior_W(W51,0,model_data = model_data5),-(Q_matrix(0,model_data5) + hessian_log_likelihood(W51,model_data5)))
  expect_equal(hessian_log_posterior_W(W71,c(0,0),model_data = model_data7),-(Q_matrix(c(0,0),model_data7) + hessian_log_likelihood(W71,model_data7)))

  expect_s4_class(hessian_log_posterior_W(W11,0,model_data = model_data1),'CsparseMatrix')
  expect_s4_class(hessian_log_posterior_W(W21,0,model_data = model_data2),'CsparseMatrix')
  expect_s4_class(hessian_log_posterior_W(W31,0,model_data = model_data3),'CsparseMatrix')
  expect_s4_class(hessian_log_posterior_W(W51,0,model_data = model_data5),'CsparseMatrix')
  expect_s4_class(hessian_log_posterior_W(W71,c(0,0),model_data = model_data7),'CsparseMatrix')
})

# Log-posterior for theta. Will recreate (inefficiently) using the mvtnorm::dmvnorm() function.
# Better version of the dmvnorm function that doesn't require inverting the precision matrix.
better_dmvnorm <- function(W,mean,Q) {
  p <- length(W)
  as.numeric(-(p/2)*log(2*pi) + (1/2)*determinant(Q,logarithm = TRUE)$modulus - (1/2)*crossprod(W-mean,crossprod(Q,W-mean)))
}
lpt_inefficient <- function(theta,W,model_data) {
  Q <- Q_matrix(theta,model_data) # Densify for dmvnorm (rolls eyes)
  C <- hessian_log_likelihood(W,model_data)
  model_data$theta_logprior(theta) +
    better_dmvnorm(W = W,mean = rep(0,length(W)),Q = Q) + # Mean 0- prior
    log_likelihood(W,model_data) -
    better_dmvnorm(W = W,mean = W,Q = Q + C)
}



test_that("Log posterior for theta is computed correctly",{
  expect_equal(log_posterior_theta(0,W11,model_data1),lpt_inefficient(0,W11,model_data1))
  expect_equal(log_posterior_theta(-5,W11,model_data1),lpt_inefficient(-5,W11,model_data1))
  expect_equal(log_posterior_theta(0,W12,model_data1),lpt_inefficient(0,W12,model_data1))
  expect_equal(log_posterior_theta(-5,W12,model_data1),lpt_inefficient(-5,W12,model_data1))
  expect_equal(log_posterior_theta(0,W13,model_data1),lpt_inefficient(0,W13,model_data1))
  expect_equal(log_posterior_theta(-5,W13,model_data1),lpt_inefficient(-5,W13,model_data1))
  expect_equal(log_posterior_theta(0,W14,model_data1),lpt_inefficient(0,W14,model_data1))
  expect_equal(log_posterior_theta(-5,W14,model_data1),lpt_inefficient(-5,W14,model_data1))


  expect_equal(log_posterior_theta(0,W21,model_data2),lpt_inefficient(0,W21,model_data2))
  expect_equal(log_posterior_theta(-5,W21,model_data2),lpt_inefficient(-5,W21,model_data2))
  expect_equal(log_posterior_theta(0,W22,model_data2),lpt_inefficient(0,W22,model_data2))
  expect_equal(log_posterior_theta(-5,W22,model_data2),lpt_inefficient(-5,W22,model_data2))
  expect_equal(log_posterior_theta(0,W23,model_data2),lpt_inefficient(0,W23,model_data2))
  expect_equal(log_posterior_theta(-5,W23,model_data2),lpt_inefficient(-5,W23,model_data2))
  expect_equal(log_posterior_theta(0,W24,model_data2),lpt_inefficient(0,W24,model_data2))
  expect_equal(log_posterior_theta(-5,W24,model_data2),lpt_inefficient(-5,W24,model_data2))

  expect_equal(log_posterior_theta(0,W31,model_data3),lpt_inefficient(0,W31,model_data3))
  expect_equal(log_posterior_theta(-5,W31,model_data3),lpt_inefficient(-5,W31,model_data3))
  expect_equal(log_posterior_theta(0,W32,model_data3),lpt_inefficient(0,W32,model_data3))
  expect_equal(log_posterior_theta(-5,W32,model_data3),lpt_inefficient(-5,W32,model_data3))
  expect_equal(log_posterior_theta(0,W33,model_data3),lpt_inefficient(0,W33,model_data3))
  expect_equal(log_posterior_theta(-5,W33,model_data3),lpt_inefficient(-5,W33,model_data3))
  expect_equal(log_posterior_theta(0,W34,model_data3),lpt_inefficient(0,W34,model_data3))
  expect_equal(log_posterior_theta(-5,W34,model_data3),lpt_inefficient(-5,W34,model_data3))

  expect_equal(log_posterior_theta(0,W51,model_data5),lpt_inefficient(0,W51,model_data5))
  expect_equal(log_posterior_theta(-5,W51,model_data5),lpt_inefficient(-5,W51,model_data5))
  expect_equal(log_posterior_theta(0,W52,model_data5),lpt_inefficient(0,W52,model_data5))
  expect_equal(log_posterior_theta(-5,W52,model_data5),lpt_inefficient(-5,W52,model_data5))
  expect_equal(log_posterior_theta(0,W53,model_data5),lpt_inefficient(0,W53,model_data5))
  expect_equal(log_posterior_theta(-5,W53,model_data5),lpt_inefficient(-5,W53,model_data5))
  expect_equal(log_posterior_theta(0,W54,model_data5),lpt_inefficient(0,W54,model_data5))
  expect_equal(log_posterior_theta(-5,W54,model_data5),lpt_inefficient(-5,W54,model_data5))

  expect_equal(log_posterior_theta(c(0,0),W71,model_data7),lpt_inefficient(c(0,0),W71,model_data7))
  expect_equal(log_posterior_theta(c(-5,-5),W71,model_data7),lpt_inefficient(c(-5,-5),W71,model_data7))
  expect_equal(log_posterior_theta(c(0,0),W72,model_data7),lpt_inefficient(c(0,0),W72,model_data7))
  expect_equal(log_posterior_theta(c(-5,-5),W72,model_data7),lpt_inefficient(c(-5,-5),W72,model_data7))
  expect_equal(log_posterior_theta(c(0,0),W73,model_data7),lpt_inefficient(c(0,0),W73,model_data7))
  expect_equal(log_posterior_theta(c(-5,-5),W73,model_data7),lpt_inefficient(c(-5,-5),W73,model_data7))
  expect_equal(log_posterior_theta(c(0,0),W74,model_data7),lpt_inefficient(c(0,0),W74,model_data7))
  expect_equal(log_posterior_theta(c(-5,-5),W74,model_data7),lpt_inefficient(c(-5,-5),W74,model_data7))



})



