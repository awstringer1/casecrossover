context("Likelihood")
library(casecrossover)
library(dplyr) # For tests
source("prep-sample-data.R")

# Prep data for log-lik
dataprep1 <- prep_data_for_log_lik(W = 1:model_data1$Wd,model_data = model_data1) # Linear, case == 1
dataprep2 <- prep_data_for_log_lik(W = 1:model_data2$Wd,model_data = model_data2) # Linear, case != 1
dataprep3 <- prep_data_for_log_lik(W = 1:model_data3$Wd,model_data = model_data3) # Smooth
dataprep5 <- prep_data_for_log_lik(W = 1:model_data5$Wd,model_data = model_data5) # Smooth and linear
test_that("Data is prepared for likelihood correctly",{
  # Linear, case == 1
  expect_length(dataprep1,2)
  expect_equal(names(dataprep1),c("1","2"))
  expect_equal(dataprep1[[1]],c(1,2))
  expect_equal(dataprep1[[2]],3)
  # Linear, case != 1
  expect_length(dataprep2,2)
  expect_equal(names(dataprep2),c("1","2"))
  expect_equal(dataprep2[[1]],c(1,2))
  expect_equal(dataprep2[[2]],3)
  # Smooth
  expect_length(dataprep3,2)
  expect_equal(names(dataprep3),c("1","2"))
  expect_equal(dataprep3[[1]],c(1,2))
  expect_equal(dataprep3[[2]],3)
  # Smooth and linear
  expect_length(dataprep5,2)
  expect_equal(names(dataprep5),c("1","2"))
  expect_equal(dataprep5[[1]],c(1,2))
  expect_equal(dataprep5[[2]],3)

  # Case day extraction
  expect_equal(get_casedays(model_data1),c(1,1,1))
  expect_equal(get_casedays(model_data2),c(2,2,2))
})


# Compute log-likelihood
# With and without weights
#
# Manually computed values:
W1 <- rep(0,model_data1$Wd)
W2 <- c(rep(0,model_data1$Wd - 1),100000) # The last value, for beta, shouldn't affect the likelihood.
W3 <- 1:model_data1$Wd # Nonzero parameter values
W4 <- rnorm(model_data1$Wd) # Random values- ensure this is robust!

# No weights
ll1_1 <- -1 * ( log(1 + exp(-W1[1]) + exp(-W1[2])) + log(1 + exp(-W1[3])))
ll1_2 <- -1 * ( log(1 + exp(-W2[1]) + exp(-W2[2])) + log(1 + exp(-W2[3])))
ll1_3 <- -1 * ( log(1 + exp(-W3[1]) + exp(-W3[2])) + log(1 + exp(-W3[3])))
ll1_4 <- -1 * ( log(1 + exp(-W4[1]) + exp(-W4[2])) + log(1 + exp(-W4[3])))


# With weights
ll2_1 <- -1 * ( 2 * log(1 + exp(-W1[1]) + exp(-W1[2])) + 2 * log(1 + exp(-W1[3])))
ll2_2 <- -1 * ( 2 * log(1 + exp(-W2[1]) + exp(-W2[2])) + 2 * log(1 + exp(-W2[3])))
ll2_3 <- -1 * ( 2 * log(1 + exp(-W3[1]) + exp(-W3[2])) + 2 * log(1 + exp(-W3[3])))
ll2_4 <- -1 * ( 2 * log(1 + exp(-W4[1]) + exp(-W4[2])) + 2 * log(1 + exp(-W4[3])))


test_that("Log-likelihood is computed correctly",{
  # No weights
  expect_equal(log_likelihood(W1,model_data1),ll1_1)
  expect_equal(log_likelihood(W2,model_data1),ll1_2)
  expect_equal(log_likelihood(W3,model_data1),ll1_3)
  expect_equal(log_likelihood(W4,model_data1),ll1_4)
  # With weights
  expect_equal(log_likelihood(W1,model_data2),ll2_1)
  expect_equal(log_likelihood(W2,model_data2),ll2_2)
  expect_equal(log_likelihood(W3,model_data2),ll2_3)
  expect_equal(log_likelihood(W4,model_data2),ll2_4)
})

# Gradient. Use same W's

d1_1_1 <- 1 + sum(exp(-W1[1:2]))
d1_1_2 <- 1 + sum(exp(-W1[3]))
d1_2_1 <- 1 + sum(exp(-W2[1:2]))
d1_2_2 <- 1 + sum(exp(-W2[3]))
d1_3_1 <- 1 + sum(exp(-W3[1:2]))
d1_3_2 <- 1 + sum(exp(-W3[3]))
d1_4_1 <- 1 + sum(exp(-W4[1:2]))
d1_4_2 <- 1 + sum(exp(-W4[3]))

g1_1 <- c(exp(-W1[1:2])/d1_1_1,exp(-W1[3])/d1_1_2,0)
g1_2 <- c(exp(-W2[1:2])/d1_2_1,exp(-W2[3])/d1_2_2,0)
g1_3 <- c(exp(-W3[1:2])/d1_3_1,exp(-W3[3])/d1_3_2,0)
g1_4 <- c(exp(-W4[1:2])/d1_4_1,exp(-W4[3])/d1_4_2,0)

g2_1 <- 2 * g1_1
g2_2 <- 2 * g1_2
g2_3 <- 2 * g1_3
g2_4 <- 2 * g1_4

test_that("Gradient is computed correctly",{
  # No weights
  expect_equal(grad_log_likelihood(W1,model_data1),g1_1)
  expect_equal(grad_log_likelihood(W2,model_data1),g1_2)
  expect_equal(grad_log_likelihood(W3,model_data1),g1_3)
  expect_equal(grad_log_likelihood(W4,model_data1),g1_4)
  # Weights
  expect_equal(grad_log_likelihood(W1,model_data2),g2_1)
  expect_equal(grad_log_likelihood(W2,model_data2),g2_2)
  expect_equal(grad_log_likelihood(W3,model_data2),g2_3)
  expect_equal(grad_log_likelihood(W4,model_data2),g2_4)
})


# Hessian
# Compute the structure and the values, and combine them

hs <- list(i = as.integer(c(0,1,0,1,2)),p = as.integer(c(0,2,4,5,5)))

hx1_1_1 <- g1_1[1]*(1 - g1_1[1])
hx1_1_2 <- -g1_1[1]*g1_1[2]
hx1_2_2 <- g1_1[2]*(1 - g1_1[2])
hx1_3_3 <- g1_1[3]*(1 - g1_1[3])
hx1 <- c(hx1_1_1,hx1_1_2,hx1_2_2,hx1_3_3)

h1 <- new("dgCMatrix")
h1@i <- hs$i
h1@p <- hs$p
h1@Dim <- as.integer(rep(length(hs$p) - 1,2))
h1@x <- numeric(length(h1@i))
h1 <- as(h1,'symmetricMatrix')
h1@x <- hx1
h1 <- as(h1,'dgCMatrix')

hx2_1_1 <- g1_2[1]*(1 - g1_2[1])
hx2_1_2 <- -g1_2[1]*g1_2[2]
hx2_2_2 <- g1_2[2]*(1 - g1_2[2])
hx2_3_3 <- g1_2[3]*(1 - g1_2[3])
hx2 <- c(hx2_1_1,hx2_1_2,hx2_2_2,hx2_3_3)

h2 <- new("dgCMatrix")
h2@i <- hs$i
h2@p <- hs$p
h2@Dim <- as.integer(rep(length(hs$p) - 1,2))
h2@x <- numeric(length(h2@i))
h2 <- as(h2,'symmetricMatrix')
h2@x <- hx2
h2 <- as(h2,'dgCMatrix')

hx3_1_1 <- g1_3[1]*(1 - g1_3[1])
hx3_1_2 <- -g1_3[1]*g1_3[2]
hx3_2_2 <- g1_3[2]*(1 - g1_3[2])
hx3_3_3 <- g1_3[3]*(1 - g1_3[3])
hx3 <- c(hx3_1_1,hx3_1_2,hx3_2_2,hx3_3_3)

h3 <- new("dgCMatrix")
h3@i <- hs$i
h3@p <- hs$p
h3@Dim <- as.integer(rep(length(hs$p) - 1,2))
h3@x <- numeric(length(h3@i))
h3 <- as(h3,'symmetricMatrix')
h3@x <- hx3
h3 <- as(h3,'dgCMatrix')

hx4_1_1 <- g1_4[1]*(1 - g1_4[1])
hx4_1_2 <- -g1_4[1]*g1_4[2]
hx4_2_2 <- g1_4[2]*(1 - g1_4[2])
hx4_3_3 <- g1_4[3]*(1 - g1_4[3])
hx4 <- c(hx4_1_1,hx4_1_2,hx4_2_2,hx4_3_3)

h4 <- new("dgCMatrix")
h4@i <- hs$i
h4@p <- hs$p
h4@Dim <- as.integer(rep(length(hs$p) - 1,2))
h4@x <- numeric(length(h4@i))
h4 <- as(h4,'symmetricMatrix')
h4@x <- hx4
h4 <- as(h4,'dgCMatrix')

# With the weights
hwx1 <- 2 * hx1
hwx2 <- 2 * hx2
hwx3 <- 2 * hx3
hwx4 <- 2 * hx4
hw1 <- 2 * h1
hw2 <- 2 * h2
hw3 <- 2 * h3
hw4 <- 2 * h4



test_that("Hessian is computed correctly",{
  # Function works as called
  expect_s4_class(hessian_log_likelihood(rep(0,3),model_data1),'CsparseMatrix')
  expect_s4_class(hessian_log_likelihood(rep(0,3),model_data1,hs),'CsparseMatrix')

  # x computed correctly
  expect_equal(hessian_log_likelihood_x(W1,model_data1),hx1)
  expect_equal(hessian_log_likelihood_x(W2,model_data1),hx2)
  expect_equal(hessian_log_likelihood_x(W3,model_data1),hx3)
  expect_equal(hessian_log_likelihood_x(W4,model_data1),hx4)


  # Hessian computed correctly
  expect_equal(hessian_log_likelihood(W1,model_data1),h1)
  expect_equal(hessian_log_likelihood(W2,model_data1),h2)
  expect_equal(hessian_log_likelihood(W3,model_data1),h3)
  expect_equal(hessian_log_likelihood(W4,model_data1),h4)

  ## With weights
  # Function works as called
  expect_s4_class(hessian_log_likelihood(rep(0,3),model_data2),'CsparseMatrix')
  expect_s4_class(hessian_log_likelihood(rep(0,3),model_data2,hs),'CsparseMatrix')

  # x computed correctly
  expect_equal(hessian_log_likelihood_x(W1,model_data2),hwx1)
  expect_equal(hessian_log_likelihood_x(W2,model_data2),hwx2)
  expect_equal(hessian_log_likelihood_x(W3,model_data2),hwx3)
  expect_equal(hessian_log_likelihood_x(W4,model_data2),hwx4)


  # Hessian computed correctly
  expect_equal(hessian_log_likelihood(W1,model_data2),hw1)
  expect_equal(hessian_log_likelihood(W2,model_data2),hw2)
  expect_equal(hessian_log_likelihood(W3,model_data2),hw3)
  expect_equal(hessian_log_likelihood(W4,model_data2),hw4)

})


