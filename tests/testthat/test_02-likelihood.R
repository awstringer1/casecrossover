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
})


# Compute log-likelihood
# With and without weights
#
# Manually computed values:
W1 <- rep(0,model_data1$Wd)
W2 <- c(rep(0,model_data1$Wd),100000) # The last value, for beta, shouldn't affect the likelihood.
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
