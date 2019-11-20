context("Post-hoc")
library(casecrossover)
library(dplyr) # For tests
# source("prep-sample-data.R")


# Log-posterior of theta
# The log-posterior function already tested, so just test that the output is consistent
# Also, test the computation of the sigma log-posterior.

test_that("Theta and sigma log-posteriors computed correctly",{
  expect_s3_class(logpost1,"tbl_df")
  expect_equal(logpost1$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost1$sigma[[1]]),logpost1$solution[[1]],model_data1))
  expect_equal(logpost1$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost1$sigma[[2]]),logpost1$solution[[2]],model_data1))
  expect_equal(logpost1$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost1$sigma[[3]]),logpost1$solution[[3]],model_data1))

  expect_s3_class(logpost2,"tbl_df")
  expect_equal(logpost2$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost2$sigma[[1]]),logpost2$solution[[1]],model_data2))
  expect_equal(logpost2$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost2$sigma[[2]]),logpost2$solution[[2]],model_data2))
  expect_equal(logpost2$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost2$sigma[[3]]),logpost2$solution[[3]],model_data2))

  expect_s3_class(logpost3,"tbl_df")
  expect_equal(logpost3$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost3$sigma[[1]]),logpost3$solution[[1]],model_data3))
  expect_equal(logpost3$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost3$sigma[[2]]),logpost3$solution[[2]],model_data3))
  expect_equal(logpost3$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost3$sigma[[3]]),logpost3$solution[[3]],model_data3))

  expect_s3_class(logpost4,"tbl_df")
  expect_equal(logpost4$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost4$sigma[[1]]),logpost4$solution[[1]],model_data4))
  expect_equal(logpost4$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost4$sigma[[2]]),logpost4$solution[[2]],model_data4))
  expect_equal(logpost4$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost4$sigma[[3]]),logpost4$solution[[3]],model_data4))

  expect_s3_class(logpost5,"tbl_df")
  expect_equal(logpost5$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost5$sigma[[1]]),logpost5$solution[[1]],model_data5))
  expect_equal(logpost5$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost5$sigma[[2]]),logpost5$solution[[2]],model_data5))
  expect_equal(logpost5$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost5$sigma[[3]]),logpost5$solution[[3]],model_data5))

  expect_s3_class(logpost6,"tbl_df")
  expect_equal(logpost6$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost6$sigma[[1]]),logpost6$solution[[1]],model_data6))
  expect_equal(logpost6$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost6$sigma[[2]]),logpost6$solution[[2]],model_data6))
  expect_equal(logpost6$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost6$sigma[[3]]),logpost6$solution[[3]],model_data6))

  expect_s3_class(logpost7,"tbl_df")
  expect_equal(logpost7$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost7$sigma[[1]]),logpost7$solution[[1]],model_data7))
  expect_equal(logpost7$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost7$sigma[[2]]),logpost7$solution[[2]],model_data7))
  expect_equal(logpost7$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost7$sigma[[3]]),logpost7$solution[[3]],model_data7))

  expect_s3_class(logpost8,"tbl_df")
  expect_equal(logpost8$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost8$sigma[[1]]),logpost8$solution[[1]],model_data8))
  expect_equal(logpost8$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost8$sigma[[2]]),logpost8$solution[[2]],model_data8))
  expect_equal(logpost8$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost8$sigma[[3]]),logpost8$solution[[3]],model_data8))

  expect_s3_class(logpost9,"tbl_df")
  expect_equal(logpost9$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost9$sigma[[1]]),logpost9$solution[[1]],model_data9))
  expect_equal(logpost9$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost9$sigma[[2]]),logpost9$solution[[2]],model_data9))
  expect_equal(logpost9$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost9$sigma[[3]]),logpost9$solution[[3]],model_data9))

  expect_s3_class(logpost10,"tbl_df")
  expect_equal(logpost10$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost10$sigma[[1]]),logpost10$solution[[1]],model_data10))
  expect_equal(logpost10$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost10$sigma[[2]]),logpost10$solution[[2]],model_data10))
  expect_equal(logpost10$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost10$sigma[[3]]),logpost10$solution[[3]],model_data10))

  expect_s3_class(logpost11,"tbl_df")
  expect_equal(logpost11$sigma_logposterior[1],log_posterior_sigma(as.numeric(logpost11$sigma[[1]]),logpost11$solution[[1]],model_data11))
  expect_equal(logpost11$sigma_logposterior[2],log_posterior_sigma(as.numeric(logpost11$sigma[[2]]),logpost11$solution[[2]],model_data11))
  expect_equal(logpost11$sigma_logposterior[3],log_posterior_sigma(as.numeric(logpost11$sigma[[3]]),logpost11$solution[[3]],model_data11))
})



test_that("Normalizing the posterior works as expected, simple functions",{
  expect_equal(round(normalize_log_posterior(f1vals,x1),2),0)
  expect_equal(normalize_log_posterior(f1vals,x1),log(mvQuad::quadrature(f1exp,x1)))

  expect_equal(round(normalize_log_posterior(f2vals,x2),2),0)
  expect_equal(normalize_log_posterior(f2vals,x2),log(mvQuad::quadrature(f2exp,x2)))

  expect_equal(round(normalize_log_posterior(f3vals,x3),2),0)
  expect_equal(normalize_log_posterior(f3vals,x3),log(mvQuad::quadrature(f3exp,x3)))
})

test_that("Normalizing the posterior works as expected, actual model objects",{
  expect_equal(sum(exp(logpost_norm1$theta_logposterior) * mvQuad::getWeights(attributes(logpost1)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm2$theta_logposterior) * mvQuad::getWeights(attributes(logpost2)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm3$theta_logposterior) * mvQuad::getWeights(attributes(logpost3)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm4$theta_logposterior) * mvQuad::getWeights(attributes(logpost4)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm5$theta_logposterior) * mvQuad::getWeights(attributes(logpost5)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm6$theta_logposterior) * mvQuad::getWeights(attributes(logpost6)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm7$theta_logposterior) * mvQuad::getWeights(attributes(logpost7)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm8$theta_logposterior) * mvQuad::getWeights(attributes(logpost8)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm9$theta_logposterior) * mvQuad::getWeights(attributes(logpost9)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm10$theta_logposterior) * mvQuad::getWeights(attributes(logpost10)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm11$theta_logposterior) * mvQuad::getWeights(attributes(logpost11)$thetagrid)),1)
  expect_equal(sum(exp(logpost_norm13$theta_logposterior) * mvQuad::getWeights(attributes(logpost13)$thetagrid)),1)
})

# Obtaining the correct indices for model terms
test_that("Obtaining indices works as expected",{
  expect_s3_class(index1,"ccindex")
  expect_s3_class(index3,"ccindex")
  expect_s3_class(index5,"ccindex")
  expect_s3_class(index7,"ccindex")
  expect_s3_class(index9,"ccindex")
  expect_s3_class(index11,"ccindex")

  expect_equal(index1$linear,c("x" = 4))
  expect_null(index1$smooth)

  expect_null(index3$linear)
  expect_equal(index3$smooth,c("x" = 4,"x" = 5))

  expect_equal(index5$linear,c("x" = 6))
  expect_equal(index5$smooth,c("x" = 4,"x" = 5))

  expect_null(index7$linear)
  expect_equal(index7$smooth,c("x" = 4,"x" = 5,"x2" = 6,"x2" = 7,"x2" = 8,"x2" = 9))

  expect_equal(index9$linear,c("x" = 10,"x" = 11,"x2" = 12,"x2" = 13,"x2" = 14))
  expect_equal(index9$smooth,c("x" = 4,"x" = 5,"x2" = 6,"x2" = 7,"x2" = 8,"x2" = 9))

  expect_equal(index11$linear,c("x" = 10,"x" = 11,"x2" = 12,"x2" = 13,"x2" = 14))
  expect_equal(index11$smooth,c("x" = 4,"x" = 5,"x2" = 6,"x2" = 7,"x2" = 8,"x2" = 9))

  expect_equal(index13$linear,c("x" = 10,"x" = 11,"x2" = 12,"x2" = 13,"x2" = 14))
  expect_equal(index13$smooth,c("x" = 4,"x" = 5,"x2" = 6,"x2" = 7,"x2" = 8,"x2" = 9))
})

# Linear combinations
test_that("Making model linear combinations works as expected",{
  # Error if not both linear and smooth
  expect_error(make_model_lincombs(model_data1))
  expect_error(make_model_lincombs(model_data3))
  expect_error(make_model_lincombs(model_data7))

  # Both linear and smooth
  expect_s4_class(make_model_lincombs(model_data5),"sparseMatrix")
  expect_s4_class(make_model_lincombs(model_data9),"sparseMatrix")

  expect_equal(make_model_lincombs(model_data5)[ ,1],c(0,0,0,1,0,2))
  expect_equal(make_model_lincombs(model_data5)[ ,2],c(0,0,0,0,1,3))

  expect_equal(make_model_lincombs(model_data9)[ ,1],c(0,0,0,1,0,0,0,0,0,2^(1:2),0,0,0))
  expect_equal(make_model_lincombs(model_data9)[ ,2],c(0,0,0,0,1,0,0,0,0,3^(1:2),0,0,0))
  expect_equal(make_model_lincombs(model_data9)[ ,3],c(0,0,0,0,0,1,0,0,0,0,0,0^(1:3)))
  expect_equal(make_model_lincombs(model_data9)[ ,4],c(0,0,0,0,0,0,1,0,0,0,0,1^(1:3)))
  expect_equal(make_model_lincombs(model_data9)[ ,5],c(0,0,0,0,0,0,0,1,0,0,0,6^(1:3)))
  expect_equal(make_model_lincombs(model_data9)[ ,6],c(0,0,0,0,0,0,0,0,1,0,0,8^(1:3)))

  expect_equal(make_model_lincombs(model_data11)[ ,1],c(0,0,0,1,0,0,0,0,0,2^(1:2),0,0,0))
  expect_equal(make_model_lincombs(model_data11)[ ,2],c(0,0,0,0,1,0,0,0,0,3^(1:2),0,0,0))
  expect_equal(make_model_lincombs(model_data11)[ ,3],c(0,0,0,0,0,1,0,0,0,0,0,0^(1:3)))
  expect_equal(make_model_lincombs(model_data11)[ ,4],c(0,0,0,0,0,0,1,0,0,0,0,1^(1:3)))
  expect_equal(make_model_lincombs(model_data11)[ ,5],c(0,0,0,0,0,0,0,1,0,0,0,6^(1:3)))
  expect_equal(make_model_lincombs(model_data11)[ ,6],c(0,0,0,0,0,0,0,0,1,0,0,8^(1:3)))

  expect_equal(make_model_lincombs(model_data13)[ ,1],c(0,0,0,1,0,0,0,0,0,2^(1:2),0,0,0))
  expect_equal(make_model_lincombs(model_data13)[ ,2],c(0,0,0,0,1,0,0,0,0,3^(1:2),0,0,0))
  expect_equal(make_model_lincombs(model_data13)[ ,3],c(0,0,0,0,0,1,0,0,0,0,0,0^(1:3)))
  expect_equal(make_model_lincombs(model_data13)[ ,4],c(0,0,0,0,0,0,1,0,0,0,0,1^(1:3)))
  expect_equal(make_model_lincombs(model_data13)[ ,5],c(0,0,0,0,0,0,0,1,0,0,0,6^(1:3)))
  expect_equal(make_model_lincombs(model_data13)[ ,6],c(0,0,0,0,0,0,0,0,1,0,0,8^(1:3)))
})

# Linear constraints
test_that("Linear constraints are converted to matrix format correctly",{
  expect_error(make_linear_constraints(model_data1))
  expect_equal(make_linear_constraints(model_data3),0)
  expect_equal(make_linear_constraints(model_data5),0)
  expect_equal(make_linear_constraints(model_data7),0)
  expect_equal(make_linear_constraints(model_data9),0)

  expect_equal(make_linear_constraints(model_data11)@x,1)
  expect_equal(make_linear_constraints(model_data11)@i,7) # Row and column indices are 0-based
  expect_equal(make_linear_constraints(model_data11)@j,0)
  expect_equal(make_linear_constraints(model_data11)@Dim,c(model_data11$Wd,1))

  expect_equal(make_linear_constraints(model_data13)@x,c(1,1))
  expect_equal(make_linear_constraints(model_data13)@i,c(3,7)) # Row and column indices are 0-based
  expect_equal(make_linear_constraints(model_data13)@j,c(0,1))
  expect_equal(make_linear_constraints(model_data13)@Dim,c(model_data11$Wd,2))
})

# Final means and variances
test_that("Post-hoc quantities are computed as expected",{
  # No linear constraints or combinations
  expect_gt(posthoc1$variance,0)
  expect_null(posthoc1$lincombvars)

  expect_gt(posthoc2$variance,0)
  expect_null(posthoc2$lincombvars)

  expect_true(all(posthoc3$variance>0))
  expect_equal(length(posthoc3$mean),length(posthoc3$variance))
  expect_null(posthoc3$lincombvars)

  expect_true(all(posthoc4$variance>0))
  expect_equal(length(posthoc4$mean),length(posthoc4$variance))
  expect_null(posthoc4$lincombvars)

  expect_true(all(posthoc5$variance>0))
  expect_equal(length(posthoc5$mean),length(posthoc5$variance))
  expect_equal(ncol(make_model_lincombs(model_data5)),length(posthoc5$lincombvars))

  expect_true(all(posthoc6$variance>0))
  expect_equal(length(posthoc6$mean),length(posthoc6$variance))
  expect_equal(ncol(make_model_lincombs(model_data6)),length(posthoc6$lincombvars))

  expect_true(all(posthoc7$variance>0))
  expect_equal(length(posthoc7$mean),length(posthoc7$variance))
  expect_null(posthoc7$lincombvars)

  expect_true(all(posthoc8$variance>0))
  expect_equal(length(posthoc8$mean),length(posthoc8$variance))
  expect_null(posthoc8$lincombvars)

  expect_true(all(posthoc9$variance>0))
  expect_equal(length(posthoc9$mean),length(posthoc9$variance))
  expect_equal(ncol(make_model_lincombs(model_data9)),length(posthoc9$lincombvars))

  expect_true(all(posthoc10$variance>0))
  expect_equal(length(posthoc10$mean),length(posthoc10$variance))
  expect_equal(ncol(make_model_lincombs(model_data10)),length(posthoc10$lincombvars))

  # Ones with additional constraints:
  expect_true(all(posthoc11$variance>0))
  expect_equal(length(posthoc11$mean),length(posthoc11$variance))
  expect_equal(ncol(make_model_lincombs(model_data11)),length(posthoc11$lincombvars))
  expect_equal(round(posthoc11$mean[make_linear_constraints(model_data11)@i+1],3),0)

  expect_true(all(posthoc13$variance>0))
  expect_equal(length(posthoc13$mean),length(posthoc13$variance))
  expect_equal(ncol(make_model_lincombs(model_data13)),length(posthoc13$lincombvars))
  expect_equal(round(posthoc13$mean[make_linear_constraints(model_data13)@i+1],3),c(0,0))
})

