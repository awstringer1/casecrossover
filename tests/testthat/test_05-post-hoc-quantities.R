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

# Posterior normalization
# Try some density functions, which integrate to 1 (so log(normconst) == 0)
# Test equality to 3 decimal places.
# The beta distribution is a very difficult case; test to 1 decimal place.
x1 <- seq(-10,10,by = .01) # Real line
x2 <- seq(0,50,by = .01) # Positive reals
x3 <- seq(0.00001,1-.00001,by = 0.00001) # (0,1)

test_that("Normalizing the posterior works as expected",{
  expect_equal(round(normalize_log_posterior(dnorm(x1,log=TRUE),x1),3),0)
  expect_equal(round(normalize_log_posterior(dgamma(x2,1,1,log=TRUE),x2),3),0)
  expect_equal(round(normalize_log_posterior(dgamma(x2,2,3,log=TRUE),x2),3),0)
  expect_equal(round(normalize_log_posterior(dbeta(x3,1,1,log=TRUE),x3),3),0)
  expect_equal(round(normalize_log_posterior(dbeta(x3,2,3,log=TRUE),x3),1),0)
  expect_equal(round(normalize_log_posterior(dbeta(x3,.5,1,log=TRUE),x3),1),0)
})

f1 <- function(x,y,z) mvtnorm::dmvnorm(x = c(x,y,z),log = TRUE)
f2 <- function(x,y) dgamma(x,2,3,log=TRUE) + dgamma(y,3,2,log=TRUE)

ttg1 <- expand.grid(seq(-5,5,by=.5),seq(-5,5,by=.5),seq(-5,5,by=.5))
pp1 <- ttg1 %>% rowwise() %>% mutate(fx = f1(Var1,Var2,Var3)) %>% pull(fx)
tt1 <- list()
for (i in 1:nrow(ttg1)) {
  tt1[[i]] <- as.numeric(ttg1[i, ])
}

ttg2 <- expand.grid(seq(0,5,by=.1),seq(0,5,by=.1))
pp2 <- ttg2 %>% rowwise() %>% mutate(fx = f2(Var1,Var2)) %>% pull(fx)
tt2 <- list()
for (i in 1:nrow(ttg2)) {
  tt2[[i]] <- as.numeric(ttg2[i, ])
}

test_that("Normalizing the posterior works as expected, multiple dimensions",{
  expect_equal(round(normalize_log_posterior(pp1,tt1),1),0)
  expect_equal(round(normalize_log_posterior(pp2,tt2),1),0)
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
})

