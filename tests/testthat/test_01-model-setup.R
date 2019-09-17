context("Model setup")
library(casecrossover)
library(dplyr) # For tests
source("prep-sample-data.R")


test_that("Placeholder test",{
  expect_equal(1,1)
})

test_that("Control parameters correctly specified",{
  expect_equal(names(cc_default_control()),names(cc_control()))
})

PARSELENGTH <- 5
test_that("Formula parsing works as expected",{
  expect_error(parse_formula(y ~ f(x)))
  expect_error(parse_formula(y ~ x + s(z))) # No strata
  expect_length(parse_formula(y ~ x + z + s(zz) + strata(id)),PARSELENGTH)
  expect_length(parse_formula(y ~ x + strata(id)),PARSELENGTH)
  expect_length(parse_formula(y ~ s(z) + strata(id)),PARSELENGTH)
  expect_length(parse_formula(y ~ x + s(z) + strata(id)),PARSELENGTH)
})

test_that("Priors created as expected",{
  expect_true(validate_prior_distribution(pc_prior(u = 3,alpha = .5),verbose = FALSE)) # Test all the supported ones are valid...
  expect_false(validate_prior_distribution(c(1,2,3),verbose = FALSE))
  expect_false(validate_prior_distribution(list(name = "alexprior"),verbose = FALSE))
  expect_false(validate_prior_distribution(list(name = "pc.prec",params = c(u = 1,a = 2)),verbose = FALSE))
  expect_false(validate_prior_distribution(list(name = "pc.prec",params = c(alpha = 2)),verbose = FALSE))
})

uu <- (1:10)[sample(1:10)]
test_that("Linear constraints created as expected",{
  expect_error(create_linear_constraints(uu,c(1,2,56,8)))
  expect_equal(create_linear_constraints(uu,1)[[1]]$constraint[[1]],sparseVector(1,1,10))
  expect_equal(create_linear_constraints(uu,c(1,2,3))[[1]]$constraint[[1]],sparseVector(1,1,10))
  expect_equal(create_linear_constraints(uu,c(1,2,3))[[1]]$constraint[[2]],sparseVector(1,2,10))
  expect_equal(create_linear_constraints(uu,c(1,2,3))[[1]]$constraint[[3]],sparseVector(1,3,10))
  expect_true(validate_linear_constraints(create_linear_constraints(uu,1),verbose = FALSE))
  expect_true(validate_linear_constraints(create_linear_constraints(uu,c(1,2,3)),verbose = FALSE))
  expect_false(validate_linear_constraints(c(1,2,3),verbose = FALSE))
  expect_false(validate_linear_constraints(list(list(u = 1,constraint = list(c(1,2,3)))),verbose = FALSE))
  expect_false(validate_linear_constraints(list(u = 1,constraint = list(c(1,2,3))),verbose = FALSE))
  expect_false(validate_linear_constraints(list(constraint = list(c(1,2,3))),verbose = FALSE))
  expect_false(validate_linear_constraints(list(list(constraint = list(c(1,2,3)))),verbose = FALSE))
  expect_false(validate_linear_constraints(list(list(u = c(1,2,3,4,5),constraint = list(c(1,2,3)))),verbose = FALSE))
  expect_false(validate_linear_constraints(list(list(u = c(1,2,3),constraint = list(c(1,2,3)))),verbose = FALSE))
  expect_false(validate_linear_constraints(list(list(u = c(1,2,3,4,5),constraint = list(sparseVector(1,1,10)))),verbose = FALSE))

  # Check it works with duplicates
  expect_true(validate_linear_constraints(create_linear_constraints(c(1,1,1,2),1),verbose = FALSE))
  expect_equal(create_linear_constraints(c(1,1,1,2),1)[[1]]$constraint[[1]],sparseVector(1,1,2))
})

### Big part: testing the actual model data is created correctly ###


test_that("Model data created correctly",{
  # Make sure the function actually returns something
  expect_s3_class(model_data1,"ccmodeldata")
  expect_s3_class(model_data2,"ccmodeldata")
  # Make sure the formula parsing rejects incorrect formulae
  expect_error(model_setup(case1 ~ x,sampledata)) # No strata
  expect_error(model_setup(case1 ~ f(x))) # Wrong elements
  # Check creation of linear terms works correctly
  expect_equal(model_data1$p,1)
  expect_equal(ncol(model_data1$X),1)
  expect_equal(nrow(model_data1$X),5)
  expect_equal(ncol(model_data1$Xd),1)
  expect_equal(nrow(model_data1$Xd),3)
  expect_s4_class(model_data1$X,"CsparseMatrix")
  expect_s4_class(model_data1$Xd,"CsparseMatrix")
  expect_null(model_data1$A)
  expect_length(model_data1$model_elements$smooth,0)

  # Check creation of smooth terms works correctly
  expect_warning(model_setup(case1 ~ s(x) + strata(id),sampledata)) # No constraints specified
  expect_equal(model_data3$M,2)
  expect_equal(names(model_data3$A),"x")
  expect_length(model_data3$A,1)
  expect_s4_class(model_data3$A$x$A,"CsparseMatrix")
  expect_s4_class(model_data3$A$x$Ad,"CsparseMatrix")
  expect_equal(ncol(model_data3$A$x$A),3)
  expect_equal(nrow(model_data3$A$x$A),5)
  expect_equal(ncol(model_data3$A$x$Ad),3)
  expect_equal(nrow(model_data3$A$x$Ad),3)
  expect_length(model_data3$model_elements$linear,0)
  expect_equal(model_data3$vectorofcolumnstoremove,1)

  # Check that creating smooth and linear terms with the same variable works
  expect_equal(model_data1$p,1)
  expect_equal(ncol(model_data5$X),1)
  expect_equal(nrow(model_data5$X),5)
  expect_equal(ncol(model_data5$Xd),1)
  expect_equal(nrow(model_data5$Xd),3)
  expect_s4_class(model_data5$X,"CsparseMatrix")
  expect_s4_class(model_data5$Xd,"CsparseMatrix")

  expect_equal(model_data5$M,2)
  expect_equal(names(model_data5$A),"x")
  expect_length(model_data5$A,1)
  expect_s4_class(model_data5$A$x$A,"CsparseMatrix")
  expect_s4_class(model_data5$A$x$Ad,"CsparseMatrix")
  expect_equal(ncol(model_data5$A$x$A),3)
  expect_equal(nrow(model_data5$A$x$A),5)
  expect_equal(ncol(model_data5$A$x$Ad),3)
  expect_equal(nrow(model_data5$A$x$Ad),3)
  expect_equal(model_data5$vectorofcolumnstoremove,1)

  # Check control days and case days correctly recorded
  expect_equal(model_data1$control_days,c("1" = 2,"2" = 1))
  expect_equal(model_data2$control_days,c("1" = 2,"2" = 1))
  expect_equal(model_data3$control_days,c("1" = 2,"2" = 1))
  expect_equal(model_data4$control_days,c("1" = 2,"2" = 1))
  expect_equal(model_data5$control_days,c("1" = 2,"2" = 1))
  expect_equal(model_data6$control_days,c("1" = 2,"2" = 1))

  expect_equal(model_data1$case_days,c("1" = 1,"2" = 1))
  expect_equal(model_data2$case_days,c("1" = 2,"2" = 2))
  expect_equal(model_data3$case_days,c("1" = 1,"2" = 1))
  expect_equal(model_data4$case_days,c("1" = 2,"2" = 2))
  expect_equal(model_data5$case_days,c("1" = 1,"2" = 1))
  expect_equal(model_data6$case_days,c("1" = 2,"2" = 2))

  # Check dimensions correctly calculated
  expect_equal(model_data1$Ne,5)
  expect_equal(model_data1$Nd,3)
  expect_equal(model_data1$n,2)

  expect_equal(model_data2$Ne,5)
  expect_equal(model_data2$Nd,3)
  expect_equal(model_data2$n,2)

  expect_equal(model_data3$Ne,5)
  expect_equal(model_data3$Nd,3)
  expect_equal(model_data3$n,2)

  expect_equal(model_data4$Ne,5)
  expect_equal(model_data4$Nd,3)
  expect_equal(model_data4$n,2)

  expect_equal(model_data5$Ne,5)
  expect_equal(model_data5$Nd,3)
  expect_equal(model_data5$n,2)

  expect_equal(model_data6$Ne,5)
  expect_equal(model_data6$Nd,3)
  expect_equal(model_data6$n,2)

  ## Test that the prior is specified with the correct parameters

  # Function should always take one argument
  expect_error(model_data1$theta_logprior())
  expect_error(model_data1$theta_logprior(1,1))
  expect_error(model_data3$theta_logprior())
  expect_error(model_data3$theta_logprior(1,1))
  expect_error(model_data5$theta_logprior())
  expect_error(model_data5$theta_logprior(1,1))

  # When no smooth terms specified, function just returns 0
  expect_equal(model_data1$theta_logprior(0),0)

  # When smooth terms specified, function should return the log-prior with the correct parameters
  expect_equal(model_data3$theta_logprior(1),pcprec(1,3,.75))
  expect_equal(model_data3$theta_logprior(0),pcprec(0,3,.75))
  expect_equal(model_data3$theta_logprior(-1),pcprec(-1,3,.75))
  expect_equal(
    model_setup(case1 ~ s(x) + strata(id),sampledata,
                control = cc_control(smooth_prior = pc_prior(1,.5),
                                     linear_constraints = create_linear_constraints(u = sampledata$x,
                                                                                    whichzero = 1,
                                                                                    nm = "x"))
    )$theta_logprior(1),pcprec(1,1,.5)
  ) # Different parameters

  # Linear constraints
  expect_warning(model_setup(case1 ~ s(x) + strata(id),sampledata)) # Smooth terms but no constraints.

})




