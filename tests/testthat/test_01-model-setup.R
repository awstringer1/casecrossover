context("Model setup")
library(casecrossover)
library(dplyr) # For tests

test_that("Placeholder test",{
  expect_equal(1,1)
})

# Create some test data
sampledata <- tibble(
  id = c(1,1,1,2,2),
  case1 = c(0,0,1,0,1),
  case2 = c(0,0,2,0,2),
  x = c(1,2,3,1,2)
)

ff1 <- case1 ~ x + strata(id)
control1 <- cc_control(smooth_prior = list(name = "pc.prec",params = c(u = 3,alpha = .75)))

# model_data1 <- model_setup(case1 ~ x + strata(id),sampledata)
# model_data2 <- model_setup(case2 ~ x + strata(id),sampledata)


# test_that("Model data created correctly",{
#   expect_s3_class(model_data1,"ccmodeldata")
#   expect_s3_class(model_data2,"ccmodeldata")
# })

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
  expect_true(validate_prior_distribution(pc_prior(u = 3,alpha = .5))) # Test all the supported ones are valid... lol
  expect_false(validate_prior_distribution(c(1,2,3)))
  expect_false(validate_prior_distribution(list(name = "alexprior")))
  expect_false(validate_prior_distribution(list(name = "pc.prec",params = c(u = 1,a = 2))))
  expect_false(validate_prior_distribution(list(name = "pc.prec",params = c(alpha = 2))))
})

uu <- (1:10)[sample(1:10)]
test_that("Linear constraints created as expected",{
  expect_error(create_linear_constraints(uu,c(1,2,56,8)))
  expect_equal(create_linear_constraints(uu,1)[[1]]$constraint,sparseVector(1,1,10))
  expect_equal(create_linear_constraints(uu,c(1,2,3))[[1]]$constraint,sparseVector(1,1,10))
  expect_equal(create_linear_constraints(uu,c(1,2,3))[[2]]$constraint,sparseVector(1,2,10))
  expect_equal(create_linear_constraints(uu,c(1,2,3))[[3]]$constraint,sparseVector(1,3,10))
  expect_true(validate_linear_constraints(create_linear_constraints(uu,1)))
  expect_true(validate_linear_constraints(create_linear_constraints(uu,c(1,2,3))))
  expect_false(validate_linear_constraints(c(1,2,3)))
  expect_false(validate_linear_constraints(list(list(u = 1,constraint = c(1,2,3)))))
  expect_false(validate_linear_constraints(list(u = 1,constraint = c(1,2,3))))
  expect_false(validate_linear_constraints(list(constraint = c(1,2,3))))
  expect_false(validate_linear_constraints(list(list(constraint = c(1,2,3)))))
  expect_false(validate_linear_constraints(list(list(u = c(1,2,3,4,5),constraint = c(1,2,3)))))
  expect_false(validate_linear_constraints(list(list(u = c(1,2,3),constraint = c(1,2,3)))))
  expect_false(validate_linear_constraints(list(list(u = c(1,2,3,4,5),constraint = sparseVector(1,1,10)))))
})

