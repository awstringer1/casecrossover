context("Model setup")
library(casecrossover)
library(dplyr) # For tests

# Source utility functions
devtools::load_all()

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

