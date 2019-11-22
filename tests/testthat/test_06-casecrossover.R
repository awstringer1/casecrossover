context("Final function")
library(casecrossover)

OUTPUTLENGTH <- 4

test_that("Case crossover works as expected!",{
  # Function works and returns an object of the correct class
  expect_s3_class(cc1,"cc_fit")
  expect_s3_class(cc2,"cc_fit")
  expect_s3_class(cc3,"cc_fit")
  expect_s3_class(cc4,"cc_fit")
  expect_s3_class(cc5,"cc_fit")
  expect_s3_class(cc6,"cc_fit")
  expect_s3_class(cc7,"cc_fit")
  expect_s3_class(cc8,"cc_fit")
  # expect_s3_class(cc9,"cc_fit")
  expect_s3_class(cc10,"cc_fit")
  expect_s3_class(cc11,"cc_fit")
  expect_s3_class(cc13,"cc_fit")

  # Function returns a list of the correct length
  expect_length(cc1,OUTPUTLENGTH)
  expect_length(cc2,OUTPUTLENGTH)
  expect_length(cc3,OUTPUTLENGTH)
  expect_length(cc4,OUTPUTLENGTH)
  expect_length(cc5,OUTPUTLENGTH)
  expect_length(cc6,OUTPUTLENGTH)
  expect_length(cc7,OUTPUTLENGTH)
  expect_length(cc8,OUTPUTLENGTH)
  # expect_length(cc9,OUTPUTLENGTH)
  expect_length(cc10,OUTPUTLENGTH)
  expect_length(cc11,OUTPUTLENGTH)
  expect_length(cc13,OUTPUTLENGTH)
})
