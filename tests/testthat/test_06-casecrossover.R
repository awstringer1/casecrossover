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
  expect_s3_class(cc9,"cc_fit")
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
  expect_length(cc9,OUTPUTLENGTH)
  expect_length(cc10,OUTPUTLENGTH)
  expect_length(cc11,OUTPUTLENGTH)
  expect_length(cc13,OUTPUTLENGTH)
})

test_that("S3 generics work as expected",{
  expect_s3_class(summary(cc1),"cc_summary")
  expect_s3_class(summary(cc2),"cc_summary")
  expect_s3_class(summary(cc3),"cc_summary")
  expect_s3_class(summary(cc4),"cc_summary")
  expect_s3_class(summary(cc5),"cc_summary")
  expect_s3_class(summary(cc6),"cc_summary")
  expect_s3_class(summary(cc7),"cc_summary")
  expect_s3_class(summary(cc8),"cc_summary")
  expect_s3_class(summary(cc9),"cc_summary")
  expect_s3_class(summary(cc10),"cc_summary")
  expect_s3_class(summary(cc11),"cc_summary")
  expect_s3_class(summary(cc13),"cc_summary")

  expect_s3_class(plot(cc1),c("cc_plot","list"))
  expect_s3_class(plot(cc2),c("cc_plot","list"))
  expect_s3_class(plot(cc3),c("cc_plot","list"))
  expect_s3_class(plot(cc4),c("cc_plot","list"))
  expect_s3_class(plot(cc5),c("cc_plot","list"))
  expect_s3_class(plot(cc6),c("cc_plot","list"))
  expect_s3_class(plot(cc7),c("cc_plot","list"))
  expect_s3_class(plot(cc8),c("cc_plot","list"))
  expect_s3_class(plot(cc9),c("cc_plot","list"))
  expect_s3_class(plot(cc10),c("cc_plot","list"))
  expect_s3_class(plot(cc11),c("cc_plot","list"))
  expect_s3_class(plot(cc13),c("cc_plot","list"))
})
