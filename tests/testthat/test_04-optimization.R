context("Optimization")
library(casecrossover)
library(dplyr) # For tests
# source("prep-sample-data.R")

# Optimization for a single theta

test_that("Optimization works for a single theta",{
  expect_type(opt_single_1_1,"list")
  expect_type(opt_single_1_2,"list")
  expect_type(opt_single_2_1,"list")
  expect_type(opt_single_2_2,"list")
  expect_type(opt_single_3_1,"list")
  expect_type(opt_single_3_2,"list")
  expect_type(opt_single_4_1,"list")
  expect_type(opt_single_4_2,"list")
  expect_type(opt_single_5_1,"list")
  expect_type(opt_single_5_2,"list")
  expect_type(opt_single_6_1,"list")
  expect_type(opt_single_6_2,"list")
  expect_type(opt_single_7_1,"list")
  expect_type(opt_single_7_2,"list")
  expect_type(opt_single_8_1,"list")
  expect_type(opt_single_8_2,"list")
  expect_type(opt_single_9_1,"list")
  expect_type(opt_single_9_2,"list")
  expect_type(opt_single_10_1,"list")
  expect_type(opt_single_10_2,"list")
  expect_type(opt_single_11_1,"list")
  expect_type(opt_single_11_2,"list")


  expect_true(opt_single_1_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_1_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_2_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_2_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_3_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_3_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_4_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_4_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_5_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_5_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_6_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_6_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_7_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_7_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_8_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_8_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_9_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_9_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_10_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_10_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_11_1$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))
  expect_true(opt_single_11_2$optimizer$status %in% c("Success","Radius of trust region is less than stop.trust.radius"))

})

# Parallel optimization of theta



test_that("Parallel optimization of theta works",{
  # Theta formatting
  expect_error(optimize_all_thetas_parallel(c(1,2),model_data1)) # Not a NIGrid
  expect_error(optimize_all_thetas_parallel(list(c("a"),c("b")),model_data1))
  expect_error(optimize_all_thetas_parallel(list(c(1),c(1,2)),model_data1))

  expect_s3_class(opt_1,"tbl_df")
  expect_s3_class(opt_2,"tbl_df")
  expect_s3_class(opt_3,"tbl_df")
  expect_s3_class(opt_4,"tbl_df")
  expect_s3_class(opt_5,"tbl_df")
  expect_s3_class(opt_6,"tbl_df")
  expect_s3_class(opt_7,"tbl_df")
  expect_s3_class(opt_8,"tbl_df")
  expect_s3_class(opt_9,"tbl_df")
  expect_s3_class(opt_10,"tbl_df")
  expect_s3_class(opt_11,"tbl_df")


  expect_true(all(purrr::map(opt_1$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_2$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_3$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_4$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_5$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_6$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_7$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_8$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_9$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_10$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_11$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))

})



