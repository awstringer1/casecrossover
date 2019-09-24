context("Optimization")
library(casecrossover)
library(dplyr) # For tests
source("prep-sample-data.R")

# Optimization for a single theta

optcontrol <- cc_control()$opt_control # Defaults
optcontrol$report.freq = 0
optcontrol$report.level = 0


sink("./tmpoutput07687")
opt_single_1_1 <- optimize_latentfield_trustoptim(0,model_data1,optcontrol = optcontrol)
opt_single_1_2 <- optimize_latentfield_trustoptim(-10,model_data1,optcontrol = optcontrol)

opt_single_2_1 <- optimize_latentfield_trustoptim(0,model_data2,optcontrol = optcontrol)
opt_single_2_2 <- optimize_latentfield_trustoptim(-10,model_data2,optcontrol = optcontrol)

opt_single_3_1 <- optimize_latentfield_trustoptim(0,model_data3,optcontrol = optcontrol)
opt_single_3_2 <- optimize_latentfield_trustoptim(-10,model_data3,optcontrol = optcontrol)

opt_single_4_1 <- optimize_latentfield_trustoptim(0,model_data4,optcontrol = optcontrol)
opt_single_4_2 <- optimize_latentfield_trustoptim(-10,model_data4,optcontrol = optcontrol)

opt_single_5_1 <- optimize_latentfield_trustoptim(0,model_data5,optcontrol = optcontrol)
opt_single_5_2 <- optimize_latentfield_trustoptim(-10,model_data5,optcontrol = optcontrol)

opt_single_6_1 <- optimize_latentfield_trustoptim(0,model_data6,optcontrol = optcontrol)
opt_single_6_2 <- optimize_latentfield_trustoptim(-10,model_data6,optcontrol = optcontrol)

opt_single_7_1 <- optimize_latentfield_trustoptim(c(0,0),model_data7,optcontrol = optcontrol)
opt_single_7_2 <- optimize_latentfield_trustoptim(c(-10,-10),model_data7,optcontrol = optcontrol)

opt_single_8_1 <- optimize_latentfield_trustoptim(c(0,0),model_data8,optcontrol = optcontrol)
opt_single_8_2 <- optimize_latentfield_trustoptim(c(-10,-10),model_data8,optcontrol = optcontrol)
file.remove("./tmpoutput07687")
sink()


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

})

# Parallel optimization of theta

thetagrid1 <- list(
  c(0),c(1),c(-1)
)

thetagrid2 <- list(
  c(0,0),c(1,1),c(-1,-1)
)


sink("./tmpoutput07687")
opt_1 <- optimize_all_thetas_parallel(thetagrid1,model_data1,optcontrol = optcontrol)

opt_2 <- optimize_all_thetas_parallel(thetagrid1,model_data2,optcontrol = optcontrol)

opt_3 <- optimize_all_thetas_parallel(thetagrid1,model_data3,optcontrol = optcontrol)

opt_4 <- optimize_all_thetas_parallel(thetagrid1,model_data4,optcontrol = optcontrol)

opt_5 <- optimize_all_thetas_parallel(thetagrid1,model_data5,optcontrol = optcontrol)

opt_6 <- optimize_all_thetas_parallel(thetagrid1,model_data6,optcontrol = optcontrol)

opt_7 <- optimize_all_thetas_parallel(thetagrid2,model_data7,optcontrol = optcontrol)

opt_8 <- optimize_all_thetas_parallel(thetagrid2,model_data8,optcontrol = optcontrol)
file.remove("./tmpoutput07687")
sink()

test_that("Parallel optimization of theta works",{
  # Theta formatting
  expect_error(optimize_all_thetas_parallel(c(1,2),model_data1)) # Not a list
  expect_error(optimize_all_thetas_parallel(list(c("a"),c("b")),model_data1)) # Not a list of numeric vectors
  expect_error(optimize_all_thetas_parallel(list(c(1),c(1,2)),model_data1)) # Different lengths

  expect_s3_class(opt_1,"tbl_df")
  expect_s3_class(opt_2,"tbl_df")
  expect_s3_class(opt_3,"tbl_df")
  expect_s3_class(opt_4,"tbl_df")
  expect_s3_class(opt_5,"tbl_df")
  expect_s3_class(opt_6,"tbl_df")
  expect_s3_class(opt_7,"tbl_df")
  expect_s3_class(opt_8,"tbl_df")

  expect_true(all(purrr::map(opt_1$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_2$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_3$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_4$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_5$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_6$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_7$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
  expect_true(all(purrr::map(opt_8$optimizer,"status") %in% c("Success","Radius of trust region is less than stop.trust.radius")))
})



