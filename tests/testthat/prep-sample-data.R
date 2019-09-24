# Create example data for use in testthat tests.



sampledata <- tibble(
  id = c(1,1,1,2,2),
  case1 = c(0,0,1,0,1),
  case2 = c(0,0,2,0,2),
  x = c(1,2,3,1,2),
  x2 = c(6,2,8,1,0)
)

ff1 <- case1 ~ x + strata(id)
controlsmooth <- cc_control(smooth_prior = list(pc_prior(3,.75)),
                            linear_constraints = create_linear_constraints(u = sampledata$x,
                                                                           whichzero = 1,
                                                                           nm = "x"))

controlsmooth2 <- cc_control(
  smooth_prior = list(
    pc_prior(3,.75),
    pc_prior(5,.1)
  ),
  linear_constraints = c(
    create_linear_constraints(u = sampledata$x,
                              whichzero = 1,
                              nm = "x"),
    create_linear_constraints(u = sampledata$x2,
                              whichzero = 2,
                              nm = "x2")
  )
)

model_data1 <- model_setup(case1 ~ x + strata(id),sampledata)
model_data2 <- model_setup(case2 ~ x + strata(id),sampledata)

model_data3 <- model_setup(case1 ~ s(x) + strata(id),sampledata,controlsmooth)
model_data4 <- model_setup(case2 ~ s(x) + strata(id),sampledata,controlsmooth)

model_data5 <- model_setup(case1 ~ x + s(x) + strata(id),sampledata,controlsmooth)
model_data6 <- model_setup(case2 ~ x + s(x) + strata(id),sampledata,controlsmooth)

model_data7 <- model_setup(case1 ~ s(x) + s(x2) + strata(id),sampledata,controlsmooth2)
model_data8 <- model_setup(case2 ~ s(x) + s(x2) + strata(id),sampledata,controlsmooth2)

model_data9 <- model_setup(case1 ~ s(x) + s(x2) + poly(x,2) + poly(x2,3) + strata(id),sampledata,controlsmooth2)
model_data10 <- model_setup(case2 ~ s(x) + s(x2) + poly(x,2) + poly(x2,3) + strata(id),sampledata,controlsmooth2)



# Optimization

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


# Add log-posteriors

logpost1 <- add_log_posterior_values(opt_1,model_data1)
logpost2 <- add_log_posterior_values(opt_2,model_data2)
logpost3 <- add_log_posterior_values(opt_3,model_data3)
logpost4 <- add_log_posterior_values(opt_4,model_data4)
logpost5 <- add_log_posterior_values(opt_5,model_data5)
logpost6 <- add_log_posterior_values(opt_6,model_data6)
logpost7 <- add_log_posterior_values(opt_7,model_data7)
logpost8 <- add_log_posterior_values(opt_8,model_data8)



# Proper indexing
index1 <- get_indices(model_data1)
index3 <- get_indices(model_data3)
index5 <- get_indices(model_data5)
index7 <- get_indices(model_data7)
index9 <- get_indices(model_data9)


