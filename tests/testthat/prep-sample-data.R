# Create example data for use in testthat tests.



sampledata <- tibble(
  id = c(1,1,1,2,2),
  case1 = c(0,0,1,0,1),
  case2 = c(0,0,2,0,2),
  x = c(1,2,3,1,2)
)

ff1 <- case1 ~ x + strata(id)
controlsmooth <- cc_control(smooth_prior = pc_prior(3,.75),
                            linear_constraints = create_linear_constraints(u = sampledata$x,
                                                                           whichzero = 1,
                                                                           nm = "x"))

model_data1 <- model_setup(case1 ~ x + strata(id),sampledata)
model_data2 <- model_setup(case2 ~ x + strata(id),sampledata)

model_data3 <- model_setup(case1 ~ s(x) + strata(id),sampledata,controlsmooth)
model_data4 <- model_setup(case2 ~ s(x) + strata(id),sampledata,controlsmooth)

model_data5 <- model_setup(case1 ~ x + s(x) + strata(id),sampledata,controlsmooth)
model_data6 <- model_setup(case2 ~ x + s(x) + strata(id),sampledata,controlsmooth)
