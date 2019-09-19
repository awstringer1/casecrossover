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

