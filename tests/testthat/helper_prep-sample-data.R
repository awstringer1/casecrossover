# Create example data for use in testthat tests.



sampledata <- dplyr::tibble(
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

controlsmooth3 <- cc_control(
  smooth_prior = list(
    pc_prior(3,.75),
    pc_prior(5,.1)
  ),
  linear_constraints = c(
    create_linear_constraints(u = sampledata$x,
                              whichzero = 1,
                              nm = "x"),
    create_linear_constraints(u = sampledata$x2,
                              whichzero = c(2,6),
                              nm = "x2")
  )
)

controlsmooth4 <- cc_control(
  smooth_prior = list(
    pc_prior(3,.75),
    pc_prior(5,.1)
  ),
  linear_constraints = c(
    create_linear_constraints(u = sampledata$x,
                              whichzero = c(1,2),
                              nm = "x"),
    create_linear_constraints(u = sampledata$x2,
                              whichzero = c(2,6),
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

model_data9 <- model_setup(case1 ~ s(x) + s(x2) + poly(x,2,raw = TRUE) + poly(x2,3,raw=TRUE) + strata(id),sampledata,controlsmooth2)
model_data10 <- model_setup(case2 ~ s(x) + s(x2) + poly(x,2) + poly(x2,3) + strata(id),sampledata,controlsmooth2)

model_data11 <- model_setup(case1 ~ s(x) + s(x2) + poly(x,2) + poly(x2,3) + strata(id),sampledata,controlsmooth3)

model_data13 <- model_setup(case1 ~ s(x) + s(x2) + poly(x,2) + poly(x2,3) + strata(id),sampledata,controlsmooth4)




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

opt_single_9_1 <- optimize_latentfield_trustoptim(c(0,0),model_data9,optcontrol = optcontrol)
opt_single_9_2 <- optimize_latentfield_trustoptim(c(-10,-10),model_data9,optcontrol = optcontrol)

opt_single_10_1 <- optimize_latentfield_trustoptim(c(0,0),model_data10,optcontrol = optcontrol)
opt_single_10_2 <- optimize_latentfield_trustoptim(c(-10,-10),model_data10,optcontrol = optcontrol)

opt_single_11_1 <- optimize_latentfield_trustoptim(c(0,0),model_data11,optcontrol = optcontrol)
opt_single_11_2 <- optimize_latentfield_trustoptim(c(-10,-10),model_data11,optcontrol = optcontrol)

opt_single_13_1 <- optimize_latentfield_trustoptim(c(0,0),model_data13,optcontrol = optcontrol)
opt_single_13_2 <- optimize_latentfield_trustoptim(c(-10,-10),model_data13,optcontrol = optcontrol)



file.remove("./tmpoutput07687")
sink()


# thetagrid1 <- list(
#   c(0),c(1),c(-1)
# )
#
# thetagrid2 <- list(
#   c(0,0),c(1,1),c(-1,-1)
# )

thetagrid1 <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 3)
mvQuad::rescale(thetagrid1,domain = matrix(c(-1,1),ncol=1))

thetagrid2 <- mvQuad::createNIGrid(dim = 2,type = "GLe",level = 3)
mvQuad::rescale(thetagrid2,domain = matrix(c(-1,-1,1,1),ncol=2))




sink("./tmpoutput07687")
opt_1 <- optimize_all_thetas_parallel(thetagrid1,model_data1,optcontrol = optcontrol)

opt_2 <- optimize_all_thetas_parallel(thetagrid1,model_data2,optcontrol = optcontrol)

opt_3 <- optimize_all_thetas_parallel(thetagrid1,model_data3,optcontrol = optcontrol)

opt_4 <- optimize_all_thetas_parallel(thetagrid1,model_data4,optcontrol = optcontrol)

opt_5 <- optimize_all_thetas_parallel(thetagrid1,model_data5,optcontrol = optcontrol)

opt_6 <- optimize_all_thetas_parallel(thetagrid1,model_data6,optcontrol = optcontrol)

opt_7 <- optimize_all_thetas_parallel(thetagrid2,model_data7,optcontrol = optcontrol)

opt_8 <- optimize_all_thetas_parallel(thetagrid2,model_data8,optcontrol = optcontrol)

opt_9 <- optimize_all_thetas_parallel(thetagrid2,model_data9,optcontrol = optcontrol)

opt_10 <- optimize_all_thetas_parallel(thetagrid2,model_data10,optcontrol = optcontrol)

opt_11 <- optimize_all_thetas_parallel(thetagrid2,model_data11,optcontrol = optcontrol)

opt_13 <- optimize_all_thetas_parallel(thetagrid2,model_data13,optcontrol = optcontrol)

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
logpost9 <- add_log_posterior_values(opt_9,model_data9)
logpost10 <- add_log_posterior_values(opt_10,model_data10)
logpost11 <- add_log_posterior_values(opt_11,model_data11)
logpost13 <- add_log_posterior_values(opt_13,model_data13)



# Proper indexing
index1 <- get_indices(model_data1)
index3 <- get_indices(model_data3)
index5 <- get_indices(model_data5)
index7 <- get_indices(model_data7)
index9 <- get_indices(model_data9)
index11 <- get_indices(model_data11)
index13 <- get_indices(model_data13)


# Posterior normalization
# Try some density functions, which integrate to 1 (so log(normconst) == 0)
# Test equality to 3 decimal places.
# The beta distribution is a very difficult case; test to 1 decimal place.
# x1 <- seq(-10,10,by = .01) # Real line
# x2 <- seq(0,50,by = .01) # Positive reals
# x3 <- seq(0.00001,1-.00001,by = 0.00001) # (0,1)

x1 <- mvQuad::createNIGrid(dim = 3,type = "GLe",level = 20,ndConstruction = "product")
mvQuad::rescale(x1, domain = matrix(c(-10,-10,-10, 10,10,10), ncol=3))
x2 <- mvQuad::createNIGrid(dim = 2,type = "GLe",level = 20,ndConstruction = "product")
mvQuad::rescale(x2, domain = matrix(c(0,0,20,20), ncol=2))
x3 <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 20,ndConstruction = "product")
mvQuad::rescale(x3, domain = matrix(c(-10, 10), ncol=1))


f1 <- function(x) mvtnorm::dmvnorm(x,log = TRUE)
f2 <- function(x) dgamma(x[,1],2,3,log=TRUE) + dgamma(x[,2],3,2,log=TRUE)
f3 <- function(x) dnorm(x,0,1,log = TRUE)

f1exp <- function(x) exp(f1(x))
f2exp <- function(x) exp(f2(x))
f3exp <- function(x) exp(f3(x))

f1vals <- f1(mvQuad::getNodes(x1))
f2vals <- f2(mvQuad::getNodes(x2))
f3vals <- f3(mvQuad::getNodes(x3))

# Try it for the actual model objects
logpost_norm1 <- normalize_optresults_logpost(logpost1)
logpost_norm2 <- normalize_optresults_logpost(logpost2)
logpost_norm3 <- normalize_optresults_logpost(logpost3)
logpost_norm4 <- normalize_optresults_logpost(logpost4)
logpost_norm5 <- normalize_optresults_logpost(logpost5)
logpost_norm6 <- normalize_optresults_logpost(logpost6)
logpost_norm7 <- normalize_optresults_logpost(logpost7)
logpost_norm8 <- normalize_optresults_logpost(logpost8)
logpost_norm9 <- normalize_optresults_logpost(logpost9)
logpost_norm10 <- normalize_optresults_logpost(logpost10)
logpost_norm11 <- normalize_optresults_logpost(logpost11)
logpost_norm13 <- normalize_optresults_logpost(logpost13)

# Final means and variances
posthoc1 <- compute_marginal_means_and_variances(logpost1,model_data1)
posthoc2 <- compute_marginal_means_and_variances(logpost2,model_data2)
posthoc3 <- compute_marginal_means_and_variances(logpost3,model_data3)
posthoc4 <- compute_marginal_means_and_variances(logpost4,model_data4)
posthoc5 <- compute_marginal_means_and_variances(logpost5,model_data5)
posthoc6 <- compute_marginal_means_and_variances(logpost6,model_data6)
posthoc7 <- compute_marginal_means_and_variances(logpost7,model_data7)
posthoc8 <- compute_marginal_means_and_variances(logpost8,model_data8)
posthoc9 <- compute_marginal_means_and_variances(logpost9,model_data9)
posthoc10 <- compute_marginal_means_and_variances(logpost10,model_data10)
posthoc11 <- compute_marginal_means_and_variances(logpost11,model_data11)
posthoc13 <- compute_marginal_means_and_variances(logpost13,model_data13)

