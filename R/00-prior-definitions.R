# Script containing prior definitions
#

# PC PRECISION pc.prec
# Returns log of pc precision prior for theta, the log-precision
pcprec <- function(theta,u,alpha) {
  lambda <- -log(alpha)/u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
