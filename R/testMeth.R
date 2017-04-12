# Apply poisson.test to Bisulfite RNA data's coverage and non-conversion ratio.
#
#
# This function takes RNA bisulfite sequencing coverage and non-conversion
#   ratio, applies poisson.test and returns a p-value.
#
#  In:
#     X           A list containing coverage as 1st element and non-conversion
#                 ratio as 2nd element, for one C position.
#
#     lambda       ratio (Poisson parameter / coverage)
#
# Out:
#     p-value from poisson.test.
#
#
testMeth <- function(X,lambda)
{
  cov     <- X[1]
  ratio   <- X[2]
  pv.pois <- stats::poisson.test(x=as.integer(ratio*cov),T=as.integer(cov),r=lambda)$p.value
  return(pv.pois)
}
