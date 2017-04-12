# Function implementing Fisher's method to combine independent p-values
#
#
# This function takes a list of p-values as input, 
#   - determines the number of degrees of freedom (2 * number of p-values),
#   - combines the p-values using Fisher's method.
#
# In:
#        pvalues     A list of p-values (unadjusted)
#
# Out: 
#        A list of p-value combined using Fisher's method.
#
# Reference: 
#   Fisher RA (1925) Statistical Methods for Research Workers. Edinburg: Oliver and Boyd.
#   Fisher RA (1948) Questions and Answers #14. In: Mosteller F, Fisher RA (1948) The American Statistician, 2:30-31
#   url: http://www.jstor.org/stable/2681650
#
#
fisher.method <- function(pvalues)
{
  df <- 2*length(pvalues)
  stats::pchisq( -2*sum(log(pvalues),na.rm=TRUE),
                 df,
                 lower.tail=FALSE )
}
