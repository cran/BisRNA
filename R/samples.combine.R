# Combine samples p-values and ratios
#
#
# This function takes several bisulfite sequencing samples, in form of BisXP
#   objects, as inputs. It is recommended to provide at least 3 samples and
#   in any case all available, relevant samples. 
# Using RNA and C positions present in all samples, the adjusted p-values of
#   each sample are combined using Fisher's method. Median and standard error 
#   of the non-conversion ratio are calculated.
#
# In:
#        BisXP1     A BisXP object containing non-conversion ratio and p-value
#        ...        One or more additional samples, in the form of BisXP objects,
#                   e.g. BisXP2, BisXP3, BisXP4, etc.
#
# Out:
#       A data frame containing:
#         - row names corresponding to the RNA and C position present in all samples,
#         - p.adj.combined      p-value adjusted (done in the preparation of
#                               the BisXP object) and combined (done here),
#         - nonconv.ratio.med   Median of bisulfite non-conversion ratio for a
#                               specific RNA and C position,
#         - nonconv.ratio.se    Standard error of bisulfite non-conversion ratio
#                               for a specific RNA and C position.
#    
# Reference:
#   Fisher RA (1925) Statistical Methods for Research Workers. Edinburg: Oliver and Boyd.
#   Fisher RA (1948) Questions and Answers #14. In: Mosteller F, Fisher RA (1948) The American Statistician, 2:30-31
#   url: http://www.jstor.org/stable/2681650
#
# 
samples.combine <- function(BisXP1,...) {

      ## Read arguments and determine their number
      arg    <- list(BisXP1,...)
      n      <- length(arg)
      if(class(arg[[1]])!="BisXP") stop("Each argument of samples.merge must be a BisXP object. Please use BisXP.class to cast your data into a BisXP object.")


      ## Take intersection
      BisXP.merged <- data.frame(BisXP1$nonconv.ratio, BisXP1$pv.adj, row.names=BisXP1$RNA.pos)
      # Take each additional sample into account
      for (i in 2:n)
      {
        if(class(arg[[i]])!="BisXP") stop("Each argument of samples.merge must be a BisXP object. Please use BisXP.class to cast your data into a BisXP object.")
        BisXPi <- arg[[i]]
        TabXPi <- data.frame(BisXPi$nonconv.ratio, BisXPi$pv.adj, row.names=BisXPi$RNA.pos)
        BisXP.merged.new <- intersectMatrix(BisXP.merged,TabXPi)
        BisXP.merged     <- BisXP.merged.new
      }


      ## Combine p.values and non-conversion ratio
      p.adj.list        <- seq(2,2*n,2)
      ncratio.list      <- seq(1,2*n,2)
      #
      p.adj.combined    <- apply(BisXP.merged[p.adj.list],   1, fisher.method)
      nonconv.ratio.med <- apply(BisXP.merged[ncratio.list], 1, stats::median)
      nonconv.ratio.se  <- apply(BisXP.merged[ncratio.list], 1, stats::sd) / sqrt(n)

      # Cast results in a data frame
      BisXP.combined    <- data.frame(p.adj.combined,
                                      nonconv.ratio.med,
                                      nonconv.ratio.se,
                                      row.names = row.names(BisXP.merged))


      ## Return data frame containing combined p-values,
      ## median and standard error of non-conversion ratio.
      return(BisXP.combined)

}
