# Test RNA methylation based on a Poisson distribution
#
#
# This function takes RNA bisulfite sequencing data from one sample
#   and the ratio (Poisson parameter / coverage) as inputs. 
# Then, the dataset is restricted to those C positions where non-conversion 
#   ratio is larger than (Poisson parameter / coverage). 
# Finally, this function carries out poisson.test and a correction for multiple 
#   testing. 
# The output of function RNAmeth.poisson.test is a BisXP 
#   object containing RNA names, C positions, non-conversion ratios 
#   and adjusted p-values. The formatting into a BisXP object guarantees 
#   that some consistency checks are satisfied.
#
# In:
#        BisRNA     A dataframe containing RNA names, C positions, coverage and
#                   non-conversion ratios
#
#        lambda     Ratio (Poisson parameter / coverage)
# 
#        method     Adjustment method for multiple testing, either
#                   "BH" (Benjamini-Hochberg) or "IHW" (Independent
#                   Hypothesis Weighting, from R package IHW) if available.
# 
# Out:
#        BisXP object whose elements correspond to RNAs whose non-conversion ratios 
#          are higher than the ratio (Poisson parameter / coverage). The variables
#          contained in this object are the RNA names, C positions, non-conversion
#          ratios and adjusted p-values.
#
#
RNAmeth.poisson.test <- function(BisRNA,lambda,method="BH")
{

  ## Filter out RNAs whose non-conversion ratio is lower than (Poisson parameter / coverage)
  non0       <- BisRNA[BisRNA$ncratio>lambda,]

  ## Obtain p-values and adjust for multiple testing with Benjamini-Hochberg or IHW (reweighted Benjamini-Hochberg) method.
  pvalues    <- apply(as.data.frame(non0[,c("coverage","ncratio")]),
                      1,
                      testMeth,
                      lambda)
  
  if (method == "BH") {
    ## Benjamini-Hochberg method
    pvaluesAdj <- stats::p.adjust(pvalues,method = "BH")
    
  } else if (method == "IHW") {
        ## IHW reweighted Benjamini-Hochberg method
        ## Check presence of IHW before adjusting pvalues with IHW.
        if(requireNamespace("IHW", quietly=TRUE)) {
          ihw.obj <- IHW::ihw(pvalues         = pvalues,
                              covariates      = non0[,"ncratio"],
                              alpha           = 0.05,
                              adjustment_type = "BH")
          pvaluesAdj <- IHW::adj_pvalues(ihw.obj)
        } else {
        ## Revert to Benjamini-Hochberg method in case IHW is not available.
          print(" ")
          print("Warning:")
          print("As package IHW is not installed, multiple testing adjustment is done with Benjamini-Hochberg method.")
          print(" ")
          pvaluesAdj <- stats::p.adjust(pvalues,method = "BH")
          
        }
  } else {
    stop('Method chosen in RNAmeth.poisson.test should be either "BH" or "IHW" (if the latter is available).')
  }

  ## Format into BisXP object
  non0WithPv        <- cbind(non0[,c("RNA","Cpos","ncratio")],pvaluesAdj)
  names(non0WithPv) <- c("RNA","Cpos","ncratio", "pv.adj")
  BisXPint          <- class.BisXP(non0WithPv)

  return(BisXPint)

}
