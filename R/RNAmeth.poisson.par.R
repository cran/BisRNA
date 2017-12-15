# Find ratio (Poisson parameter / coverage) for one sample of bisulfite-converted RNA
#
#
# This function takes RNA bisulfite sequencing data from one sample as input, 
#   - restricts the dataset to coverage at least 10,
#   - divides data into coverage bins,
#   - applies a fit to a Poisson distribution using fitdistr from package MASS,
#     on each coverage bin,
#   - and returns the ratio (Poisson parameter / coverage) (median and 95% confidence interval).
#
# In:
#        BSrna     A dataframe containing RNA name, C position, coverage and
#                  non-conversion ratio, in this order.
#
# Out:
#       A list containing:
#          - estimate   The ratio (Poisson parameter / coverage) 
#                       (median taken over coverage bins)
#          - bca.ci     A confidence interval for the ratio (Poisson parameter / coverage) 
#                       (bootstrap confidence interval of the median, type "bca")
#
#
RNAmeth.poisson.par <- function(BSrna)
{

  # Rename variables
  names(BSrna) <- c("RNA","Cpos","coverage","ncratio")

  # Keep only if coverage > 10
  BSrna <- BSrna[BSrna$coverage>=10,]

  # Cut to meR<0.3 : assuming that values above 0.3 follow a different model
  txfit <- BSrna$ncratio[BSrna$ncratio<0.3]
  tcov  <- BSrna$coverage[BSrna$ncratio<0.3]
  
  ## Determine coverage bins
  maxtcov <- max(tcov)
  maxhico <- max(1000, 1000*ceiling(maxtcov/1000))
  if (maxtcov>9999) maxhico <- max(10000, 10000*ceiling(maxtcov/10000))
  
  hico <- graphics::hist(tcov,
               breaks=c(seq(0,50,10),seq(60,280,20),seq(300,600,50),seq(1000,maxhico,1000)),
               plot=FALSE)
  if (maxtcov>9999) {
    hico <- graphics::hist(tcov,
               breaks=c(seq(0,50,10),seq(60,280,20),seq(300,600,50),seq(1000,10000,1000), seq(10000,maxhico,10000)),
               plot=FALSE)
  }

  ## Ratio (Poisson parameter / coverage) in each coverage bin
  resuP  <- c()
  for (i in 2:length(hico$mids)-1)
  {
    tcovlocal <- hico$breaks[i+1]
    tofit     <- txfit[tcov>=hico$mids[i] & tcov<hico$mids[(i+1)]] * tcovlocal
    tofit     <- apply(cbind(tofit),1,round)
    # apply only if at least 5 points in bin
    if (length(tofit)>=5)
    {
      fit.P <- MASS::fitdistr(x=tofit,densfun = "poisson")
      resuP <- rbind(resuP,c(hico$breaks[i+1],unlist(fit.P)))
    }
  }


  ## Format table containing ratio (Poisson parameter / coverage)
  colnames(resuP)[1] <-"midCoverage"
  resuP <- data.frame(resuP)


  ## Poisson parameter / coverage
  rateP <- resuP$estimate.lambda/resuP$midCoverage
  
  ## Bootstrap replicates 
  sampleNmed <- function(tseed, data, lg)
  {
    set.seed(tseed)
    return( stats::median( sample(x=data, size=lg, replace = TRUE) ) )
  }

  nech=2000
  seed.vec <- sample(1:100000, size=nech, replace=FALSE)
  b1       <- apply(cbind(seed.vec), 1, sampleNmed, data=rateP, lg=length(rateP))

  ## Bootstrap median and percentile confidence interval
  lambda <- c()
  lambda$estimate  <- stats::median(b1)
  lambda$bootbca95 <- round(stats::quantile( b1, c(0.025,0.975) ),digits = 4)


  ## Return Poisson parameter / coverage (median and bootstrap confidence interval)
  return(lambda)

}
