\docType{package}
\name{BisRNA-package}
\alias{BisRNA}
\alias{BisRNA-package}
\title{Analysis of RNA Cytosine-5 Methylation from bisulfite sequencing experiments}
\description{
BisRNA analyses non-conversion of bisulfite-treated RNA
  of a set of samples, in order to distinguish actual 
  cytosine-5 methylation marks from artifacts. 
  Calculations are done as follows:
\itemize{
  \item each sample's ratio (Poisson parameter / coverage) is determined, assuming a null Poisson distribution (\code{\link{RNAmeth.poisson.par}}),
  \item non-conversion p-values are calculated for each sample (\code{\link{RNAmeth.poisson.test}}),
  \item combined non-conversion p-values and standard error of the non-conversion ratios
  are calculated on the intersection of the set of samples \code{\link{samples.combine}}.
}
A small non-conversion p-value points to methylation of a RNA cytosine, or
another event blocking bisulfite conversion.
}
\examples{

 ## Load data
 data(Bisdata,package="BisRNA")

 ## Ratio (Poisson parameter / coverage), derived p-values, 
 #  adjusted for multiple testing using either BH (here), or IHW if available.
 lambda1 <- RNAmeth.poisson.par(Bisdata1)$estimate
 BisXP1  <- RNAmeth.poisson.test(Bisdata1,lambda1,method="BH")
 #
 lambda2 <- RNAmeth.poisson.par(Bisdata2)$estimate
 BisXP2  <- RNAmeth.poisson.test(Bisdata2,lambda2,method="BH")
 #
 lambda3 <- RNAmeth.poisson.par(Bisdata3)$estimate
 BisXP3  <- RNAmeth.poisson.test(Bisdata3,lambda3,method="BH")

 ## Combined p-values ; median and standard error of
 ## bisulfite non-conversion ratio
 BisXP.combined <- samples.combine(BisXP1,BisXP2,BisXP3)
}

