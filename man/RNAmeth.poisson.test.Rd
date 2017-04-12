\name{RNAmeth.poisson.test}
\alias{RNAmeth.poisson.test}
\title{Test RNA Cytosine methylation based on Poisson model distribution}
\usage{
RNAmeth.poisson.test(BisRNA, lambda, method="BH")
}
\arguments{
\item{BisRNA}{    A dataframe containing RNA name, C position, coverage and
non-conversion ratio}

\item{lambda}{    Ratio (Poisson parameter / coverage)}

\item{method}{    Adjustment method for multiple testing, either
                   "BH" (Benjamini-Hochberg) or "IHW" (Independent
                   Hypothesis Weighting, from R package IHW)}
}
\value{
This function returns a BisXP object whose elements correspond to
  RNAs where non-conversion ratio is higher than 
  lambda =  ratio (Poisson parameter / coverage).
  The variables contained in this object are the 
  non-conversion ratio and the adjusted p-value.
}
\description{
\code{RNAmeth.poisson.test} tests RNA bisulfite sequencing non-conversion
  based on a sample and the ratio (Poisson parameter / coverage) for this sample.
}
\details{
This function takes RNA bisulfite sequencing data from one sample
   and the ratio (Poisson parameter / coverage) as inputs. Then, the
   dataset is restricted to those C positions where non-conversion 
   ratio is larger than (Poisson parameter / coverage). 
   Finally, it carries out poisson.test and a correction for multiple 
   testing. The output of function RNAmeth.poisson.test is a BisXP 
   object which contains RNA names, C positions, non-conversion ratios 
   and adjusted p-values. The formatting into a BisXP object guarantees 
   that some consistency checks are satisfied.
}
\examples{

 ## Load data, find out ratio Poison parameter / coverage,
 ## and produce corresponding adjusted p-values.
 data(Bisdata,package="BisRNA")
 lambda1 <- RNAmeth.poisson.par(Bisdata1)$estimate
 BisXP1  <- RNAmeth.poisson.test(Bisdata1,lambda1,method="BH")

 ## Print BisXP1
 BisXP1

 ## Display BisXP1 as a data frame
 BisXP1.df <- data.frame(BisXP1$nonconv.ratio, BisXP1$pv.adj, row.names=BisXP1$RNA.pos)
 BisXP1.df

}

