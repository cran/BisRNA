\name{RNAmeth.poisson.par}
\alias{RNAmeth.poisson.par}
\title{Find ratio (Poisson parameter / coverage) for one sample of bisulfite-converted RNA}
\usage{
RNAmeth.poisson.par(BSrna)
}
\arguments{
\item{BSrna}{A dataframe containing RNA name, C position, coverage and
non-conversion ratio, in this order.}
}
\value{
This function returns a list containing:
  \describe{
    \item{estimate}{The ratio (Poisson parameter / coverage) (median taken over coverage bins)}
    \item{bca.ci}{A confidence interval for the ratio (Poisson parameter / coverage) (bootstrap
      confidence interval of the median, type "bca")}
  }
}
\description{
\code{RNAmeth.poisson.par} determines the ratio (Poisson parameter / coverage),
  based on the assumption that bisulfite-treated RNA's non-conversion ratios
  follow a Poisson distribution.
}
\details{
This function takes RNA bisulfite sequencing data from one sample as input, restricts the dataset to RNAs for which coverage is at least 10,
  and divides data into coverage bins. Then, the Poisson parameter is estimated at each coverage bin using fitdistr from package MASS.
  The ratio (Poisson parameter / coverage) (median and 95\% confidence interval) is returned in output.
}
\examples{

 ## Load data, find out ratio (Poisson parameter / coverage).
 data(Bisdata,package="BisRNA")
 lambda1 <- RNAmeth.poisson.par(Bisdata1)$estimate
 
 ## P-values adjusted for multiple testing, using either BH (here) or IHW method.
 BisXP1  <- RNAmeth.poisson.test(Bisdata1,lambda1,method="BH")

 ## Print BisXP1
 BisXP1

 ## Display BisXP1 as a data frame
 BisXP1.df <- data.frame(BisXP1$nonconv.ratio, BisXP1$pv.adj, row.names=BisXP1$RNA.pos)
 BisXP1.df

}

