\name{samples.combine}
\alias{samples.combine}
\title{Combine samples p-values and ratios}
\usage{
samples.combine(BisXP1, ...)
}
\arguments{
\item{BisXP1}{    A BisXP object containing non-conversion ratio and p-value}

\item{...}{       One or more additional samples, in the form of BisXP objects}
}
\value{
This function returns a data frame whose row names correspond to the
  RNA and C position which are present in all samples, and the following
  variables:
  \describe{
    \item{p.adj.combined}{p-value adjusted (done in the preparation of
      the BisXP object) and combined (done here)}
    \item{nonconv.ratio.med}{Median of bisulfite non-conversion ratio for a
      specific RNA and C positions}
    \item{nonconv.ratio.se}{Standard error of bisulfite non-conversion ratio
      for a specific RNA and C positions}
   }
}
\description{
\code{samples.combine} takes RNA bisulfite sequencing samples as input,
  combines their p-values using Fisher's method, and produces median and
  standard error of the bisulfite non-conversion ratios.
}
\details{
This function takes several bisulfite sequencing samples, in form of BisXP
  objects, as inputs. It is recommended to provide at least 3 samples and
  in any case all available, relevant samples. Using RNA and C positions
  present in all samples, the adjusted p-values of each sample are combined
  using Fisher's method. Median and standard error of the non-conversion
  ratio are also given in output.
}
\examples{

 ## Load data
 data(Bisdata,package="BisRNA")

 ## Obtain the ratio (Poisson parameter / coverage), and p-values
 ## adjusted for multiple testing using BH (here) or IHW method.
 lambda1 <- RNAmeth.poisson.par(Bisdata1)$estimate
 BisXP1  <- RNAmeth.poisson.test(Bisdata1,lambda1,method="BH")
 lambda2 <- RNAmeth.poisson.par(Bisdata2)$estimate
 BisXP2  <- RNAmeth.poisson.test(Bisdata2,lambda2,method="BH")
 lambda3 <- RNAmeth.poisson.par(Bisdata3)$estimate
 BisXP3  <- RNAmeth.poisson.test(Bisdata3,lambda3,method="BH")

 ## Combine samples
 BisXP.combined <- samples.combine(BisXP1,BisXP2,BisXP3)

}
\references{
Fisher RA (1925) Statistical Methods for Research Workers. Edinburg: Oliver and Boyd.

Fisher RA (1948) Questions and Answers #14. In: Mosteller F, Fisher RA (1948) The American Statistician, 2:30-31
  \url{http://www.jstor.org/stable/2681650}
}

