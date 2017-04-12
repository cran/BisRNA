\name{fisher.method}
\alias{fisher.method}
\title{Function implementing Fisher's method to combine independent p-values}
\usage{
fisher.method(pvalues)
}
\arguments{
\item{pvalues}{A list of p-values}
}
\value{
The p-value combined using Fisher's method.
}
\description{
\code{fisher.method} takes a list of independent p-values and combines them
  using Fisher's method.
}
\details{
This function takes a list of p-values as input, determines the number of
  degrees of freedom (2 * number of p-values), combines the p-values using
  Fisher's method and returns the combined p-value in output.
}
\examples{
list.pv    <- c(0.0001,0.0142,0.0150)
p.combined <- fisher.method(list.pv)

}
\references{
  Fisher RA (1925) Statistical Methods for Research Workers. Edinburg: Oliver and Boyd.

  Fisher RA (1948) Questions and Answers #14. In: Mosteller F, Fisher RA (1948) The American Statistician, 2:30-31 \url{http://www.jstor.org/stable/2681650}
}

