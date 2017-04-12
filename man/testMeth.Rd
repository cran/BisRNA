\name{testMeth}
\alias{testMeth}
\title{Apply poisson.test to BS coverage and non-conversion ratio.}
\usage{
testMeth(X, lambda)
}
\arguments{
\item{X}{A list containing coverage as 1st element and non-conversion
ratio as 2nd element, for one C position.}

\item{lambda}{ratio (Poisson parameter / coverage)}
}
\value{
This function returns the p-value from poisson.test.
}
\description{
\code{testMeth} applies poisson.test to a sample from bisulfite-treated RNA experiment taking into account the ratio (Poisson parameter / coverage) characteristic of this sample.
}
\details{
This function takes RNA bisulfite sequencing coverage and non-conversion
  ratio, applies poisson.test and outputs the p-value.
}

