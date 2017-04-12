\name{read.BisXP}
\alias{read.BisXP}
\title{Read RNA bisulfite experiment data and cast it into a BisXP object}
\usage{
read.BisXP(filename)
}
\arguments{
\item{filename}{Address of the file containing data from a bisulfite experiment in 4 columns separated by a tabulation, with header on the first line:
  \itemize{
  \item RNA name
  \item C position (integer, in [1,+Inf])
  \item bisulfite non-conversion ratio      (numeric, in [0,1])
  \item pvalue.adj (numeric, in [0,1])}}
}
\value{
If the input table is correct, then the output will be a BisXP object
  corresponding to a consolidated data frame with rows labelled after a
  RNA_C.position pattern.
  }
\description{
  \code{read.BisXP} reads a table containing RNA name, C position,
  bisulfite non-conversion ratio and adjusted p-value from a file
  and casts it into a BisXP object.
}
\details{
This function takes a file name as input, reads the bisulfite data table
  contained in this file, performs checks, and casts the data into a BisXP object,
  which contain input data with rows labelled after a RNA_C.position pattern.
}
\examples{
 RNAs     <- c("NM_00001","NM_00001","NM_00002")
 Cpos     <- c(1,5,1)
 ncratio  <- c(0.1,0.5,0.3)
 pv.adj   <- c(0.001,0.1,0.3)
 BSdata <- data.frame(RNAs, Cpos, ncratio, pv.adj)
}

