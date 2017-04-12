\name{class.BisXP}
\alias{class.BisXP}
\title{Cast bisulfite experiment data into a BisXP object}
\usage{
class.BisXP(BisData)
}
\arguments{
\item{BisData}{A data frame with 4 columns:
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
\code{class.BisXP} creates a BisXP object from a table containing RNA name,
  C position, bisulfite non-conversion ratio and adjusted p-value.
}
\details{
This function takes bisulfite experiment data as an input table, performs checks, and casts the table
  into a BisXP object, with rows labelled after a RNA_C.position pattern.
}
\examples{
 RNA     <- c("NM_00001","NM_00001","NM_00002")
 Cpos     <- as.integer(c(1,5,1))
 ncratio  <- c(0.1,0.5,0.3)
 pv.adj   <- c(0.001,0.1,0.3)
 BSdata   <- data.frame(RNA, Cpos, ncratio, pv.adj, stringsAsFactors = FALSE)
 bsXP     <- class.BisXP(BSdata)
}

