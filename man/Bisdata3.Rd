\docType{data}
\name{Bisdata3}
\alias{Bisdata3}
\title{Synthetic sample 3 of RNA bisulfite sequencing}
\format{Data frame with 330 to 345 rows and the 4 following variables:
\describe{
  \item{RNA}{RNA name or identifier}
  \item{Cpos}{Cytosine position in RNA}
  \item{coverage}{Number of available reads at Cpos}
  \item{ncratio}{Bisulfite non-conversion ratio at Cpos}
}}
\usage{
data(Bisdata)
}
\description{
A dataset containing the RNA names, C position, coverage and 
non-conversion ratio for synthetic sample 3 of a RNA
bisulfite sequencing experiment.
}
\examples{
data(Bisdata)
ls(pattern="Bisdata")
head(Bisdata3)

}

