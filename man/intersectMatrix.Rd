\name{intersectMatrix}
\alias{intersectMatrix}
\title{Take intersection of 2 tables}
\usage{
intersectMatrix(Tab1, Tab2)
}
\arguments{
\item{Tab1}{A matrix or data frame with defined row.names}

\item{Tab2}{A matrix or data frame with defined row.names}
}
\value{
A matrix with rows common to both Tab1 and Tab2,
   and concatenated columns.
}
\description{
\code{intersectMatrix} takes 2 matrices and outputs their
  intersection based on common row.names.
}
\details{
This function takes 2 matrices as input, determines the intersection of
 their row names, and returns a single matrix containing the rows in the
 intersection and concatenated columns of the initial matrices.
}

