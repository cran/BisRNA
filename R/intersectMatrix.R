# Take intersection of 2 tables
#
#
# This function takes 2 matrices as input, determines the intersection based on
#  their row names, and outputs a single matrix containing the rows in the
#  intersect and all columns of the initial matrices.
#
# In:
#        Tab1     A matrix or data frame with defined row.names
#        Tab2     A matrix or data frame with defined row.names
#
# Out:
#        A matrix with rows common to both Tab1 and Tab2, and concatenated columns.
#
#
intersectMatrix <- function(Tab1,Tab2) {
  targets       <- intersect(row.names(Tab1), row.names(Tab2))
  TabMatrix  <- cbind(Tab1[targets,],Tab2[targets,])
  return(TabMatrix)}
