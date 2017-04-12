# Here is the documented source code for several functions:
#
# - class.BisXP
# - read.BisXP
#


# Cast bisulfite experiment data into a BisXP object
#
#
# This function takes a table as input, performs checks, and casts the table
#   into a BisXP object containing the input data and with rows labelled
#   after a RNA_C.position name.
#
# In:
#       BisData      A data frame containing 4 variables:
#                        column1. RNA name   (character)
#                        column2. C position (integer, in [[1,+Inf]])
#                        column3. ratio      (numeric, in [0,1])
#                        column4. pvalue.adj (numeric, in [0,1])
#
# Out:  
#       A BisXP object corresponding to a consolidated data frame 
#       with rows labelled after a RNA_C.position name.
#
class.BisXP <- function(BisData) {
  # check matrix or data frame
  if (class(BisData)!="data.frame") stop("class.BisXP requires a data frame as input")
  # check that there are 4 columns
  if (ncol(BisData)!=4) stop("class.BisXP's input should have 4 columns: RNA, C.position, non-conversion ratio, adjusted p-value")
  # check colum 1..4 has class character..float
  if (class(BisData[,1])!="character" && class(BisData[,1])!="factor") stop("Column 1 should contain RNA names, as character or as factor")
  if (class(BisData[,2])!="integer" || min(BisData[,2])<=0) stop("Column 2 should contain cytosine position - starts at 1, integer")
  if (class(BisData[,3])!="numeric" || min(BisData[,3])<0. || max(BisData[,3])>1.) stop("Column 3 should contain non-conversion ratio - numeric >=0 and <=1")
  if (class(BisData[,4])!="numeric" || min(BisData[,4])<0. || max(BisData[,4])>1.) stop("Column 4 should contain adjusted p-values - numeric >=0 and <=1")
  # check number of digits needed to print C position
  if (max(BisData[,2])>=1e8) stop("C position > 100,000,000, please increase this limit, in class.BisXP.R")
  # Row labels
  RNA.pos <- apply(BisData, 1, function(x) paste(x[1],sprintf("%08d", as.integer(x[2])),sep="_"))
  # BisXP object
  Bisx <- data.frame(RNA.pos= RNA.pos,
                     nonconv.ratio= BisData[,3],
                     pv.adj= BisData[,4],
                     stringsAsFactors = FALSE)
  class(Bisx) <- "BisXP"
  return(Bisx)
}


# Read bisulfite experiment data and cast it into a BisXP object
#
#
# This function takes a file name as input, reads the bisulfite data table
#   that this file should contain, performs checks, and casts the data into
#   a BisXP object, which contains the input data and with rows labelled
#   after a RNA_C.position pattern.
#
# In:
#       filename   Address of file containing data in 4 columns separated 
#                  by a tabulation, with header on the first line:
#                  column1. RNA name   (character)
#                  column2. C position (integer, in [[1,+Inf]])
#                  column3. ratio      (numeric, in [0,1])
#                  column4. pvalue.adj (numeric, in [0,1])
#
# Out:  A BisXP object corresponding to a consolidated data frame with
#       rows labelled after a RNA_C.position pattern.
#
#
read.BisXP <- function(filename) {
  BisData  <- utils::read.delim(filename)
  Bisx     <- class.BisXP(BisData)
}
