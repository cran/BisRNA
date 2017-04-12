#
# Tests of all routines, in this order:
# ------------------------------------
#
# - class.BisXP
# - read.BisXP
# - fisher.method
# - intersectMatrix
# - RNAmeth.poisson.par
# - RNAmeth.poisson.test
# - samples.combine
#
#

### Test class.BisXP
RNA     <- c("NM_00001","NM_00001","NM_00002")
Cpos     <- as.integer(c(1,5,1))
Cpos2    <- as.integer(c(1,100000005,1))
ncratio  <- c(0.1,0.5,0.3)
pv.adj   <- c(0.001,0.1,0.3)
BSdata   <- data.frame(RNA, Cpos, ncratio, pv.adj, stringsAsFactors = FALSE)
BSdata2  <- as.matrix(BSdata)
BSdata3  <- data.frame(RNA, Cpos, ncratio, stringsAsFactors = FALSE)
BSdata4  <- data.frame(Cpos, RNA, ncratio, pv.adj, stringsAsFactors = FALSE)
BSdata5  <- data.frame(RNA, ncratio, Cpos, pv.adj, stringsAsFactors = FALSE)
BSdata6  <- data.frame(RNA, Cpos, Cpos, pv.adj, stringsAsFactors = FALSE)
BSdata7  <- data.frame(RNA, Cpos, ncratio, Cpos, stringsAsFactors = FALSE)
BSdata8  <- data.frame(RNA, Cpos2, ncratio, pv.adj, stringsAsFactors = FALSE)
bsXP     <- class.BisXP(BSdata)

testthat::test_that("Does class.BisXP yield BisXP object or throw an error?", {
  testthat::expect_equal(class(bsXP), "BisXP")
  testthat::expect_error(class.BisXP(BSdata2),"class.BisXP requires a data frame as input")
  testthat::expect_error(class.BisXP(BSdata3),"class.BisXP's input should have 4 columns: RNA, C.position, non-conversion ratio, adjusted p-value")
  testthat::expect_error(class.BisXP(BSdata4),"Column 1 should contain RNA names, as character or as factor")
  testthat::expect_error(class.BisXP(BSdata5),"Column 2 should contain cytosine position - starts at 1, integer")
  testthat::expect_error(class.BisXP(BSdata6),"Column 3 should contain non-conversion ratio - numeric >=0 and <=1")
  testthat::expect_error(class.BisXP(BSdata7),"Column 4 should contain adjusted p-values - numeric >=0 and <=1")
  testthat::expect_error(class.BisXP(BSdata8),"C position > 100,000,000, please increase this limit, in class.BisXP.R")
  }
  )

### Test read.BisXP
#write.table(BSdata,"inst/BSdata.tabular",quote = FALSE, sep = '\t' ,row.names = FALSE)
fichier <- system.file("testdata","BSdata.tabular", package = "BisRNA")
BSdata9 <- read.BisXP(fichier)

testthat::test_that("Does read.BisXP correctly read BisXP object?", {
  testthat::expect_equal(class(BSdata9), "BisXP")
  }
  )

### Test Fisher's method
# Check result is identical to third-party calculation:
# URL: https://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/
# (Date: 12.09.2016 11:25)
#x = c(1e-3,1e-3,1e-3,1)
#fishersMethod(x)
#[1] 1.719731e-06
x  <- c(1e-3,1e-3,1e-3,1)
px <- fisher.method(x)

testthat::test_that("Does fisher.method works correctly?", {
  testthat::expect_equal(px,1.719731e-06)
  }
  )

### Test intersectMatrix 
BSinter <- intersectMatrix(as.data.frame(BSdata,
                              row.names=c("NM_00001_1","NM_00001_5","NM_00002_1")),
                           as.data.frame(BSdata8,
                              row.names=c("NM_00001_1","NM_00001_10005","NM_00002_1")) )
                              
testthat::test_that("Does intersectMatrix work?", {
  testthat::expect_equal(nrow(BSinter),2)
  }
  )

### Test RNAmeth.poisson.par
RNA      <- c("NM1")
Cpos     <- 1:1000
coverage <- c(4:1003)
set.seed(1)
ncratio  <- rpois(n = 1000,lambda = c(4:1003)*0.1) / 4:1003
BSrna    <- data.frame(RNA, Cpos, coverage, ncratio, stringsAsFactors = FALSE)
lambdaetc <- RNAmeth.poisson.par(BSrna)

testthat::test_that("Does RNAmeth.poisson.par produce the actual rate?", {
  testthat::expect_lt(lambdaetc$bootbca95[1],0.1)
  testthat::expect_gt(lambdaetc$bootbca95[2],0.1)
  }
  )


### Test RNAmeth.poisson.test
# RNAmeth.poisson.test is simply a wrapper of poisson.test with some cutoffs 
# and specific names, therefore the calculation itself is not checked, 
# but only the format.

data(Bisdata,package="BisRNA")
lambda1 <- RNAmeth.poisson.par(Bisdata1)$estimate
BisXP1  <- RNAmeth.poisson.test(Bisdata1,lambda1,method="BH")

testthat::test_that("Does RNAmeth.poisson.test produce the correct output?", {
  testthat::expect_equal(class(BisXP1),"BisXP")
  }
  )

### Test samples.combine
# samples.combine makes use of fisher.method and intersect.matrix,
# which where tested previously. So, here also, only the format is tested

data(Bisdata,package="BisRNA")
lambda1 <- RNAmeth.poisson.par(Bisdata1)$estimate
BisXP1  <- RNAmeth.poisson.test(Bisdata1,lambda1,method="BH")
lambda2 <- RNAmeth.poisson.par(Bisdata2)$estimate
BisXP2  <- RNAmeth.poisson.test(Bisdata2,lambda2,method="BH")
lambda3 <- RNAmeth.poisson.par(Bisdata3)$estimate
BisXP3  <- RNAmeth.poisson.test(Bisdata3,lambda3,method="BH")
BisXP.combined <- samples.combine(BisXP1,BisXP2,BisXP3)

testthat::test_that("Does samples.combine produce the expected output?", {
  testthat::expect_equal(class(BisXP.combined),"data.frame")
  testthat::expect_equal(ncol(BisXP.combined),3)
  testthat::expect_lte(nrow(BisXP.combined),length(BisXP1$RNA.pos))
  }
  )



