# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)
library(microbenchmark)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that('Check correctness of soft-thresholding',{
  expect_equal(soft(10,5), soft_c(10,5)) # First test for fixed values
  
  x <- rnorm(1)
  expect_equal(soft(x,0), soft_c(x,0)) # Second test for random values
})

m1_1 <- microbenchmark(soft(10,5),
                       soft_c(10,5),
                       times = 1000)
m1_1

x <- rnorm(1)
m1_2 <- microbenchmark(soft(x,0),
                       soft_c(x,0),
                       times = 1000)
m1_2



# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)

