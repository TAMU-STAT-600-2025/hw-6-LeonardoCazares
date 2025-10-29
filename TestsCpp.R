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

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
Xtilde <- matrix(c(1, 2,
                   3, 4), nrow = 2, ncol = 2, byrow = TRUE)
Ytilde <- c(1, -1)
beta   <- c(0.5, -0.25)
lambda <- 0.1

Xtilde_big <- matrix(1:2500,
                     nrow = 50, ncol = 50, byrow = TRUE)
Ytilde_big <- -24:25
beta_big   <- rnorm(50)
lambda_big <- 0.1

test_that('Check correctness of lasso objective function', {
  out_r   <- lasso(Xtilde, Ytilde, beta, lambda)
  out_cpp <- lasso_c(Xtilde, Ytilde, beta, lambda)
  
  expect_equal(out_cpp, out_r, tolerance = 1e-12)
  
  out_r   <- lasso(Xtilde_big, Ytilde_big, beta_big, lambda_big)
  out_cpp <- lasso_c(Xtilde_big, Ytilde_big, beta_big, lambda_big)
  
  expect_equal(out_cpp, out_r, tolerance = 1e-12)
})

m2_1 <- microbenchmark(lasso(Xtilde, Ytilde, beta, lambda),
                       lasso_c(Xtilde, Ytilde, beta, lambda),
                       times = 1000)
m2_1

m2_2 <- microbenchmark(lasso(Xtilde_big, Ytilde_big, beta_big, lambda_big),
                       lasso_c(Xtilde_big, Ytilde_big, beta_big, lambda_big),
                       times = 100)
m2_2

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Small datasets
X_small <- matrix(c(1, 2,
              2, 3,
              3, 4),
            nrow = 3, byrow = TRUE)

# Center and scale data
Xtilde_small <- scale(X_small, center = TRUE, scale = TRUE)
Y <- c(1, -1, 2)
Ytilde_small <- as.numeric(scale(Y, center = TRUE, scale = FALSE))

lambda_small <- 0.5

test_that("fitLASSOstandardized_c matches R implementation on a small examples", {
  
  # R version
  fit_r <- fitLASSOstandardized(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda = lambda_small,
    beta_start = NULL,
    eps = 1e-4
  )
  beta_r  <- fit_r$beta

  # Rccp Arma. version
  beta_cpp <- fitLASSOstandardized_c(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda = lambda_small,
    beta_start = numeric(0),
    eps = 1e-4
  )
  
  expect_equal(as.numeric(beta_cpp),
               as.numeric(beta_r),
               tolerance = 1e-6)

})

set.seed(2025)
n <- 40
p <- 5

# Generate random bigger X
X_big <- matrix(rnorm(n * p), nrow = n, ncol = p)
Xtilde_big <- scale(X_big, center = TRUE, scale = TRUE)

# Generate centered Y (no scale)
Y_big <- rnorm(n)
Ytilde_big <- as.numeric(scale(Y_big, center = TRUE, scale = FALSE))

lambda_big <- 0.2

# Initial beta
beta_start0 <- rep(0, p)

test_that("fitLASSOstandardized_c matches R implementation on random standardized data.", {
  
  # R version
  fit_r <- fitLASSOstandardized(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda = lambda_big,
    beta_start = beta_start0,
    eps = 1e-4
  )
  beta_r  <- fit_r$beta

  # Rccp Arma. version
  beta_cpp <- fitLASSOstandardized_c(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda = lambda_big,
    beta_start = beta_start0,
    eps = 1e-4
  )
  
  # 1) Coefficients should match
  expect_equal(as.numeric(beta_cpp),
               as.numeric(beta_r),
               tolerance = 1e-6)
})

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

m3_1 <- microbenchmark(
  fitLASSOstandardized(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda = lambda_small,
    beta_start = NULL,
    eps = 1e-4
  ),
  beta_cpp <- fitLASSOstandardized_c(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda = lambda_small,
    beta_start = numeric(0),
    eps = 1e-4),
  times = 1000)

m3_1

m3_2 <- microbenchmark(
  fitLASSOstandardized(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda = lambda_big,
    beta_start = NULL,
    eps = 1e-4
  ),
  beta_cpp <- fitLASSOstandardized_c(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda = lambda_big,
    beta_start = numeric(0),
    eps = 1e-4),
  times = 1000)

m3_2

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

lambda_seq <- c(1.0, 0.5, 0.1)

test_that("fitLASSOstandardized_seq_c matches R version on the small dataset", {

  # R version
  fit_r <- fitLASSOstandardized_seq(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda_seq = lambda_seq,
    n_lambda = NULL,    # ignored because we supply lambda_seq
    eps = 1e-4
  )
  beta_mat_r <- fit_r$beta_mat   # p x length(lambda_seq)
  
  # Rcpp Arma. version
  beta_mat_cpp <- fitLASSOstandardized_seq_c(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda_seq = lambda_seq,
    eps = 1e-4
  )
  
  # Compare beta's matrices
  expect_equal(as.numeric(beta_mat_cpp),
               as.numeric(beta_mat_r),
               tolerance = 1e-6)
})

test_that("fitLASSOstandardized_seq_c matches R version on the normal standardized data.", {
  
  # R version
  fit_r <- fitLASSOstandardized_seq(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda_seq = lambda_seq,
    n_lambda = NULL,    # ignored because we supply lambda_seq
    eps = 1e-4
  )
  beta_mat_r <- fit_r$beta_mat   # p x length(lambda_seq)
  
  # Rcpp Arma. version
  beta_mat_cpp <- fitLASSOstandardized_seq_c(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda_seq = lambda_seq,
    eps = 1e-4
  )
  
  # Compare beta's matrices
  expect_equal(as.numeric(beta_mat_cpp),
               as.numeric(beta_mat_r),
               tolerance = 1e-6)
})

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################


m4_1 <- microbenchmark(
  fitLASSOstandardized_seq(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda_seq = lambda_seq,
    n_lambda = NULL,    # ignored because we supply lambda_seq
    eps = 1e-4
  ),
  fitLASSOstandardized_seq_c(
    Xtilde = Xtilde_small,
    Ytilde = Ytilde_small,
    lambda_seq = lambda_seq,
    eps = 1e-4
  ),
  times = 1000)

m4_1

m4_2 <- microbenchmark(
  fitLASSOstandardized_seq(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda_seq = lambda_seq,
    n_lambda = NULL,    # ignored because we supply lambda_seq
    eps = 1e-4
  ),
  fitLASSOstandardized_seq_c(
    Xtilde = Xtilde_big,
    Ytilde = Ytilde_big,
    lambda_seq = lambda_seq,
    eps = 1e-4
  ),
  times = 1000)

m4_2

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

