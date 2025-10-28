#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Compute |a| - lambda
  double difference = std::fabs(a) - lambda;
  
  // If difference <= 0 return zero
  if (difference <= 0.0) {
    return 0.0;
  }
  
  // Otherwise return sign(a) * difference if difference > 0
  if (a > 0.0) {
    return difference;
  } else {
    return -difference;
  }

}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
 
  // Get dimentions
  double n = Xtilde.n_rows;
  double p = Xtilde.n_cols;

  // Define differences
  arma::colvec r = Ytilde - Xtilde * beta;

  // Loss term
  double loss = arma::dot(r, r) / (2.0 * n);    
  
  // L1 penalty
  double pen = lambda * arma::sum(arma::abs(beta));
  
  // Lasso objective function
  return loss + pen;  
  
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
}