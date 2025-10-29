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
  
  // Obtain the design matrix dimentions
  double n = Xtilde.n_rows;
  double p = Xtilde.n_cols;

  // Initialize the beta given or not
  arma::colvec beta;
  if (beta_start.n_elem == 0) {
    beta = arma::zeros<arma::colvec>(p);
  } else {
    beta = beta_start;
  }

  // Get the current residual
  arma::colvec r = Ytilde - Xtilde * beta;
  
  // Compute the current lasso objective function
  double fmin_old = lasso_c(Xtilde, Ytilde, beta, lambda);
  
  // Maximum number of iterations for coordinate descent, and the iteration variable
  int max_iter = 1000;
  int it = 0;
  
  // Coordinate-descent
  while (true) {
    // New iteration
    it += 1;
    
    // Loop over coordinates of beta
    for (int j = 0; j < p; ++j) {
      // Current column x_j
      arma::colvec xj = Xtilde.col(j);
      
      // rho_j = sum( x_j * (r + x_j * beta[j]) ) / n
      double rho_j = arma::dot(xj, (r + xj * beta(j))) / static_cast<double>(n);
      
      double lamj = lambda;
      
      // Soft-threshold update for coordinate j
      double bj_new = soft_c(rho_j, lamj);
      
      // If coefficient changed, update residual and beta
      if (bj_new != beta(j)) {
        // Update
        r = r - xj * (bj_new - beta(j));
        beta(j) = bj_new;
      }
    }
    
    // Compute new objective
    double fmin = lasso_c(Xtilde, Ytilde, beta, lambda);
    
    // Stop if improvement is smaller than eps OR max_iter reached
    if ( (fmin_old - fmin) < eps || it >= max_iter ) {
      break;
    }
    
    // Otherwise continue
    fmin_old = fmin;
  }
  
  // Return beta (solution vector)
  return beta; 

}

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  
  // Get dimensions
  int p = Xtilde.n_cols;                
  int L = lambda_seq.n_elem;            
  
  // Matrix to store solutions
  arma::mat beta_mat(p, L, arma::fill::zeros);
  
  // Beta with zeros
  arma::colvec beta_start = arma::zeros(p);
  
  // Loop over each lambda in the sequence
  for (int k = 0; k < L; ++k) {
    double lam = lambda_seq(k);
    
    // Fit LASSO for this lambda using coordinate descent
    arma::colvec beta_k = fitLASSOstandardized_c(
      Xtilde,
      Ytilde,
      lam,
      beta_start,
      eps
    );
    
    // Store current beta in the beta matrix
    beta_mat.col(k) = beta_k;
    
    // Take the last beta as the new starter value
    beta_start = beta_k;
  }
  
  // Return p x L matrix of coefficients
  return beta_mat;
}
  