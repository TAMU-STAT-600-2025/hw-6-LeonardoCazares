# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y) {
  # Check all arrays are matrices
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  # Get dimentions of data
  n <- nrow(X)
  p <- ncol(X)
  
  ## [ToDo] Center Y
  Ymean <- mean(Y) # Mean on the vector Y
  Ytilde <- Y - Ymean # Centered vector
  
  ## [ToDo] Center and scale X
  Xmeans <- colMeans(X) # Vector of means for all the columns
  Xtilde <- X - matrix(Xmeans, n, p, byrow = TRUE) # Centered X
  
  weights <- sqrt(colSums(Xtilde^2) / n) # Sample variance for all the columns
  Xtilde <- Xtilde %*% diag(1 / weights) # Normalized X
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(
    Xtilde = Xtilde,
    Ytilde = Ytilde,
    Ymean = Ymean,
    Xmeans = Xmeans,
    weights = weights
  ))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda) {
  return(sign(a) * max(abs(a) - lambda, 0)) # Soft-threshholding value
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda) {
  # Check all the arrays are numeric
  Xtilde <- as.matrix(Xtilde)
  Ytilde <- as.numeric(Ytilde)
  beta <- as.numeric(beta)
  
  # Get dimentions
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  # Check dimentional compatibility among: Xilde, Ytilde and beta
  if (length(Ytilde) != n) {
    stop('Xtilde and Ytilde are incompatible in dimentions.')
  }
  if (length(beta) != p) {
    stop('Xtilde and beta are incompatible in dimentions.')
  }
  
  r <- Ytilde - as.vector(Xtilde %*% beta)
  loss <- as.numeric(crossprod(r)) / (2 * n)
  
  pen <- as.numeric(lambda * sum(abs(beta)))
  
  return(loss + pen)
  
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde,
                                 Ytilde,
                                 lambda,
                                 beta_start = NULL,
                                 eps = 0.001) {
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  
  # Get dimentions from Xtilde and Ytilde
  Xtilde <- as.matrix(Xtilde)
  Ytilde <- as.numeric(Ytilde)
  nXtilde <- nrow(Xtilde)
  nYtilde <- length(Ytilde)
  p <- ncol(Xtilde)
  
  # Check dimentional compatibility for nXtilde and nYtilde
  if (nXtilde != nYtilde) {
    stop('n is not the same between Xtilde and Ytilde')
  }
  
  #[ToDo]  Check that lambda is non-negative
  if (any(lambda < 0)) {
    # Check lambda sign
    stop('lambda is negative')
  }
  
  #[ToDo]  Check for starting point beta_start.
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if (is.null(beta_start)) {
    beta <- rep(0, p) # Initialize beta of 0's
  } else {
    if (length(beta_start) != p) {
      stop('Initial beta is not compatible with Xtilde') # Check dimentional correctness for beta
    }
    beta <- as.numeric(beta_start)
  }
  
  #[ToDo]  Coordinate-descent implementation.
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and4 not have another iteration
  
  n <- nXtilde # Number of samples
  r <- Ytilde - as.vector(Xtilde %*% beta) # Current residual
  # Objective function value for the initial value of beta
  fmin_old <- lasso(Xtilde, Ytilde, beta, lambda) # obj_val(r, beta)
  
  # Number of maximum iterations for the coordinate training
  max_iter <- 1000L
  it <- 0L # Current iteration
  
  # Repeat the coordinate training until the maximum number of iterations has been completed or the threshold has been surpassed
  repeat {
    it <- it + 1L # New iteration
    
    # Compute the coordinate descent for each entry in beta
    for (j in 1:p) {
      xj <- Xtilde[, j] # Current column of the design matrix
      # rho_j uses current residual; with standardization denom is n
      rho_j <- sum(xj * (r + xj * beta[j])) / n # First argument in the soft function
      lamj  <- if (length(lambda) == 1L)
        lambda
      else
        lambda[j] # Check the nature of lambda
      
      bj_new <- soft(rho_j, lamj) # Compute the new value of the j-th entry of beta
      if (bj_new != beta[j]) {
        r <- r - xj * (bj_new - beta[j]) # Update of the residual given the bj_new entry and xj column
        beta[j] <- bj_new # Update the beta
      }
    }
    
    # Update values for the objective function and the difference with the former one.
    fmin <- lasso(Xtilde, Ytilde, beta, lambda)
    # Check if the threshold has been overcomed or the max_iter reached.
    if ((fmin_old - fmin) < eps || it >= max_iter)
      break
    fmin_old <- fmin
  }
  
  # Return
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}


# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde,
                                     Ytilde,
                                     lambda_seq = NULL,
                                     n_lambda = 60,
                                     eps = 0.001) {
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  
  # Check class correctness
  Xtilde <- as.matrix(Xtilde)
  Ytilde <- as.numeric(Ytilde)
  
  # Get Xtilde dimentions and check compatibility
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  if (length(Ytilde) != n) {
    stop('Xtilde and Ytilde have different n')
  }
  
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  
  # If lambda_seq is given
  user_supplied <- !is.null(lambda_seq)
  if (user_supplied) {
    lambda_seq <- as.numeric(lambda_seq) # Check class of lambda_seq
    lambda_seq <- lambda_seq[lambda_seq >= 0] # Keep non-negative values
    lambda_seq <- sort(unique(lambda_seq), decreasing = TRUE) # Decreasing sort
    if (length(lambda_seq) == 0) {
      # If none lambda_seq is given
      warning('None of the supplied values satisfy the requirement')
      user_supplied <- FALSE
    }
  }
  
  
  # If lambda_seq is not supplied, calculate lambda_max
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  if (!user_supplied) {
    # If none lambda_seq is given
    corr <- as.numeric(crossprod(Xtilde, Ytilde)) / n
    lambda_max <- max(abs(corr)) # Compute max lambda
    if (!is.finite(lambda_max) || lambda_max <= 0) {
      lambda_max <- 1e-8 # In case the lambda_max is non-positive or infinite
    }
    lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda)) # Use the given hint
    lambda_seq <- sort(lambda_seq, decreasing = TRUE) # Sort
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda
  # (make sure supplied eps is carried over).
  # Use warm starts strategy discussed in class for setting the starting values.
  L <- length(lambda_seq) # Number of lambdas
  beta_mat <- matrix(0, nrow = p, ncol = L) # Matrix for betas
  fmin_vec <- numeric(L) # Vector for objective function
  
  beta_start <- rep(0, p)  # warm start for the first (largest) lambda
  for (k in seq_len(L)) {
    fit <- fitLASSOstandardized(
      Xtilde,
      Ytilde,
      lambda = lambda_seq[k],
      beta_start = beta_start,
      eps = eps
    ) # Apply fitLASSOstandardized for different lambdas
    beta_mat[, k] <- fit$beta # Keep optimized beta
    fmin_vec[k] <- fit$fmin # Keep current objective function
    beta_start <- fit$beta  # Warm start for the next lambda
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(
    lambda_seq = lambda_seq,
    beta_mat = beta_mat,
    fmin_vec = fmin_vec
  ))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,
                     Y,
                     lambda_seq = NULL,
                     n_lambda = 60,
                     eps = 0.001) {
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  
  std <- standardizeXY(X, Y) # Standardize data
  Xtilde  <- std$Xtilde # Centered and normalized X
  Ytilde  <- std$Ytilde # Centerd Y
  Xmeans  <- std$Xmeans # Means for X
  Ymean   <- std$Ymean # Mean for Y
  weights <- std$weights # Sample Std for each colum of centered X
  p <- ncol(Xtilde) # Number of features in data
  
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  
  # Apply fitLASSOstandardized_seq to the data
  fit_seq <- fitLASSOstandardized_seq(
    Xtilde = Xtilde,
    Ytilde = Ytilde,
    lambda_seq = lambda_seq,
    n_lambda = n_lambda,
    eps = eps
  )
  lambda_seq <- fit_seq$lambda_seq # Sequence of lambdas
  beta_std <- fit_seq$beta_mat # p x L, coefficients for standardized X
  
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  
  L <- length(lambda_seq) # Number of lambda functions
  beta_mat <- matrix(0, nrow = p, ncol = L) # A matrix for each beta
  nz <- weights > 0 # Keep cases where
  if (any(nz)) {
    beta_mat[nz, ] <- sweep(beta_std[nz, , drop = FALSE], 1, weights[nz], "/")
  }
  beta0_vec <- as.numeric(Ymean - crossprod(Xmeans, beta_mat))
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(
    lambda_seq = lambda_seq,
    beta_mat = beta_mat,
    beta0_vec = beta0_vec
  ))
}


###### FALTA COMENTAR
# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,
                    Y,
                    lambda_seq = NULL,
                    n_lambda = 60,
                    k = 5,
                    fold_ids = NULL,
                    eps = 0.001) {
  # [ToDo] Fit Lasso on original data using fitLASSO
  
  # fitLASSO to the original data
  full_fit <- fitLASSO(X,
                       Y,
                       lambda_seq = lambda_seq,
                       n_lambda = n_lambda,
                       eps = eps)
  lambda_seq <- full_fit$lambda_seq
  beta_mat <- full_fit$beta_mat
  beta0_vec <- full_fit$beta0_vec
  
  n <- nrow(as.matrix(X)) # Number of samples
  p <- ncol(as.matrix(X)) # Number of features
  L <- length(lambda_seq) # Number of lambdas
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if (is.null(fold_ids)) {
    fold_ids <- sample(rep(seq_len(k), length.out = n)) # If fold_ids is NULL create random k-fold indexes
  }
  else {
    if (length(fold_ids) != n) {
      stop("fold_ids must have length n.") # Not all samples can be covered within at least one k-fold
    }
    fold_ids <- as.integer(fold_ids) # Verify numeric folds
    k <- max(fold_ids) # Get number of folds (k folds)
  }
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  mse_mat <- matrix(NA_real_, nrow = L, ncol = k)  # MSE matrix for each L value (row) and k value (column).
  
  for (fold in seq_len(k)) {
    idx_val <- which(fold_ids == fold) # Get the indexes that belongs to the the current fold
    idx_tr  <- setdiff(seq_len(n), idx_val) # Get the rest of indexes for training
    
    # Train fitLASSO on the remaining data
    fit_k <- fitLASSO(X[idx_tr, , drop = FALSE],
                      Y[idx_tr],
                      lambda_seq = lambda_seq,
                      n_lambda = n_lambda,
                      eps = eps)
    
    Xv <- as.matrix(X[idx_val, , drop = FALSE]) # Test data X
    yv <- as.numeric(Y[idx_val]) # Test data Y
    pred <- Xv %*% fit_k$beta_mat # Prediction on the test data
    pred <- sweep(pred, 2, fit_k$beta0_vec, "+") # Add intercept terms
    resid <- matrix(yv, nrow = length(yv), ncol = L) - pred # Differences between preds. and ground truth
    mse_mat[, fold] <- colMeans(resid^2) # MSE per lambda for this fold
  }
  
  cvm  <- rowMeans(mse_mat) # CV(lambda), i.e. average MSE over folds
  cvse <- apply(mse_mat, 1, sd) / sqrt(k) # SE_CV(lambda), i.a. std among CV(lambda)
  
  # [ToDo] Find lambda_min
  idx_min <- which.min(cvm)
  lambda_min <- lambda_seq[idx_min] # Choose the lambda with the smalles CV(lambda)
  
  # [ToDo] Find lambda_1SE
  thresh <- cvm[idx_min] + cvse[idx_min] # Mean CV(\lambda_{min}) plus the it's std.
  candidates <- which(cvm <= thresh) # Lambdas with a smaller CV(lambda)
  # lambda_seq is sorted from largest to smallest; pick the largest lambda within 1SE
  lambda_1se <- lambda_seq[min(candidates)]
  
  ## Return output
  ## Output from fitLASSO on the whole data
  ## lambda_seq - the actual sequence of tuning parameters used
  ## beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  ## beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  ## fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  ## lambda_min - selected lambda based on minimal rule
  ## lambda_1se - selected lambda based on 1SE rule
  ## cvm - values of CV(lambda) for each lambda
  ## cvse - values of SE_CV(lambda) for each lambda
  return(
    list(
      lambda_seq = lambda_seq,
      beta_mat = beta_mat,
      beta0_vec = beta0_vec,
      fold_ids = fold_ids,
      lambda_min = lambda_min,
      lambda_1se = lambda_1se,
      cvm = cvm,
      cvse = cvse
    )
  )
}

