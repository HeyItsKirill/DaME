#' @title Multivariate survival data simulation from Clayton-Oakes model
#' 
#' @description Generates simulated k-variate survival data from a Clayton-Oakes model.
#' The marginals are exponentially distributed with rate parameters 1. 
#' The joint distribution follows from the Clayton copula.
#' Censoring times are exponentially distributed with vector of rate parameters lambdaC. 
#' 
#' @param n Number of observations to generate.
#' @param theta Clayton copula parameter (-1 <= theta < inf).
#' @param lambdaC Vector of exponentially-distributed rate parameters for censoring times.
#' (0 = censored, 1 = uncensored)
#' 
#' @return An n x 2k data frame with the following items:
#' \describe{
#' \item{X1, X2, ..., Xk:}{Survival times for k dimensions.}
#' \item{Delta1, Delta2, ..., Deltak:}{Censoring indicators for k dimensions.}
#' }  
#' 
#' @examples 
#' X <- genClaytonk(200, 3, c(0.1, 0.2, 0.3, 0.4))
#' 
genClaytonk <- function(n, theta, lambdaC){
  
  # Check if the input parameters are valid. 
  
  if(theta < -1){
    stop("theta must be >= -1")
  }
  if(min(lambdaC) < 0){
    stop("lambdaC must be nonnegative")
  }
  
  # Generate random uniform points from the unit hypercube.
  
  u = matrix(runif(n * length(lambdaC)), nrow = n)
  
  # Initialize matrices for failure times (failt) and censoring times (centt).
  
  failt <- matrix(0, nrow = n, ncol = length(lambdaC)); failt[, 1] = -1 * log(runif(n))
  centt <- matrix(0, nrow = n, ncol = length(lambdaC))
  
  # Generate survival data for each variable.
  
  Y = matrix(NA, nrow = n, ncol = length(lambdaC))
  
  for(i in 1:length(lambdaC)){
    
    # Compute the failure times based on Copula parameter cases. 
    
    if (theta == 0){
      failt[, i] <- -1 * log(u[ ,i])
    }
    
    else{
      cum_u = matrix(0, nrow = n, ncol = length(lambdaC))
      
      for (j in 1:length(lambdaC)) { 
        cum_u[, j] <- (u[, j])^((-1 * theta) * (1 + ((j - 1) * theta))^(-1))
      }
      
      if (i > 1) {
        cum_prod <- 1
        for (k in 1:(i - 1)) {
          cum_prod <- cum_prod * cum_u[, k]
        }
        
        failt[, i] <- theta^(-1) * log(1 - cum_prod * (1 - cum_u[, i]))
      }
    }
    
    # Compute censoring times.
    
    if(lambdaC[i] != 0){
      centt[, i] <- rexp(n, lambdaC[i])
    } 
    else{
      centt[, i] <- max(failt[, i]) + 1
    }
    
  }
  
  X = pmin(failt, centt)
  
  Delta <- 1 * (failt < centt)
  
  # Create the data frame of survival time and censoring indicator.
  
  data = data.frame(X, Delta)
  
  colnames(data) <- c(paste("X", 1:length(lambdaC), sep = ""), 
                      paste("Delta", 1:length(lambdaC), sep = ""))
  
  data[] <- lapply(data, function(x) as.numeric(as.character(x)))
  
  
  return(data)
  
}

  
