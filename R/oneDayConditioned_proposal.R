#' Bootstrap Proposals


oneDayConditioned_proposal <- function(y, N_M){
  if(sum(N_M <= 0) > 0){
    stop("Invalid population levels given (population must be > 0)")
  }
  
  # S -> I -> R
  stoch <- matrix(c(-1, 1, 0,
                    0, -1, 1), nrow = 2, ncol = 3, byrow = T)
  
  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.
  
  
  
  # ==== Set up (Only needs to be done once per population) ====
  # Number of Metapopulations
  M <- length(N_M)
  # Total Populaiton Size
  N <- sum(N_M)
  sim <- function(X_t, t, theta){
    
    beta_G <- theta[1]
    beta_L <- theta[2]
    gamma <- theta[3]

    X <- array(dim = c(M, 3, 2))
    X[,,1] <- X_t
    
    mat <- matrix(NA, ncol = M, nrow = M - 1)
    for(i in 1:M){
      mat[,i] <- X_t[-i,2]
    }
    
    globalInfPressure <- beta_G*colSums(mat)/N
    
    localInfPressure <- beta_L*X_t[, 2]/N_M
    
    infProb <- 1 - exp(-(globalInfPressure + localInfPressure))
    
    # S -> I
    noInf <- rbinom_trunc(1, size = X_t[, 1], prob = infProb, a = y[, t] - 1)
    
    # I -> R
    noRem <- rbinom(M, size = X_t[, 2], prob = 1 - exp(-gamma))
    
    # Update Household States
    X_t <- X_t +
      noInf*matrix(stoch[1, ], nrow = M, ncol = 3, byrow = T) +
      noRem*matrix(stoch[2, ], nrow = M, ncol = 3, byrow = T)
    
    X[,,2] <- X_t
    
    return(X)
  }
  
  llh <- function(X, t, theta){
    beta_G <- theta[1]
    beta_L <- theta[2]
    gamma <- theta[3]
    
    mat <- matrix(NA, ncol = M, nrow = M - 1)
    for(i in 1:M){
      mat[,i] <- X[-i,2, 1]
    }
    
    
    globalInfPressure <- beta_G*colSums(mat)/N
    
    localInfPressure <- beta_L*X[, 2, 1]/N_M
    
    infProb <- 1 - exp(-(globalInfPressure + localInfPressure))
    noInf <- X[, 1, 1] - X[, 1, 2]
    noRem <- X[, 3, 2] - X[, 3, 1] 
    # x = #Susceptibles t_1 - #Susceptibles t_2
    llh <- sum(dbinom_trunc(noInf, size = X[, 1, 1], prob = infProb, a = y[, t] - 1, log = T))
    
    # x = #Removed t_2 - #Removed t_1
    # x = #Infected t_1 - #Infected t_2 - noInfections
    llh <- llh + sum(dbinom(noRem, size = X[, 2, 1], prob = 1 - exp(-gamma), log = T))
    return(llh)
  }
  
  return(list(sim = sim, llh = llh))
  
}