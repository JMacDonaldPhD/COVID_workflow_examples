#' @name testingObsModel
#' @title Case Ascertainment observational model for a discrete-time meta-population SIR model.
#' @description 
#' Generates a testing observation model. There is assumed
#' to be a fixed probability \alpha of an individual being recruited for testing, independent
#' of time and individual. There is then an associated sensitivity and specificity for the test.
#'
#' @param epidemic Underlying epidemic trajectory
#' @return 
#' Returns a sampling function and a log-density calculation function.
#' @export 
testingObsModel <- function(epidemic){
  
  sample <- function(alpha, sens, spec){
    noDays <- nrow(epidemic$state$S)
    noMetapops <- ncol(epidemic$state$S)
    
    day_metapop_no_tests <- function(i, n){
      X_iS <- epidemic$state$S[n, i]
      X_iI <- epidemic$state$I[n, i]
      X_iR <- epidemic$state$R[n, i]
      
      n_T <- rbinom(3, size = c(X_iS, X_iI, X_iR), prob = alpha)
      return(n_T)
    }
    
    day_metapop_pos_test <- function(i, n, n_T){
      T_pos <- c(rbinom(1, size = n_T[1], prob = 1 - spec), rbinom(1, size = n_T[2], prob = sens),
                 rbinom(1, size = n_T[3], prob = 1 - spec))
      return(T_pos)
    }
    
    one_metapop <- function(i){
      days <- 1:noDays
      n_T <- matrix(sapply(X = days, FUN = function(X) day_metapop_no_tests(i, X)), nrow = noDays,
                    ncol = 3, byrow = T)
      # Want dimensions to be noDays x noStates
      # Also test that N_T[1] < epidemic$state$S for all entries
      test1 <- sum(n_T[, 1] <= epidemic$state$S[, i]) == noDays
      test2 <- sum(n_T[, 2] <= epidemic$state$I[, i]) == noDays
      test3 <- sum(n_T[, 3] <= epidemic$state$R[, i]) == noDays
      
      if(!(test1 & test2 & test3)){
        stop("At least one draw not valid. Matrix not formatted correctly")
      }
      
      if(!all.equal(dim(n_T), c(noDays, 3))){
        stop("Dimension of N_T not correct")
      }
      
      mat <- cbind(n_T, days)
      T_pos <- t(apply(X = mat, MARGIN = 1,
                       function(X) day_metapop_pos_test(i, X[length(X)], X[-length(X)])))
      
      return(list(n_T = rowSums(n_T), T_pos = rowSums(T_pos)))
    }
    
    metapops <- 1:noMetapops
    #debug(one_metapop)
    sample <- lapply(X = metapops, FUN = one_metapop)
    
    return(sample)
  }
  
  llh <- function(sample, alpha, sens, spec){
    
    noDays <- nrow(epidemic$state$S)
    noMetapops <- ncol(epidemic$state$S)
    
    # Function for calculating log density of one metapop, one day
    metapop_day <- function(i, n){
      T_pos <- sample[[i]]$T_pos[n]
      n_T <- sample[[i]]$n_T[n]
      X_iI <- epidemic$state$I[n, i]
      X_iS <- epidemic$state$S[n, i]
      X_iR <- epidemic$state$R[n, i]
      x_D_vec <- max(0, n_T - (X_iS + X_iR)):min(X_iI, n_T)
      
      # 1. Probability of test N_T out of epidemic$state
      dn_T <- function(n_T){
        log_d <- dbinom(n_T, size = sum(X_iS, X_iI, X_iR), prob = alpha, log = T)
        # if(is.infinite(log_d)){
        #   print("Number of tests Density Returning -Inf")
        # }
        return(log_d)
      }
      log_dn_T <- dn_T(n_T)
      # 2. Probability of getting n_D tests of diseased individuals
      dx_D <- function(n_T, x_D){
        X_iS <- epidemic$state$S[n, i]
        X_iI <- epidemic$state$I[n, i]
        X_iR <- epidemic$state$R[n, i]
        
        d <- dhyper(x_D, X_iI, X_iS + X_iR, n_T, log = T)
        if(is.infinite(d)){
          print("d:", d)
          print(x_D_vec)
          print(c("Number of Tests:", n_T))
          print(c("Number of infected tested:", x_D))
          print(c("X(S, I, R):", c(X_iS, X_iI, X_iR)))
          stop("Number of diseased tests Density Returning 0")
        }
        
        return(d)
      }
      
      # 3. Probability of getting R_pos positive tests from n_T tests, n_D of which
      #    are tests of diseased individuals.
      dT_pos <- function(T_pos, n_T, x_D){
        # What are the possible number of positive tests that came from diseased
        # individuals
        x_Dpos_vec <- max(0, n_T - x_D):min(x_D, n_T)
        dens <- function(x_Dpos){
          d1 <- dbinom(x_Dpos, x_D, prob = sens, log = F)
          d2 <- dbinom(T_pos - x_Dpos, n_T - x_D, prob = 1 - spec, log = F)
          # if(d1 == 0 | d2 == 0){
          #   print(c("n_T:", n_T))
          #   print(c("n_D:", x_D))
          #   print(c("n_D+:", x_Dpos))
          #   print(c("n_T+:", T_pos))
          # }
          return(d1*d2)
        }
        d <- sum(sapply(X = x_Dpos_vec, FUN = dens))
        return(d)
      }
      #debug(dT_pos)
      log_d <- log_sum_exp(sapply(X = x_D_vec, FUN = function(X) dx_D(n_T, X) + log(dT_pos(T_pos, n_T, X))))
      return(log_dn_T + log_d)
    }
    # Function fro calculation log density of one metapop, all days
    metapop <- function(i){
      days <- 1:noDays
      log_d <- sum(sapply(X = days, FUN = function(X) metapop_day(i, X)))
      return(log_d)
    }
    # Function for calculation of log density of all metapop, all days
    
    # Return]
    metapops <- 1:noMetapops
    llh <- sum(sapply(X = metapops, FUN = metapop))
    return(llh)
  }
  return(list(sample = sample, llh = llh))
}
