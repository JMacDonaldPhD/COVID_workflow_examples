#' @name BS_PF
#' @title Bootstrap Particle Filter
#' @description 
#' Generates a bootstrap particle filter.
#' @param y Observed epidemic data
#' @param X_0 Initial epidemic state, which is assumed to be known.
#' @param obsFrame Generator function for observational model.
#' @param epiModel epidemic model
#' @return 
#' Returns log-likelihood estimate.
#' 
#' @export
BS_PF <- function(y, X_0, obsFrame, epiModel){
  
  X_t <- X_0
  noDays <- ncol(y)
  particle_placeholder <- array(X_0, dim = c(nrow(X_0), ncol(X_0), noDays + 1))
  #particles <- rep(list(array(X_0, dim = c(nrow(X_0),ncol(X_0), noDays + 1))), K)
  logLikeEst <- 0
  ESS <- c()
  
  propogate <- function(particle, t, theta, alpha){
    X <- epiModel$dailyProg(particle[,,t], theta[1], theta[2], theta[3])
    
    # # Calculate weights of simulation
    # obsModel <- obsFrame(X)
    # logw_star <- obsModel$llh(y[,t], alpha)
    # 
    particle[,,t + 1] <- X[,,2]
    return(particle)
  }
  
  log_weight <- function(particle, t, alpha){
    obsModel <- obsFrame(particle[,,t:(t+1)])
    logw_star <- obsModel$llh(y[,t], alpha)
    return(logw_star)
  }
  
  particleFilter <- function(K, theta, alpha){
    particles <- rep(list(array(X_0, dim = c(nrow(X_0),ncol(X_0), noDays + 1))), K)
    for(t in 1:noDays){
      logw_star <- c()
      # for(k in 1:K){
      #   # Simulate Forward One Day
      #   X <- epiModel$dailyProg(particles[[k]][,,i], theta[1], theta[2], theta[3])
      #   
      #   # Calculate weights of simulation
      #   obsModel <- obsFrame(X)
      #   logw_star[k] <- obsModel$llh(y[,i], alpha)
      #   
      #   particles[[k]][,,i + 1] <- X[,,2]
      #   
      # }
      
      particles <- lapply(particles, FUN = propogate, t = t, theta = theta, alpha = alpha)
      logw_star <- sapply(particles, FUN = log_weight, t = t, alpha = alpha)
      
      
      # Normalise weights
      if(sum(!is.infinite(logw_star)) == 0){
        logLikeEst <- -Inf
        break
      }
      logw_star_min <- min(logw_star[!is.infinite(logw_star)])
      logw_star <- logw_star - logw_star_min
      w_star <- exp(logw_star)
      
      logLikeEst <- logLikeEst + (log(mean(w_star)) + logw_star_min)
      
      w <- w_star/sum(w_star)
      # if(sum(is.na(w)) > 0){
      #   print(theta)
      # }
      
      ESS[t] <- 1/sum((w^2))
      
      if(t != noDays){
        # Resample (IF K IS LARGE MAYBE QUICKER TO USE `cumsum()` AND `runif()` INSTEAD)
        resample_ind <- sample(1:K, size = K, replace = T, prob = w)
        particles <- particles[resample_ind]
      }
    }
    #return(list(logLikeEst = logLikeEst, ESS = ESS, particles = particles))
    return(list(logLikeEst = logLikeEst, ESS = ESS))
  }
  
  return(particleFilter)

}



