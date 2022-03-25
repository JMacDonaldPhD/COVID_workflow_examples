# Adapt Pseudo-Marginal MCMC

# Adaptive Scheme for lambda and V

adapt_PsMCMC <- function(init, epiModel, obsFrame, epiSample, I0, alpha, logPrior, lambda0, V0, delta = 0.05, 
                         N, noIts){

  # Set up functions
  IS_estimator <- ImportanceSampler(epiSample, epiModel, obsFrame, parallel = F) # Likelihood estimator
  k <- length(init)
  lambda <- lambda0
  # Estimate Likelihood for initial parameters
  logLikeCurr <- -Inf
  while(is.infinite(logLikeCurr)){
    logLikeCurr <- IS_estimator(N, list(I0, init, alpha), func = log_mean_exp)$estimates
  }
  curr <- init
  accept <- 0
  draws <- matrix(ncol = k + 1, nrow = noIts)
  for(i in 1:noIts){
    
    adapt <- accept > 10 & (runif(1, 0, 1) > delta)
    # Propose new parameter
    
    if(adapt){
      V <- var(draws[1:(i-1), 1:k])
      prop <- abs(curr + mvtnorm::rmvnorm(1, mean = rep(0, k), sigma = lambda*V))
    } else{
      prop <- abs(curr + mvtnorm::rmvnorm(1, mean = rep(0, k), sigma = lambda0*V0))
    }

    # Estimate Likelihood
    logLikeProp <- IS_estimator(N, list(I0, prop, alpha), func = log_mean_exp)$estimates
    
    if(!is.infinite(logLikeProp)){
      logAccProb <- (logLikeProp + logPrior(prop)) - (logLikeCurr + logPrior(curr)) 
      #print(logAccPRob)
      if(log(runif(1, 0, 1)) < logAccProb){
        curr <- prop
        logLikeCurr <- logLikeProp 
        accept <- accept + 1
        
        if(adapt){
          lambda <- lambda + 0.93*(lambda/sqrt(i))
        }
      } else{
        if(adapt){
          lambda <- lambda - 0.07*(lambda/sqrt(i))
        }
      }
    } else{
      if(adapt){
        lambda <- lambda - 0.07*(lambda/sqrt(i))
      }
    }
    draws[i, ] <- c(curr, logLikeCurr)
    #print(lambda)
  }
  return(list(lambda = lambda, V = var(draws[,1:k])))
}

