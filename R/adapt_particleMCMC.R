#' Adapt Particle MCMC

adapt_particleMCMC <- function(init, epiModel, obsFrame, y, X0, alpha, logPrior, lambda0, V0, 
                               delta = 0.05, K, noIts){
  
  # Set up functions
  particleFilter <- BS_PF(y, X0, obsFrame, epiModel)
  k <- length(init)
  lambda <- lambda0
  # Estimate Likelihood for initial parameters
  logLikeCurr <- -Inf
  while(is.infinite(logLikeCurr)){
    logLikeCurr <- particleFilter(K, init, alpha)
  }
  curr <- init
  accept <- 0
  lambda_vec <- lambda
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
    logLikeProp <- particleFilter(K, prop, alpha)

    if(!is.infinite(logLikeProp)){
      logAccProb <- (logLikeProp + logPrior(prop)) - (logLikeCurr + logPrior(curr)) 
      #print(logAccPRob)
      if(log(runif(1, 0, 1)) < logAccProb){
        curr <- prop
        logLikeCurr <- logLikeProp 
        accept <- accept + 1
        
        if(adapt){
          lambda <- lambda + 0.90*(lambda/sqrt(i))
        }
      } else{
        if(adapt){
          lambda <- lambda - 0.10*(lambda/sqrt(i))
        }
      }
    } else{
      if(adapt){
        lambda <- lambda - 0.10*(lambda/sqrt(i))
      }
    }
    draws[i, ] <- c(curr, logLikeCurr)
    lambda_vec[i + 1] <- lambda
    #print(lambda)
  }
  par(mfrow = c(1,1))
  plot(lambda_vec, type = 'l')
  return(list(lambda = lambda, lambda_vec = lambda_vec, V = var(draws[,1:k])))
}
