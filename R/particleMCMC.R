#' Particle MCMC


particleMCMC <- function(init, epiModel, obsFrame, y, X0, alpha, logPrior, lambda, V, K, noIts){
  # Set up functions
  particleFilter <- BS_PF(y, X0, obsFrame, epiModel)
  k <- length(init)
  # Estimate Likelihood for initial parameters
  logLikeCurr <- -Inf
  while(is.infinite(logLikeCurr)){
    logLikeCurr <- particleFilter(K, init, alpha)
  }
  curr <- init
  accept <- 0
  draws <- matrix(ncol = k + 1, nrow = noIts)
  for(i in 1:noIts){
    # Propose new parameters
    prop <- abs(curr + mvnfast::rmvn(1, mu = rep(0, k), sigma = lambda*V))
    
    # Estimate Likelihood
    logLikeProp <- particleFilter(K, prop, alpha)

    if(!is.infinite(logLikeProp)){
      logAccProb <- (logLikeProp + logPrior(prop)) - (logLikeCurr + logPrior(curr)) 
      #print(logAccPRob)
      if(log(runif(1, 0, 1)) < logAccProb){
        curr <- prop
        logLikeCurr <- logLikeProp 
        accept <- accept + 1
      }
    }
    draws[i, ] <- c(curr, logLikeCurr)
  }
  return(list(draws = draws, acceptRate = accept/noIts))
} 