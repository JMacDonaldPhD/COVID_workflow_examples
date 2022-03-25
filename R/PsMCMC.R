#' Pseudo-Marginal MCMC Sampler for epidemic models

PsMCMC <- function(init, epiModel, obsFrame, epiSample, I0, alpha, logPrior, lambda, V, N, noIts){
  # Set up functions
  IS_estimator <- ImportanceSampler(epiSample, epiModel, obsFrame) # Likelihood estimator
  k <- length(init)
  # Estimate Likelihood for initial parameters
  logLikeCurr <- -Inf
  while(is.infinite(logLikeCurr)){
    logLikeCurr <- IS_estimator(N, list(I0, init, alpha), func = log_mean_exp)$estimates
  }
  curr <- init
  accept <- 0
  draws <- matrix(ncol = k + 1, nrow = noIts)
  for(i in 1:noIts){
    # Propose new parameters
    prop <- abs(curr + mvnfast::rmvn(1, mu = rep(0, k), sigma = lambda*V))
    
    # Estimate Likelihood
    logLikeProp <- IS_estimator(N, list(I0, prop, alpha), log_mean_exp)$estimates
    
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