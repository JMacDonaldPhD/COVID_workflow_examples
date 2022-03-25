#' Importance Sampling MCMC for non-centered epidemic models



ISMCMC <- function(init, U_init, epiModel, obsFrame, epiSample, I0, alpha, logPrior, lambda, V, S, noIts){
  # Set up functions
  k <- length(init)
  r <- length(U_init)
  # Estimate Likelihood for initial parameters
  epi_init <- epiModel$sim(list(I0, init, U_init))
  obsModelInit <- obsFrame(epi_init)
  logLikeCurr <- obsModel$llh(epiSample, alpha)
  while(is.infinite(logLikeCurr)){
    U_init <- matrix(runif(r, 0, 1), ncol = ncol(U_init), nrow = nrow(U_init))
    epi_init <- epiModel$sim(list(I0, init, U_init))
    obsModelInit <- obsFrame(epi_init)
    logLikeCurr <- obsModelInit$llh(epiSample, alpha)
  }
  curr <- init
  U_curr <- U_init
  epi_curr <- epi_init
  obsModelCurr <- obsModelInit
  
  accept.theta <- 0
  accept.S <- 0
  draws <- matrix(ncol = k + 1, nrow = noIts)
  for(i in 1:noIts){
    # Propose new epidemic parameters
    prop <- abs(curr + mvnfast::rmvn(1, mu = rep(0, k), sigma = lambda*V))

    # Estimate Likelihood
    
    epiProp <- epiModel$sim(list(I0, prop, U_curr))
    obsModelProp <- obsFrame(epiProp)
    logLikeProp <- obsModelProp$llh(epiSample, alpha)
    
    
    if(!is.infinite(logLikeProp)){
      logAccProb <- (logLikeProp + logPrior(prop)) - (logLikeCurr + logPrior(curr)) 
      #print(logAccPRob)
      if(log(runif(1, 0, 1)) < logAccProb){
        curr <- prop
        logLikeCurr <- logLikeProp 
        accept.theta <- accept.theta + 1
      }
    }
    
    U_prop <- U_curr
    U_prop[sample(1:r, S, replace = F)] <- runif(S, 0, 1)
    
    # Estimate Likelihood
    epiProp <- epiModel$sim(list(I0, curr, U_prop))
    obsModelProp <- obsFrame(epiProp)
    logLikeProp <- obsModelProp$llh(epiSample, alpha)
    
    if(!is.infinite(logLikeProp)){
      logAccProb <- (logLikeProp) - (logLikeCurr) 
      #print(logAccPRob)
      if(log(runif(1, 0, 1)) < logAccProb){
        U_curr <- U_prop
        logLikeCurr <- logLikeProp 
        accept.S <- accept.S + 1
      }
    }
    draws[i, ] <- c(curr, logLikeCurr)
  }
  return(list(draws = draws, acceptRate = c(accept.theta/noIts, accept.S/noIts)))
}