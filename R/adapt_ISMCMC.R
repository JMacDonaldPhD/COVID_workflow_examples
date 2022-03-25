#' Adapt RWM/IS MCMC



adapt_ISMCMC <- function(init, U_init, epiModel, obsFrame, epiSample, I0, alpha, logPrior, 
                         lambda0, V0, S0, noIts, delta = 0.05){
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
  
  S <- S0
  lambda <- lambda0
  
  accept <- 0
  draws <- matrix(ncol = k + 1, nrow = noIts)
  for(i in 1:noIts){
    adapt <- accept > 10 & (runif(1, 0, 1) > delta)
    
    # Propose new epidemic parameters
    if(adapt){
      V <- var(draws[1:(i-1), 1:k])
      prop <- abs(curr + mvnfast::rmvn(1, mu = rep(0, k), sigma = lambda*V))
    } else{
      prop <- abs(curr + mvtnorm::rmvnorm(1, mean = rep(0, k), sigma = lambda0*V0))
    }

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
        accept <- accept + 1
        if(adapt){
          lambda <- lambda + 3*(lambda/i^(1/3))
        } 
      } else{
        if(adapt){
          lambda <- lambda - 1*(lambda/i^(1/3))
        }
      }
    }
    

    # Propose new epidemic parameters
    if(adapt){
      U_prop <- U_curr
      U_prop[sample(1:r, S, replace = F)] <- runif(S, 0, 1)
    } else{
      U_prop <- U_curr
      U_prop[sample(1:r, S0, replace = F)] <- runif(S0, 0, 1)
    }
    
    # Estimate Likelihood
    epiProp <- epiModel$sim(list(I0, curr, U_prop))
    obsModelProp <- obsFrame(epiProp)
    logLikeProp <- obsModelProp$llh(epiSample, alpha)
    
    
    if(!is.infinite(logLikeProp)){
      logAccProb <- (logLikeProp + logPrior(prop)) - (logLikeCurr + logPrior(curr)) 
      #print(logAccPRob)
      if(log(runif(1, 0, 1)) < logAccProb){
        U_curr <- U_prop
        logLikeCurr <- logLikeProp 
        accept <- accept + 1
        if(adapt){
          S <- S + rbinom(1, 3, prob = 0.75/sqrt(i))
          #ceiling(0.75*(S/sqrt(i)))
          if(S > r){
            S <- r
          }
        } 
      } else{
        if(adapt){
          S <- S - rbinom(1, 1, prob = 0.25/sqrt(i))
          #ceiling(0.25*(S/sqrt(i)))
          if(S < 1){
            S <- 1
          }
        }
      }
    }
    draws[i, ] <- c(curr, logLikeCurr)
  }
  return(list(lambda = lambda, V = var(draws[,1:k]), S = S))
}