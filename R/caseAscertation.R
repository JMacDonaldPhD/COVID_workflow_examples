#' Case Ascertation observational model for a discrete-time meta-population SIR model.

#' @param alpha Probability that an infection is detected on the same day that it
#'              occurs.
#' @return Returns a sampling function and a density calculation function.
caseAscObsModel <- function(epidemic){
  #' Sampling function for Case Ascertation regime.
  #'@param epidemic The output of the discrete-time, meta-population epidemic model
  #'@return Returns a cases ascertation sample from the given epidemic
  sample <- function(alpha){
    noDays <- nrow(epidemic$state$S) - 1
    noMetapops <- ncol(epidemic$state$S)
    
    #a <- length(epidemic$daily_inf)
    #obs_days <- seq(from = 1, to = noDays, by = freq)
    #daily_inf <- epidemic$daily_inf[obs_days,]
    a <- length(epidemic$daily_inf)
    ascCases <- matrix(rbinom(a, epidemic$daily_inf[1:a], prob = alpha),
                       nrow = noDays, ncol = noMetapops, byrow = FALSE)
    

    if(sum(ascCases <= epidemic$daily_inf) != a){
      stop("Configuration of ascCases differs from epidemic$daily_inf.")
    }
    #ascCases <- cbind(obs_days, ascCases)
    
    return(ascCases)
  }
  
  #' @param sample A sample generated from the testing regimes sample function
  #' @param alpha  The probability that an individual is recruited for testing
  #'               on any given day, regardless of disease status.
  #' @return       Returns the log-likelihood of the sample given the epidemic and 
  #'               observational parameter alpha.
  llh <- function(sample, alpha){
    a <- length(sample)
    
    #obs_days <- seq(1, a, by = freq)
    #obs_days <- sample[,1]
    #daily_inf <- epidemic$daily_inf[obs_days,]
    #sample <- sample[, -1]
    llh <- sum(dbinom(sample[1:a], epidemic$daily_inf[1:a], prob = alpha, log = T))
    return(llh)
  }
  
  return(list(sample = sample, llh = llh))
}