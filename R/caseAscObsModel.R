#' @name caseAscObsModel
#' @title Case Ascertainment observational model for a discrete-time meta-population SIR model.
#' @description 
#' Generates a case ascertainment observation model. There is assumed
#' to be a fixed probability \alpha of observing an infection, independent
#' of time and individual.
#'
#' @param X Underlying epidemic trajectory
#' @return 
#' Returns a sampling function and a log-density calculation function.
#' @export 
caseAscObsModel <- function(X){
  # Sampling function for Case Ascertation regime.
  #@param epidemic The output of the discrete-time, meta-population epidemic model
  #@return Returns a cases ascertation sample from the given epidemic

  dimX <- dim(X)
  noDays <- dimX[3] - 1
  noMetapops <- dimX[1]

  # Calculate Daily infections
  I <- matrix(nrow = noMetapops, ncol = noDays)
  for(t in 1:noDays){
    I[,t] <- X[, 1, t] - X[, 1, t + 1]
  }
  sample <- function(alpha){


    # Calculate Daily Infs
    a <- length(I)
    ascCases <- matrix(rbinom(a, I[1:a], prob = alpha),
                       nrow = noMetapops, ncol = noDays, byrow = FALSE)


    # if(sum(ascCases <= I) != a){
    #   stop("Configuration of ascCases differs from epidemic$daily_inf.")
    # }
    #ascCases <- cbind(obs_days, ascCases)

    return(ascCases)
  }

  # @param sample A sample generated from the testing regimes sample function
  # @param alpha  The probability that an individual is recruited for testing
  #               on any given day, regardless of disease status.
  # @return       Returns the log-likelihood of the sample given the epidemic and
  #               observational parameter alpha.
  llh <- function(sample, alpha){
    a <- length(sample)

    #obs_days <- seq(1, a, by = freq)
    #obs_days <- sample[,1]
    #daily_inf <- epidemic$daily_inf[obs_days,]
    #sample <- sample[, -1]
    llh <- sum(dbinom(sample[1:a], I[1:a], prob = alpha, log = T))
    return(llh)
  }

  return(list(sample = sample, llh = llh))
}
