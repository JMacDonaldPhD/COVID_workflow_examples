#' @name r/dcaseAscObsModel
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
#' 

rAscCases <- function(alpha, X){
  dimX <- dim(X)
  noDays <- dimX[3] - 1
  noMetapops <- dimX[1]
  
  # Calculate Daily infections
  #I <- matrix(nrow = noMetapops, ncol = noDays)
  I <- X[, 1, 2] - X[, 1, 1]
  
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

dAscCase <- function(y, alpha, X){
  # dimX <- dim(X)
  # noDays <- dimX[3] - 1
  # noMetapops <- dimX[1]
  if(is.na(dim(X)[4])){
    I <- X[, 1, 1] - X[, 1, 2]
    a <- length(y)
    llh <- sum(dbinom(y[1:a], I[1:a], prob = alpha, log = T))
    return(llh)
  } else{
    #K <- dim(X)[3]
    K <- dim(X)[4]
    I <- X[, 1, 1, ] - X[, 1, 2, ]
    #I <- X[, 1, , 1] - X[, 1, , 2]
    a <- length(y)
    b <- length(I)
    vec_density <- dbinom(rep(y[1:a], K), I[1:b], prob = alpha, log = T)
    mat_density <- matrix(vec_density, nrow = a, ncol = K, byrow = F)
    #print(dim(mat_density))
    llh <- colSums(mat_density)
    #print(length(llh))
    return(llh)
  }
  # Calculate Daily infections
  #I <- matrix(nrow = noMetapops, ncol = noDays)
  #I <- X[, 1, 1, ] - X[, 1, 2, ]
  #I <- X[, 1, 1] - X[, 1, 2]  
  

  #obs_days <- seq(1, a, by = freq)
  #obs_days <- sample[,1]
  #daily_inf <- epidemic$daily_inf[obs_days,]
  #sample <- sample[, -1]

  #llh <- sum(dbinom(y[1:a], I[1:a], prob = alpha, log = T))
  
  #return(llh)
}