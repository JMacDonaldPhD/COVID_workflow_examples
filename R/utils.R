# Utilities

#' Transform Uniform Random Variable into a Binomial Random Variable

unif_to_binom <- function(U, trials, prob){
  bool <- U < pbinom(0:trials, size = trials, prob = prob)
  return(min(which(bool) - 1))
}


#' A function which adjusts average likelihood calculation to 
#' avoid numerical error which leads to false zeroes
#' THIS METHOD DOES NOT WORK IN SOME SITUATIONS
#' May make a adjustments such that exp(a_i + b) -> Inf
#' @param a vector of log-likelihood values.
#' @return The log(likelihood-estimate)
log_sum_exp <- function(a){
  N <- length(a)
  if(sum(is.infinite(a)) == N){
    return(-Inf)
  }
  
  a_min <- min(a[!is.infinite(a)])
  logsumexp <- log(sum(exp(a - a_min))) + a_min
  return(logsumexp)
}

sum_exp <- function(a){
  N <- length(a)
  
  a_min <- min(a[!is.infinite(a)])
  sumexp <- sum(exp(a - a_min))
  return(sumexp)
}

log_mean_exp <- function(a){
  N <- length(a)
  
  logsumexp <- log_sum_exp(a)
  
  return(logsumexp - log(N))
}
