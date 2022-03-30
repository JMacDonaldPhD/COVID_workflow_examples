#' @name adapt_IS
#' @title Importance Sampler Adaptation
#' @description 
#' Adapts the number of particles used in a given Importance
#' sampler to give adequate performance.
#' @param N0 Initial number of particles.
#' @param IS_estimator A function which carries out Importance sampling.
#' @param param A list of epidemic and observational parameters to be passed 
#'              on to `IS_estimator`
#' @param no_samples Number of IS estimates calculated to estimate the variance
#'                   of the IS estimator. 
#'              
#' @return 
#' Returns number of particles which gives adequate variance of the log-estimate.
#' 
#' @export
adapt_IS <- function(N0, IS_estimator, param, no_samples = 1e4){
  v <- Inf
  N <- N0
  while(v > 1){

    est_samples <- replicate(no_samples, IS_estimator(N, param)$estimates)
    v <- var(est_samples)
    if(is.na(v)){
      v <- Inf
      N <- N + 1
    } else{
      if(v > 3){
        N <- ceiling(N*v/sqrt(3))
      } else{
        return(N)
      }
    }
    print(c(N, v))
  }

}