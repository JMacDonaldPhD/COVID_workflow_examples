#' @name adapt_BS_PF
#' @title Bootstrap Particle Filter Adaptation
#' @description 
#' Adapts the number of particles used in a given bootstrap particle
#' filter to give adequate performance.
#' @param K_0 Initial number of particles.
#' @param particleFilter A function which carries out Particle Filter sampling.
#' @param theta Epidemic parameters to be passed on to `particleFilter`
#' @param alpha Observational parameters to be passed on to `particleFilter`
#' @return 
#' Returns number of particles which gives adequate variance of the log-estimate.
#' 
#' @export
adapt_BS_PF <- function(K0 = 10, particleFilter,theta, alpha){
  
  K <- K0
  varLogEst <- Inf
  while(varLogEst > 1){
    varLogEst <- var(replicate(1000, particleFilter(K, theta, alpha)$logLikeEst))
    if(is.na(varLogEst)){
      varLogEst <- Inf
      K <- K + 10
    } else{
      if(varLogEst > 1){
        K <- ceiling(K*varLogEst)
      } else{
        return(K)
      }
    }
    
    print(c(K, varLogEst))
  }
}