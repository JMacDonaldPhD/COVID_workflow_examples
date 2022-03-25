# Adapt Bootstrap Particle Filter

Adapt_BS_PF <- function(K0 = 10, particleFilter,theta, alpha){
  
  K <- K0
  varLogEst <- Inf
  while(varLogEst > 1){
    varLogEst <- var(replicate(1000, particleFilter(K, theta, alpha)))
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