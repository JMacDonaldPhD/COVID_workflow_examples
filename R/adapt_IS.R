#' Adapt Monte Carlo Estimate

#library(foreach)
#library(doParallel)

adapt_IS <- function(N0, IS_estimator, param, no_samples = 1e4, func){
  v <- Inf
  N <- N0
  while(v > 3){
    #registerDoParallel(parallel::detectCores() - 1)
    # est_samples <-
    #   foreach(i = 1:no_samples, .combine = c) %dopar% {
    #     MC_estimator(N, param)
    #   }
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