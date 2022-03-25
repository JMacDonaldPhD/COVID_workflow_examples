#' Plot MCMC


plotMCMC <- function(draws,theta, expressions){
  par(mfrow = c(ncol(draws), 3), cex = 1.5)
  
  for(i in 1:ncol(draws)){
    plot(draws[,i], type = "l", ylab = expressions[i])
    acf(draws[,i])
    plot(density(draws[,i]), col = "red")
    if(!missing(theta)){
      abline(v = theta[i], col = 'blue', lty = 2)
    }
  }
  

}
