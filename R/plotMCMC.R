#' @name plotMCMC
#' @title Plot MCMC Samples and Diagnostics
#' @description 
#' Takes MCMC samples and plots traceplots, ACFs and density
#' graphs for each parameter.
#' @param draws Matrix of MCMC samples
#' @param theta True values of parameters, if epidemic data is simulated.
#' @param expressions expression objects which represent mathematical notation
#'                    for each parameter.
#' @return 
#' None
#' @export
plotMCMC <- function(draws, theta, expressions){
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
