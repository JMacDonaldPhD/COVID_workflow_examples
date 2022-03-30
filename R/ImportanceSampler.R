#' @name IS
#' @title Importance Sampler
#' @description 
#' Generate a function which calculates Monte Carlo estimate of Epidemic sample Likelihood 
#' given a set of epidemic parameters. 
#' @param sample Observed epidemic data.
#' @param epiModel epidemic model.
#' @param obsFrame Generator function for assumed observational model.
#' @return 
#' An Importance Sampler function.
#' @export
IS <- function(y, obsFrame, epiModel, parallel = F){
  sample_w <- function(param){
    alpha <- param[[3]]
    # EpiParam Structure Check
    # type_check <- typeof(epiParam) == epiModel$param_dim$type
    # dim_check <- prod(sapply(X = 1:length(epiParam), 
    #                         function(X) length(epiParam[[X]]) == epiModel$param_dim$dims[[X]]))
    # if(!type_check | !type_check){
    #   stop("Check Epidemic Parameter Structure!")
    # }
    param <- list(param[[1]], param[[2]])
    
    X <- epiModel$sim(param)
    obsModel <- obsFrame(X)
    logw <- obsModel$llh(y, alpha)
    return(logw)
  }
  
  # ESS <- function(log_weights){
  #   
  # }
  
  est <- function(N, param){
    logw <- replicate(n = N, sample_w(param))
    
    #est <-  lapply(X = func, FUN = function(X) return(X(w)))
    est <- log_mean_exp(logw)
    #w <- exp(logw)
    #w_norm <- w/sum(w)
    #ESS <- 1/(w_norm)^2
    
    return(list(estimates = est, log_weights = logw))
  }

  # if(parallel){
  #   est <- function(N, param){
  #     registerDoParallel(parallel::detectCores() - 1)
  #     llh <- foreach(i = 1:N, .combine = c) %dopar%{
  #       sample_llh(param)
  #     }
  #     log_est <- log_sum_exp(llh)
  #     #log_est <- log(mean(exp(llh)))
  #     return(log_est)
  #   }
  # } else{
  #   est <- function(N, param){
  #     llh <- replicate(n = N, sample_llh(param))
  #     log_est <- log_sum_exp(llh)
  #     #log_est <- log(mean(exp(llh)))
  #     return(log_est)
  #   }
  # }

  return(est)
}


