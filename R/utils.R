# Utilities
# Transform Uniform Random Variable into a Binomial Random Variable
unif_to_binom <- function(U, trials, prob){
  bool <- U < pbinom(0:trials, size = trials, prob = prob)
  return(min(which(bool) - 1))
}


# Calculating product of probabilities in a more stable way
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

# ===== Truncated Binomial Draws ====

rbinom_trunc <- function(n = 1, a, size, prob){
  return(truncdist::rtrunc(n, spec = "binom", a = a, size = size, prob = prob))
}

rbinom_trunc <- Vectorize(rbinom_trunc, vectorize.args = c("a", "size", "prob"))

dbinom_trunc <- function(x, a, size, prob, log = T){
  if(log){
    return(log(truncdist::dtrunc(x, spec = "binom", a = a, size = size, prob = prob)))
  } else{
    return(truncdist::dtrunc(x, spec = "binom", a = a, size = size, prob = prob))
  }
}

dbinom_trunc <- Vectorize(dbinom_trunc, vectorize.args = c("x", "a", "size", "prob"))


# Convert 2 strings into the distribution and function desired
string2dist <- function(fun = "d", dist = NULL){
  
  if(!is.null(dist)){
    if(dist == "AC"){
      dist <- "AscCase"
    } else{
      stop("Provide a valid distribution")
    }
  } else{
    dist <- "AscCase"
  }
  return(eval(parse(text = paste0(c(fun, dist), collapse = ""))))
}

# Change objects within VmetaSIR function factory
amend_K <- function(K, VepiModel){
  epi_env <- rlang::fn_env(VepiModel$sim)
  
  if(epi_env$K != K){
   epi_env$K <- K
   epi_env$stoch <- Matrix::bdiag(rep(list(epi_env$stoch_base), K))
  }
}



