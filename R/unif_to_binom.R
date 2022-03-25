#' Transform Uniform Random Variable into a Binomial Random Variable

unif_to_binom <- function(U, trials, prob){
  bool <- U < pbinom(0:trials, size = trials, prob = prob)
  return(min(which(bool) - 1))
}
