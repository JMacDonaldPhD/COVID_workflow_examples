#' @name llh
#' @title epiModel$llh
#'
#' @description 
#' An output from an epidemic model generator, used to calculate the log-likelihood
#' of an epidemic trajectory, given the epidemic model defined. This help page serves 
#' only to give the general behaviour of such a function, and to show the typical 
#' inputs/outputs. As such, behaviour may vary depending on the epidemic model generator 
#' used.
#' 
#' @usage 
#' 
#' @param X
#' Epidemic trajectory. Structure and shape of the state object will
#' depend on the epidemic model generated (a trivial difference occurs
#' when look at individual level model vs. population level models)
#' 
#' @param theta 
#' The epidemic parameters. Structure and shape of this object will
#' depend on the epidemic model generated
#' 
#' @return 
#' Typically, this function will return a single value, the log-likelihood
#' of the epidemic trajectory, given the defined epidemic model and a set
#' of epidemic parameter.
#' 
NULL