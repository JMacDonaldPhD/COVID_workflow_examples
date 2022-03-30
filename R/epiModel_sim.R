#' @name sim
#' @title epiModel$sim
#'
#' 
#' @description 
#' An output from an epidemic model generator, used to simulate from the
#' defined epidemic model. This help page serves only to give the general
#' behaviour of such a function, and to show the typical inputs/outputs. As
#' such, behaviour may vary depending on the epidemic model generator used.
#' 
#' @param theta 
#' The epidemic parameters. Structure and shape of this object will
#' depend on the epidemic model generated
#' 
#' @param X0
#' The initial state of the population at the assumed starting point
#' of the simulation. Structure and shape of the state object will
#' depend on the epidemic model generated (a trivial difference occurs
#' when look at individual level model vs. population level models)
#' 
#' @return 
#' Typically, this function will return the states of the epidemic at
#' all simulated time points or an equivalent representation of the
#' trajectory of the epidemic.
#' 
NULL