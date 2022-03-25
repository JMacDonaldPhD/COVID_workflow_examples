#' Revised Metapopulation SIR Model

#' Returns functions to simulate or calculate the probability of an
#' Metapopulation SIR process, given a parameter set.
#' @param N_M The size of each metapopulation. The length of this
#'            vector will determine how many metapopulations there
#'            are.
#' @param startTime An arbitrary baseline for the "start" of the
#'                  epidemic. Defaults to zero.
#' @param endTime   At what time does simulation of the epidemic
#'                  stop. By default, endTime is infinite meaning
#'                  the epidemic is simulated until it dies out.
#' @param PRINT     Do you want to print diagnostic information?
metaSIR_simple <- function(N_M, endTime){

  if(sum(N_M <= 0) > 0){
    stop("Invalid population levels given (population must be > 0)")
  }

  # S -> I -> R
  stoch <- matrix(c(-1, 1, 0,
                    0, -1, 1), nrow = 2, ncol = 3, byrow = T)

  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.



  # ==== Set up (Only needs to be done once per population) ====
  # Number of Metapopulations
  M <- length(N_M)
  # Total Populaiton Size
  N <- sum(N_M)


  #' Progresses the epidemic one day forward and returns the state after
  #' the progression along with how many of each event happened in each
  #' metapopulation.
  #' @param StateX A vector of length three holding information about
  #'               the total number of individuals in each epidemic
  #'               state across all metapopulations.
  #' @param Mstate A matrix of dimension Mx4, where M is the number of
  #'               metapopulations. Each row holds information about
  #'               the number of individuals in each epidemic state in
  #'               a single metapopulation.
  #' @param beta_G Global infection parameter. The rate at which an
  #'               infective from one metapopulation infects a
  #'               susceptible in another metapopulation
  #' @param beta_L Local infection parameter. The rate at which an
  #'               infective from one metapopulation infects a
  #'               susceptible from the same metapopulation.
  #' @param gamma  Removal rate parameter. The rate at which one
  #'               moves from the infective state to the removal
  #'               state.
  #' @return       Returns the propogate StateX and Mstate along
  #'               with the number of infections and removals
  #'               which occured in each metapopulation
  #'               (Infections, Removals)
  dailyProg <- function(X_t, beta_G, beta_L, gamma){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual
    X <- array(dim = c(M, 3, 2))
    X[,,1] <- X_t
    mat <- matrix(NA, ncol = M, nrow = M - 1)
    for(i in 1:M){
      mat[,i] <- X_t[-i,2]
    }
    globalInfPressure <- beta_G*colSums(mat)/N

    localInfPressure <- beta_L*X_t[, 2]/N_M

    infProb <- 1 - exp(-(globalInfPressure + localInfPressure))

    # S -> I
    noInf <- rbinom(M, size = X_t[, 1], prob = infProb)

    # I -> R
    noRem <- rbinom(M, size = X_t[, 2], prob = 1 - exp(-gamma))

    # Update Household States
    X_t <- X_t +
      noInf*matrix(stoch[1, ], nrow = M, ncol = 3, byrow = T) +
      noRem*matrix(stoch[2, ], nrow = M, ncol = 3, byrow = T)
    X[,,2] <- X_t
    return(X)
  }

  dailyProg <- compiler::cmpfun(dailyProg)
  #debug(sim)

  #' Simulate an SIR epidemic within a metapopulation structure,
  #' given initial number of infectives in each metapopulation
  #' and the epidemic parameters.
  #' @param I0    The number of infectives in each metapopulation
  #'              at startTime (see `metaSIR`). Length should
  #'              match the number of metapopulations set at
  #'              model generation.
  #' @param theta The epidemic parameters which are used in
  #'              `dailyProg` to draw the number of infections
  #'              and removals. Three parameters need to be
  #'              given corresponding to beta_G, beta_L and
  #'              gamma in that order (see `dailyProg`).
  #' @return      Returns three objects. epidemic is a list
  #'              object with information about the epidemic
  #'              state at every simulated timepoint. daily_inf
  #'              is a TxM matrix (M is number of metapopulations,
  #'              T is the number of timepoints simulated). Each
  #'              row gives information about the number of infections
  #'              in each metapopulation at a given timepoint.
  #'              daily_rem gives the same information but about
  #'              removals.
  sim <- function(param){
    I0 <- param[[1]]
    theta <- param[[2]]
    if(sum(I0) == 0){
      stop("No epidemic present! Check I0!")
    }
    if(length(I0) != M){
      stop("Information on initial number of infectives must be given for all Meta-populations!")
    }
    if(length(theta) != 3){
      stop("Not enough parameter values given, must be 3. (beta_G, beta_L, gamma)")
    }


    X_0 <- cbind(N_M - I0, I0, rep(0, M))
    X_0 <- unname(X_0)
    X_t <- X_0

    beta_G <- theta[1]
    beta_L <- theta[2]
    gamma <- theta[3]

    X <- array(dim = c(M, 3, endTime + 1))

    time <- 0
    X[,,time + 1] <- X_0
    while(time < endTime){

      # ==== Daily Contact ====
      #debug(dailyProg)
      dayProgression <- dailyProg(X_t, beta_G, beta_L,
                                  gamma)
      #StateX <- dayProgression$StateX
      X_t <- dayProgression[,,2]
      time <- time + 1

      # ==== Store Epidemic Information ====
      X[,,time + 1] <- X_t
      if(sum(X_t[,2]) == 0 & time < endTime){
        for(i in (time + 2):(endTime + 1)){
          X[,,i] <- X_t
        }
        time <- endTime
      }
    }
    return(X)
  }



  #' Calculate the log-likelihood of an SIR epidemic which
  #' is assumed to take place in a metapopulation structure
  #' occurs, given a parameter set.
  #' @param epidemic An epidemic realisation which is in the format
  #'                 returned by `sim`.
  #' @param theta    Epidemic parameters. Three parameters need to be
  #'                 given corresponding to beta_G, beta_L and
  #'                 gamma in that order (see `dailyProg`).
  #' @return         Returns log-likelihood value corresponding to
  #'                 the given epidemic and parameter set.
  llh <- function(X_sample, theta){
    # Daily llh

    ## Daily llh per metapopulation
    one_day_one_metaPop <- function(i, n){

      # Relevant data
      X_1 <- X_sample[, , n]
      X_2 <- X_sample[, , n + 1]

      # Infections and Removals
      inf <- X_1[i, 1] - X_2[i, 1]
      rem <- X_2[i, 3] - X_1[i, 3]

      # Probability of infections
      beta_i <- (theta[1] * sum(X_1[-i, 2]/N_M[-i])) + (theta[2] * X_1[i, 2])/N
      p_inf <- 1 - exp(-beta_i)
      d_inf <- dbinom(inf, X_1[i, 1], prob = p_inf, log = T)

      # Probability of removals
      p_rem <- 1 - exp(-theta[3])
      d_rem <- dbinom(rem, X_1[i, 2], prob = p_rem, log = T)

      return(d_inf + d_rem)
    }

    M_ind <- 1:M
    day_ind <- 1:(dim(X_sample)[3] - 1)
    M_day_grid <- expand.grid(M_ind, day_ind)
    #debug(one_day_one_metaPop)
    llh <- sum(apply(X = M_day_grid, MARGIN = 1, FUN = function(X) one_day_one_metaPop(X[1], X[2])))
    return(llh)
  }

  return(list(sim = sim, llh = llh, param_dim = list(type = "list", dims = c(M, 3)),
              dailyProg = dailyProg))
}
