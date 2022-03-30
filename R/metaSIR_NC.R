
#' @name metaSIR_NC
#' @title Generate a non-centered discrete-time epidemic model
#' @description 
#' Generates a discrete-time epidemic model which constructs epidemics
#' using a stream of random variables which are independent to the 
#' epidemic parameters and state.
#' @param N_M Size of each meta-population.
#' @param endTime The time which the epidemic should be simulated for.
#' 
#' @return 
#' Returns three functions. A `sim` function to simulate from the generated model,
#' given a set of parameters. A `dailyProg` function, whose behaviour is identical
#' to `sim` but simulates for one time unit only. Finally, a `llh` function to calculate
#' he log-density for an epidemic, given a set of parameters.
#' @export
metaSIR_NC <- function(N_M, startTime, endTime, PRINT = F){
  if(sum(N_M <= 0) > 0){
    stop("Invalid population levels given (population must be > 0)")
  }
  
  # S -> I  -> R
  stoch <- matrix(c(-1, 1, 0,
                    0, -1, 1), nrow = 2, ncol = 3, byrow = T)
  if(startTime > 0){
    cutData <- 1:startTime
  }
  
  # projects epidemic one day using binomial draws
  # Either takes P0 as an argument for the intial number of presymptomatic individuals
  # OR StateX if simulating from a timepoint which is mid-epidemic.
  
  
  
  # ==== Set up (Only needs to be done once per population) ====
  # Number of Metapopulations
  M <- length(N_M)
  # Total Populaiton Size
  N <- sum(N_M)
  
  
  dim_U_i <- M*nrow(stoch)
  dim_U <- c(nrow(stoch)*(endTime - startTime), M)
  dailyProg <- function(StateX, Mstate, beta_G, beta_L, gamma, U_i){
    # Calculates the number of infected individuals exerting pressure on each susceptible individual
    mat <- matrix(NA, ncol = M, nrow = M - 1)
    for(i in 1:M){
      mat[,i] <- Mstate[-i,2]
    }
    
    # Calculates the number of infected individuals exerting pressure on each susceptible individual
    globalInfPressure <- beta_G*colSums(mat)/N
    
    localInfPressure <- beta_L*Mstate[, 2]/N
    
    infProb <- 1 - exp(-(globalInfPressure + localInfPressure))
    
    remProb <- rep(1 - exp(-gamma), M)
    
    noInf <- apply(matrix(c(U_i[1, ], Mstate[,1], infProb),  ncol = 3, nrow = M, byrow = F),
                   MARGIN = 1, FUN = function(X) unif_to_binom(X[1], X[2], X[3]))
    
    noRem <- apply(matrix(c(U_i[2, ], Mstate[,2], remProb),  ncol = 3, nrow = M, byrow = F),
                   MARGIN = 1, FUN = function(X) unif_to_binom(X[1], X[2], X[3]))
    
    
    # Update Household States
    Mstate <- Mstate +
      noInf*matrix(stoch[1, ], nrow = M, ncol = 3, byrow = T) +
      noRem*matrix(stoch[2, ], nrow = M, ncol = 3, byrow = T)
    
    # Update population state
    totInf <- sum(noInf)
    totRem <- sum(noRem)
    StateX <- StateX + c(-totInf, totInf - totRem, totRem)
    
    # Population State and Household State Check
    if(!identical(colSums(Mstate), StateX)){
      stop("Population State and Household States do not line up")
    }
    
    return(list(Mstate = Mstate, StateX = StateX, Infections = noInf, Removals = noRem))
    #, logProb = logProb))
  }
  
  dailyProg <- compiler::cmpfun(dailyProg)
  #debug(sim)
  
  sim <- function(param){
    I0 <- param[[1]]
    theta <- param[[2]]
    U <- param[[3]]
    beta_G <- theta[1]
    beta_L <- theta[2]
    gamma <- theta[3]
    
    Mstate <- cbind(N_M - I0, I0, rep(0, M))
    Mstate <- unname(Mstate)
    StateX <- c(N - sum(I0), sum(I0), 0)
    if(!identical(colSums(Mstate), StateX)){
      stop("Population State and Household States do not line up (start)")
    }
    
    # Stores summry of each epidemic state at the end of each day
    SIRsummary <- matrix(nrow = endTime + 1, ncol = 3)
    A <- matrix(nrow = endTime + 1, ncol = M)
    
    time <- 0
    SIRsummary[time + 1, ] <- StateX
    epidemic <- list(S = A, I = A, R = A)
    epidemic$S[time + 1, ] <- Mstate[,1]
    epidemic$I[time + 1, ] <- Mstate[,2]
    epidemic$R[time + 1, ] <- Mstate[,3]
    
    Infections <- A[-1, , drop = F]
    Removals <- A[-1, , drop = F]
    logProb <- 0
    U_indices <- seq(1, endTime*2, by = 2)
    while(time < endTime){
      
      # ==== Daily Contact ====
      U_i <- U[U_indices[time + 1] + c(0, 1), , drop = F]
      #debug(dailyProg)
      # if(time == 5){
      #   debug(dailyProg)
      # }

      dayProgression <- dailyProg(StateX, Mstate, beta_G, beta_L, gamma, U_i)
      StateX <- dayProgression$StateX
      Mstate <- dayProgression$Mstate
      time <- time + 1
      
      # ==== Store Epidemic Information ====
      epidemic$S[time + 1, ] <- Mstate[,1]
      epidemic$I[time + 1, ] <- Mstate[,2]
      epidemic$R[time + 1, ] <- Mstate[,3]
      SIRsummary[time + 1, ] <- dayProgression$StateX
      Infections[time, ] <- dayProgression$Infections
      Removals[time, ] <- dayProgression$Removals
      #logProb <- logProb + dayProgression$logProb

      if(StateX[2] == 0 & time < endTime){
        for(i in (time + 1):(nrow(epidemic$S) - 1)){
          epidemic$S[i + 1, ] <- Mstate[,1]
          epidemic$I[i + 1, ] <- Mstate[,2]
          epidemic$R[i + 1, ] <- Mstate[,3]
          SIRsummary[i + 1, ] <- dayProgression$StateX
          Infections[i, ] <- rep(0, nrow(Mstate))
          Removals[i, ] <- rep(0, nrow(Mstate))
        }
        time <- endTime + 1
      }
    }
    if(PRINT){
      print("Day by day household States:")
      print(epidemic)
      print("Daily Infections by Household:")
      print(Infections)
    }
    
    # Cut Data
    
    if(startTime > 0){
      epidemic$S <- epidemic$S[-cutData, ]
      epidemic$I <- epidemic$I[-cutData, ]
      epidemic$R <- epidemic$R[-cutData, ]
      SIRsummary <- SIRsummary[-cutData, ]
      Infections[-cutData[-length(cutData)], ]
      Removals[-cutData[-length(cutData)], ]
    }
    
    return(list(state = epidemic, daily_inf = Infections,
                daily_rem = Removals))
  }

  llh <- function(epidemic, theta){
    # Daily llh
    
    ## Daily llh per metapopulation
    one_day_one_metaPop <- function(i, n){
      
      X_Si <- epidemic$state$S[n, i]
      X_Ii <- epidemic$state$I[n, i]
      X_Ij <- epidemic$state$I[n, -i]
      
      inf <- epidemic$daily_inf[n, i]
      rem <- epidemic$daily_rem[n, i]
      # Probabilty of infections
      beta_i <- (theta[1] * sum(X_Ij))/N + (theta[2] * X_Ii)/N
      p_inf <- 1 - exp(-beta_i)
      d_inf <- dbinom(inf, X_Si, prob = p_inf, log = T)
      
      # Probability of removals
      p_rem <- 1 - exp(-theta[3])
      d_rem <- dbinom(rem, X_Ii, prob = p_rem, log = T)
      
      return(d_inf + d_rem)
    }
    
    M_ind <- 1:M
    day_ind <- 1:nrow(epidemic$daily_inf)
    M_day_grid <- expand.grid(M_ind, day_ind)
    #debug(one_day_one_metaPop)
    llh <- sum(apply(X = M_day_grid, MARGIN = 1, FUN = function(X) one_day_one_metaPop(X[1], X[2])))
    return(llh)
  }
  return(list(sim = sim, llh = llh))
}

