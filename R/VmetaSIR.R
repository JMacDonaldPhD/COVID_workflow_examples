
#' @export
VmetaSIR <- function(N_M, endTime, K){
  N <- sum(N_M)
  M <- length(N_M)
  stoch <- matrix(c(-1, 1, 0, 0, -1, 1), nrow = 2, ncol = 3, byrow = T)
  stoch <- Matrix::bdiag(rep(list(stoch), K))
  
  mat <- matrix(1, nrow = M, ncol = M)
  diag(mat) <- 0
  oneDaySim <- function(X_t, theta){
    
    particles <- X_t
    #particles <- array(rep(X_t, K), dim = c(dim(X_t), K))
    
    #p_t <- array(dim = c(dim(X0), K, 2))
    #p_t[, , ,1] <- particles

    globalInfectiousPressure <- theta[1]*mat%*%particles[,2,]/N # matrix MxK
    localInfectiousPressure <- theta[2]*particles[,2,]/N_M
    infProb <- 1 - exp(-(globalInfectiousPressure + localInfectiousPressure))
    
    
    probs <- rbind(infProb, matrix(1-exp(-theta[3]), nrow = M, ncol = K))
    
    events <- matrix(rbinom(2*M*K, size = particles[,1:2,][1:(2*M*K)], prob = probs[1:(2*M*K)]), 
                     nrow = M, ncol = 2*K, byrow = F)
    delta <- events%*%stoch
    particles <- particles + array(delta, dim = c(M, 3, K))
    #p_t[,,,2] <- particles
    return(particles)
  }
  
  llh <- function(X, theta){
    if(is.na(dim(X)[4])){
      
    } else{
      
      globalInfectiousPressure <- theta[1]*mat%*%particles[,2,1,]/N # matrix MxK
      localInfectiousPressure <- theta[2]*particles[,2,1,]/N_M
      infProb <- 1 - exp(-(globalInfectiousPressure + localInfectiousPressure))
      probs <- rbind(infProb, matrix(1-exp(-theta[3]), nrow = M, ncol = K))
      
      
      events <- matrix(particles[, c(1,3), 1, ] - particles[,c(1,3), 2, ], nrow = 2*M, ncol = K)*
        c(rep(1, M), rep(-1, M))
      
      dbinoms <- matrix(dbinom(events, size = X[, c(1, 2), 1, ], prob = probs, log = T), nrow = 2*M, ncol = K)
      llh <- colSums(dbinoms)
      return(llh)
      
    }
  }
  return(list(sim = oneDaySim, llh = llh))
}

