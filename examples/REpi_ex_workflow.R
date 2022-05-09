rm(list = ls())

# (Skip this section if latest REpi version is installed)
# ==== Installing REpi ====

# Will need devtools to install from github
#install.packages("devtools")

# install package
#devtools::install_github("JMacDonaldPhD/REpi", ref = "main")

# ==== Using REpi (Particle Filter Example) ====

#library(REpi)
# Random Seed for reproducibility (this will generate the dataset 'sim_sample' given in the package)
set.seed(1)

# ==== Construct Epidemic Model ====
M <- 2 # No. meta-populations
N_M <- rep(1e3, M) # Population size in each meta-population
endTime <- 5 # An end time for simulation of the epidemic


# Discrete-Time (Chain-Binomial) Epidemic model which returns 
# functions for epidemic simulation and log-density calculation.
# This simulates over 30 days so most of the epidemic is observed.
# The endTime given above will refer to the amount of timepoints
# which will be observed. This will become clear later.
epiModel <- metaSIR(N_M, endTime = 30)

# Epidemic Parameters
theta <- c(0.05, 1, 0.25) # Global Infection, Local Infection, Removal respectively
I0 <- 50 # Initial number of infectives in each population
X0 <- matrix(nrow = M, ncol = 3)
X0[1, ] <- c(N_M[1] - I0, I0, 0)
for(i in 2:M){
    X0[i, ] <- c(N_M[i], 0, 0)
}

# Simulate a realisation of the epidemic 
X_sim <- epiModel$sim(param = list(X0, theta))

# Calculate the log-density of simulated Epidemic
epiModel$llh(X_sim, theta)
# -2587.295

# ==== Construct Observation Model ====

# Constructs a Case Ascertation (Binomial sample of infections) 
# Observation model which takes an underlying epidemic as its argument,
# returns sampling and log-likelihood calculation functions.
# This might change to look more like R's 'r', 'd' etc. convention.
obsModel <- caseAscObsModel(X_sim)

# = Sampling and Log-likelihood calculation =
alpha <- 0.1 # The probability that an infection is detected on any given day.

# Simulate a sample from X_sim
y <- obsModel$sample(0.1)

# Calculate the log-density of the sample y
obsModel$llh(y, alpha)

# Reduce sample to days of interest. Generating the data 
# this way ensures that the same sample is generated
# if the random seed is the same. Then the sample can
# be truncated accordingly.
y <- y[,1:endTime, drop = F]


# Reconstruct epidemic model so simulate only the
# days of interest
epiModel <- metaSIR(N_M, endTime = endTime)

# Example of how to construct a bootstrap particle filter
particleFilter <- BS_PF(y, X0, obsFrame = caseAscObsModel, epiModel = epiModel)

# Returns Log-likelihood
particleFilter(K = 10, theta, alpha)


# Looks at the distribution of the Particle Filter Estimate
PF_sample <- replicate(1000, particleFilter(K = 5, theta, alpha)$logLikeEst)
plot(density(PF_sample), main=paste0("var(log estimate) = ", var(PF_sample[!is.infinite(PF_sample)]), collapse = ""))

# Adapt_BS_PF chooses K particles such that var of the log-likelihood
# estimate is below 1. (Only adapts by increasing K, so may end up 
# with slightly too many particles)
K <- adapt_BS_PF(K0 = 10, particleFilter, theta, alpha)

# Define prior for epidemic parameters
logPrior <- function(param){
  return(sum(dunif(param, 0, 10, log = TRUE)))
}

# Adapt the MCMC proposal parameters
# Outputs plots checking whether proposal scale parameter has stabilised (similar thing can be done 
# with covariance parameters but they are not stored as of yet)
lambda0 <- 2.38/sqrt(3)
V0 <- diag(1, length(theta))
adapt_step <- adapt_particleMCMC(init = theta, epiModel = epiModel, obsFrame = caseAscObsModel,
                                 y, X0 = X0, alpha, logPrior, lambda0, V0, K = K, noIts = 3e4)



# Run MCMC with adapted proposal parameters
MCMC_sample <- particleMCMC(init = theta, epiModel = epiModel, obsFrame = caseAscObsModel,
                            y, X0 = X0, alpha, logPrior, lambda = adapt_step$lambda,
                            V = adapt_step$V, K = K, noIts = 3e4)

MCMC_sample$acceptRate

# Save a plot of epidemic output
jpeg(filename = "testParticleMCMC.jpeg", width = 480*3, height = 480*3)
plotMCMC(MCMC_sample$draws[,1:3], theta, expressions = c(expression(beta[G]), expression(beta[L]),
                                                   expression(gamma)))
dev.off()
