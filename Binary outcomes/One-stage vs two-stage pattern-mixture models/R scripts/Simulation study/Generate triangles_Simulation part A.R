##################################################################################################
#                                                                                                #
#           Code to generate 1,000 triangles with missing binary outcome data - Part A           #
#                                                                                                #
##################################################################################################



# Initially, load the following necessary package before the simulation: truncnorm. 

library(truncnorm)

# The following function will ensure that extreme values of the chosen predictive distributions are not selected:
  
random.sample <- function(mean, sd, dist){
  if(dist == "t-scaled"){
    y <- rt.scaled(10000, df = 3, mean, sd)
  } else if(dist == "lnormal"){
    y <- rlnorm(10000, mean, sd)
  }
  value <- sample(y[y >= quantile(y, prob = 0.025) & y <= quantile(y, prob = 0.975)], 1, replace = F)  # Restrict within the 95% interval
  return(value)
}

# Description of true parameters with their pre-specified scenarios.

sim    <- 1000                         # Number of simulations
mOP    <- 1.5                          # True odds ratio on 'Old intervention versus Placebo' 
mNP    <- 2                            # True odds ratio on 'New intervention versus Placebo'
size   <- c("<50", ">100")             # Size of the trials
event  <- c("low", "frequent")         # Event frequency
tau2   <- c("small", "substantial")    # Extent of between-trial variance (tau2)
pm     <- c("moderate", "large")       # Extent of missing outcome data 

# The scenarios considered in the present simulation study:

(compile <- expand.grid(list(size = size, event = event, tau2 = tau2, pm = pm)))

# The function simulation.binmod.network generates a triangle network for binary outcome. The parameter sim indicates the number of simulations.

simulation.binmod.network <- function(mOP, mNP, size, event, tau2, pm, sim){
  options(warn = -1)
  
  
  # Consistency equation to obtain TRUE odds ratio for the comparison of interest
  mNO <- exp(log(mNP) - log(mOP))                         # True odds ratio on 'New vs Old' under consistency
  
  
  # Set matrix of 'Initial' TRUE event risks in (E)xperimental- and (C)ontrol-arm  
  pE1.OP <- pC1.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo' 
  pE1.NP <- pC1.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo' 
  pE1.NO <- pC1.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old' 
  
  
  # Set matrix of 'Initial' TRUE logits in (E)xperimental- and (C)ontrol-arm 
  logitE1.OP <- logitC1.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'  
  logitE1.NP <- logitC1.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'  
  logitE1.NO <- logitC1.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'  
  
  
  # Set matrix of TRUE logits (via 'Initial' TRUE logits) in (E)xperimental- and (C)ontrol-arm 
  logitE.OP <- logitC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'   
  logitE.NP <- logitC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'  
  logitE.NO <- logitC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old' 
  
  
  # Set matrix of TRUE event risks (the back-calculated) in Experimental- and Control-arm 
  pE.OP <- pC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo' 
  pE.NP <- pC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  pE.NO <- pC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Set matrix of OBSERVED event risks in Experimental- and Control-arm 
  poE.OP <- poC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo' 
  poE.NP <- poC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  poE.NO <- poC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Set matrix of number of OBSERVED events in Experimental- and Control-arm  
  eventE.OP <- eventC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'
  eventE.NP <- eventC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  eventE.NO <- eventC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Set matrix of logIMORs (informative missingness odds ratio) in Experimental- and Control-arm 
  deltaE.OP <- deltaC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'
  deltaE.NP <- deltaC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  deltaE.NO <- deltaC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Set matrix with coding to indicate the interventions in each comparison
  treatP.OP <- matrix(1, nrow = sim, ncol = 4)    # Placebo in 'Old vs Placebo'
  treatP.NP <- matrix(1, nrow = sim, ncol = 3)    # Placebo in 'New vs Placebo'
  treatO.NO <- matrix(2, nrow = sim, ncol = 1)    # Old in 'New vs Old'
  treatO.OP <- matrix(2, nrow = sim, ncol = 4)    # Old in 'Old vs Placebo'
  treatN.NP <- matrix(3, nrow = sim, ncol = 3)    # New in 'New vs Placebo'
  treatN.NO <- matrix(3, nrow = sim, ncol = 1)    # New in 'New vs Old'
  
  
  # Set matrix of number of arms in each trial (all two-arms)
  id.OP <- matrix(NA, nrow = sim, ncol = 4)    # 'Old vs Placebo'
  id.NP <- matrix(NA, nrow = sim, ncol = 3)    # 'New vs Placebo'
  id.NO <- matrix(NA, nrow = sim, ncol = 1)    # 'New vs Old'
  
  
  # Set matrix of number of arms in each trial (all two-arms)
  na.OP <- matrix(2, nrow = sim, ncol = 4)    # 'Old vs Placebo'
  na.NP <- matrix(2, nrow = sim, ncol = 3)    # 'New vs Placebo'
  na.NO <- matrix(2, nrow = sim, ncol = 1)    # 'New vs Old'
  
  
  # Set matrix of sample size in Experimental- and Control-arm
  nE.OP <- nC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'
  nE.NP <- nC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  nE.NO <- nC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Set matrix of number of missing outcome data (MOD) in Experimental- and Control-arm
  mE.OP <- mC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'
  mE.NP <- mC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  mE.NO <- mC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Set matrix of risk of MOD in Experimental- and Control-arm
  pmE.OP <- pmC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'
  pmE.NP <- pmC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  pmE.NO <- pmC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Set matrix of necessary quantities in the linking equation via IMOR in Experimental- and Control-arm 
  # (PMID: 25809313, Equations A.3, page 2075)
  AE.OP <- AC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'
  AE.NP <- AC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  AE.NO <- AC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  BE.OP <- BC.OP <- matrix(NA, nrow = sim, ncol = 4)    # Old/Placebo in 'Old vs Placebo'
  BE.NP <- BC.NP <- matrix(NA, nrow = sim, ncol = 3)    # New/Placebo in 'New vs Placebo'
  BE.NO <- BC.NO <- matrix(NA, nrow = sim, ncol = 1)    # New/Old in 'New vs Old'
  
  
  # Arrays of separate datasets: one for each direct comparison of the network
  dataset.OP <- array(NA, dim = c(4, 10, sim))    # 'Old vs Placebo'  
  dataset.NP <- array(NA, dim = c(3, 10, sim))    # 'New vs Placebo'
  dataset.NO <- array(NA, dim = c(1, 10, sim))    # 'NeW vs Old'
  
  
  # List with ALL datasets to be used for NMA via GeMTC
  arm.level.dataset <- array(0, dim = c(4 + 3 + 1, 10, sim), dimnames = list(NULL, c("id","r1","r2","m1","m2",
                                                                                     "n1","n2","t1","t2","na"), NULL))
   
  for(i in 1:sim){
    
    if(event == "low"){
      # Generate 'Initial' TRUE event risks in Control-am 
      pC1.OP[i, ] <- runif(4, 0.05, 0.09)    # Placebo in 'Old vs Placebo'
      pC1.NP[i, ] <- runif(3, 0.05, 0.09)    # Placebo in 'New vs Placebo'
      pC1.NO[i, ] <- runif(1, 0.10, 0.15)    # Old in 'New vs Old'
      # Obtain 'Initial' TRUE event risks in Experimental-arm 
      pE1.OP[i, ] <- (pC1.OP[i, ]*mOP)/(1 - pC1.OP[i, ] + pC1.OP[i, ]*mOP)    # Old in 'Old vs Placebo'
      pE1.NP[i, ] <- (pC1.NP[i, ]*mNP)/(1 - pC1.NP[i, ] + pC1.NP[i, ]*mNP)    # New in 'New vs Placebo'
      pE1.NO[i, ] <- (pC1.NO[i, ]*mNO)/(1 - pC1.NO[i, ] + pC1.NO[i, ]*mNO)    # New in 'New vs Old'
      
    } else if(event == "frequent"){
      # Generate 'Initial' TRUE event risks in Control-am 
      pC1.OP[i, ] <- runif(4, 0.27, 0.40)    # Placebo in 'Old vs Placebo'
      pC1.NP[i, ] <- runif(3, 0.27, 0.40)    # Placebo in 'New vs Placebo'
      pC1.NO[i, ] <- runif(1, 0.63, 0.76)    # Old in 'New vs Old'
      # Obtain 'Initial' TRUE event risks in Experimental-arm 
      pE1.OP[i, ] <- (pC1.OP[i, ]*mOP)/(1 - pC1.OP[i, ] + pC1.OP[i, ]*mOP)    # Old in 'Old vs Placebo'
      pE1.NP[i, ] <- (pC1.NP[i, ]*mNP)/(1 - pC1.NP[i, ] + pC1.NP[i, ]*mNP)    # New in 'New vs Placebo'
      pE1.NO[i, ] <- (pC1.NO[i, ]*mNO)/(1 - pC1.NO[i, ] + pC1.NO[i, ]*mNO)    # New in 'New vs Old'
    }

    
    # Obtain 'Initial' TRUE logitS in Experimental-arm 
    logitE1.OP[i, ] <- log(pE1.OP[i, ]/(1 - pE1.OP[i, ]))    # Old in 'Old vs Placebo'
    logitE1.NP[i, ] <- log(pE1.NP[i, ]/(1 - pE1.NP[i, ]))    # New in 'New vs Placebo'
    logitE1.NO[i, ] <- log(pE1.NO[i, ]/(1 - pE1.NO[i, ]))    # New in 'New vs Old'
    # Obtain 'Initial' TRUE logitS in Control-arm 
    logitC1.OP[i, ] <- log(pC1.OP[i, ]/(1 - pC1.OP[i, ]))    # Placebo in 'Old vs Placebo'
    logitC1.NP[i, ] <- log(pC1.NP[i, ]/(1 - pC1.NP[i, ]))    # Placebo in 'New vs Placebo'
    logitC1.NO[i, ] <- log(pC1.NO[i, ]/(1 - pC1.NO[i, ]))    # Old in 'New vs Old'
    
    
    ## Use predictive prior for 'all-cause mortality' for 'pharma vs. PBO/control' 
    ## (Table IV in PMID: 25475839) - Small between-trial variance
    if(tau2 == "small"){ 
      # Generate TRUE logits (via 'Initial' TRUE logitS) in Experimental-arm 
      logitE.OP[i, ] <- rnorm(4, logitE1.OP[i, ], sqrt(2*random.sample(-3.95, 1.34, "lnormal")/3))  # Old in 'Old vs Placebo'
      logitE.NP[i, ] <- rnorm(3, logitE1.NP[i, ], sqrt(2*random.sample(-3.95, 1.34, "lnormal")/3))  # New in 'New vs Placebo'
      logitE.NO[i, ] <- rnorm(1, logitE1.NO[i, ], sqrt(random.sample(-3.95, 1.34, "lnormal")/2))    # New in 'New vs Old'
      # Generate TRUE logits (via 'Initial' TRUE logitS) in Generate Control-arm 
      logitC.OP[i, ] <- rnorm(4, logitC1.OP[i, ], sqrt(random.sample(-3.95, 1.34, "lnormal")/3))    # Placebo in 'Old vs Placebo'
      logitC.NP[i, ] <- rnorm(3, logitC1.NP[i, ], sqrt(random.sample(-3.95, 1.34, "lnormal")/3))    # Placebo in 'New vs Placebo'
      logitC.NO[i, ] <- rnorm(1, logitC1.NO[i, ], sqrt(random.sample(-3.95, 1.34, "lnormal")/2))    # Old in 'New vs Old'
      
      ## Use predictive prior for 'generic healthcare setting' for 'pharma vs. PBO/control' 
      ## (Table IV in PMID: 25475839) - Substantial between-trial variance
    } else if(tau2 == "substantial"){ 
      # Generate TRUE logits (via 'Initial' TRUE logitS) in Experimental-arm 
      logitE.OP[i, ] <- rnorm(4, logitE1.OP[i, ], sqrt(2*random.sample(-2.56, 1.74, "lnormal")/3))   # Old in 'Old vs Placebo'
      logitE.NP[i, ] <- rnorm(3, logitE1.NP[i, ], sqrt(2*random.sample(-2.56, 1.74, "lnormal")/3))   # New in 'New vs Placebo'
      logitE.NO[i, ] <- rnorm(1, logitE1.NO[i, ], sqrt(random.sample(-2.56, 1.74, "lnormal")/2))     # New in 'New vs Old'
      # Generate TRUE logits (via 'Initial' TRUE logitS) in Generate Control-arm 
      logitC.OP[i, ] <- rnorm(4, logitC1.OP[i, ], sqrt(random.sample(-2.56, 1.74, "lnormal")/3))    # Placebo in 'Old vs Placebo'
      logitC.NP[i, ] <- rnorm(3, logitC1.NP[i, ], sqrt(random.sample(-2.56, 1.74, "lnormal")/3))    # Placebo in 'New vs Placebo'
      logitC.NO[i, ] <- rnorm(1, logitC1.NO[i, ], sqrt(random.sample(-2.56, 1.74, "lnormal")/2))    # Old in 'New vs Old'
    }
    
    
    # Obtain 'back-calculated' TRUE event risks in Experimental-arm 
    pE.OP[i, ] <- 1/(1 + exp(-logitE.OP[i, ]))    # Old in 'Old vs Placebo'
    pE.NP[i, ] <- 1/(1 + exp(-logitE.NP[i, ]))    # New in 'New vs Placebo'
    pE.NO[i, ] <- 1/(1 + exp(-logitE.NO[i, ]))    # New in 'New vs Old'
    # Obtain 'back-calculated' TRUE event risks in Control-arm
    pC.OP[i, ] <- 1/(1 + exp(-logitC.OP[i, ]))    # Placebo in 'Old vs Placebo'
    pC.NP[i, ] <- 1/(1 + exp(-logitC.NP[i, ]))    # Placebo in 'New vs Placebo'
    pC.NO[i, ] <- 1/(1 + exp(-logitC.NO[i, ]))    # Old in 'New vs Old'
    
    
    # Set matrix of number of arms in each trial (all two-arms)
    id.OP[i, ] <- c(1:(4))                                # 'Old vs Placebo'
    id.NP[i, ] <- c(((4) + 1):((4) + (3)))                # 'New vs Placebo'
    id.NO[i, ] <- c(((4) + (3) + 1):((4) + (3) + (1)))    # 'New vs Old'
    
    
    # Generate Sample size in Experimental- and Control-arm (equal in both arms)
    if(size == "<50"){
      nE.OP[i, ] <- nC.OP[i, ] <- round(runif(4, 12, 39), 0)    # Old/Placebo in 'Old vs Placebo'
      nE.NP[i, ] <- nC.NP[i, ] <- round(runif(3, 12, 39), 0)    # New/Placebo in 'New vs Placebo'
      nE.NO[i, ] <- nC.NO[i, ] <- round(runif(1, 15, 49), 0)    # New/Old in 'New vs Old'
    } else if(size == ">100"){
      nE.OP[i, ] <- nC.OP[i, ] <- round(runif(4, 102, 187), 0)  # Old/Placebo in 'Old vs Placebo'
      nE.NP[i, ] <- nC.NP[i, ] <- round(runif(3, 102, 187), 0)  # New/Placebo in 'New vs Placebo'
      nE.NO[i, ] <- nC.NO[i, ] <- round(runif(1, 128, 241), 0)  # New/Old in 'New vs Old'
    }

    
    ## Generate risk of MOD, while accounting for extent (low, moderate, large) and balance (balance or imbalance) of MOD 
    if(pm == "moderate"){ # Moderate & unbalanced MOD
      # Generate risk of MOD in Experimental-arm 
      pmE.OP[i, ] <- runif(4, 0.05, 0.10)    # Old in 'Old vs Placebo'
      pmE.NP[i, ] <- runif(3, 0.05, 0.10)    # New in 'New vs Placebo'
      pmE.NO[i, ] <- runif(1, 0.05, 0.10)    # New in 'New vs Old'
      # Generate risk of MOD in Control-arm 
      pmC.OP[i, ] <- runif(4, 0.11, 0.20)    # Placebo in 'Old vs Placebo'
      pmC.NP[i, ] <- runif(3, 0.11, 0.20)    # Placebo in 'New vs Placebo'
      pmC.NO[i, ] <- runif(1, 0.11, 0.20)    # Old in 'New vs Old'
      
    } else {              # Large & unbalanced MOD
      # Generate risk of MOD in Experimental-arm
      pmE.OP[i, ] <- runif(4, 0.21, 0.30)    # Old in 'Old vs Placebo'
      pmE.NP[i, ] <- runif(3, 0.21, 0.30)    # New in 'New vs Placebo'
      pmE.NO[i, ] <- runif(1, 0.21, 0.30)    # New in 'New vs Old'
      # Generate risk of MOD in Control-arm 
      pmC.OP[i, ] <- runif(4, 0.31, 0.40)    # Placebo in 'Old vs Placebo'
      pmC.NP[i, ] <- runif(3, 0.31, 0.40)    # Placebo in 'New vs Placebo'
      pmC.NO[i, ] <- runif(1, 0.31, 0.40)    # Old in 'New vs Old'
    }
    
    
    # Generate number of MOD in Experimental-arm 
    mE.OP[i, ] <- rbinom(4, nE.OP[i, ], pmE.OP[i, ])    # Old in 'Old vs Placebo'
    mE.NP[i, ] <- rbinom(3, nE.NP[i, ], pmE.NP[i, ])    # New in 'New vs Placebo'
    mE.NO[i, ] <- rbinom(1, nE.NO[i, ], pmE.NO[i, ])    # New in 'New vs Old'
    # Generate number of MOD in Control-arm 
    mC.OP[i, ] <- rbinom(4, nC.OP[i, ], pmC.OP[i, ])    # Placebo in 'Old vs Placebo'
    mC.NP[i, ] <- rbinom(3, nC.NP[i, ], pmC.NP[i, ])    # Placebo in 'New vs Placebo'
    mC.NO[i, ] <- rbinom(1, nC.NO[i, ], pmC.NO[i, ])    # Old in 'New vs Old'
    
    
    # Generate logIMOR in Experimental-arm (MNAR)
    deltaE.OP[i, ] <- rtruncnorm(4, a = log(1), b = Inf, log(2), 1)    # Old in 'Old vs Placebo'
    deltaE.NP[i, ] <- rtruncnorm(3, a = log(1), b = Inf, log(2), 1)    # New in 'New vs Placebo'
    deltaE.NO[i, ] <- rtruncnorm(1, a = log(1), b = Inf, log(2), 1)    # New in 'New vs Old'

    # Generate logIMOR in Control-arm (MNAR)
    deltaC.OP[i, ] <- rtruncnorm(4, a = -Inf, b = log(1), -log(2), 1)   # Placebo in 'Old vs Placebo'
    deltaC.NP[i, ] <- rtruncnorm(3, a = -Inf, b = log(1), -log(2), 1)   # Placebo in 'New vs Placebo'
    deltaC.NO[i, ] <- rtruncnorm(1, a = log(1), b = Inf, log(2), 1)     # Old in 'New vs Old'

    
    # Obtain A component of linkage function (Experimental-arm) (PMID: 25809313, Equations A.3, page 2075)
    AE.OP[i, ] <- (pmE.OP[i, ] - pE.OP[i, ])*(1 - exp(deltaE.OP[i, ])) - 1    # Old in 'Old vs Placebo'
    AE.NP[i, ] <- (pmE.NP[i, ] - pE.NP[i, ])*(1 - exp(deltaE.NP[i, ])) - 1    # New in 'New vs Placebo'
    AE.NO[i, ] <- (pmE.NO[i, ] - pE.NO[i, ])*(1 - exp(deltaE.NO[i, ])) - 1    # New in 'New vs Old'
    # Obtain A component of linkage function (Control-arm)
    AC.OP[i, ] <- (pmC.OP[i, ] - pC.OP[i, ])*(1 - exp(deltaC.OP[i, ])) - 1    # Placebo in 'Old vs Placebo'
    AC.NP[i, ] <- (pmC.NP[i, ] - pC.NP[i, ])*(1 - exp(deltaC.NP[i, ])) - 1    # Placebo in 'New vs Placebo'
    AC.NO[i, ] <- (pmC.NO[i, ] - pC.NO[i, ])*(1 - exp(deltaC.NO[i, ])) - 1    # Old in 'New vs Old'
    
    
    # Obtain B component of linkage function (Experimental-arm) (PMID: 25809313, Equations A.3, page 2075)
    BE.OP[i, ] <- 2*(1 - exp(deltaE.OP[i, ]))*(1 - pmE.OP[i, ])    # Old in 'Old vs Placebo'
    BE.NP[i, ] <- 2*(1 - exp(deltaE.NP[i, ]))*(1 - pmE.NP[i, ])    # New in 'New vs Placebo'
    BE.NO[i, ] <- 2*(1 - exp(deltaE.NO[i, ]))*(1 - pmE.NO[i, ])    # New in 'New vs Old'
    # Obtain B component of linkage function (Control-arm) 
    BC.OP[i, ] <- 2*(1 - exp(deltaC.OP[i, ]))*(1 - pmC.OP[i, ])    # Placebo in 'Old vs Placebo'
    BC.NP[i, ] <- 2*(1 - exp(deltaC.NP[i, ]))*(1 - pmC.NP[i, ])    # Placebo in 'New vs Placebo'
    BC.NO[i, ] <- 2*(1 - exp(deltaC.NO[i, ]))*(1 - pmC.NO[i, ])    # Old in 'New vs Old'
    
    
    ## Obtain risk of observed events in Experimental-arm (PMID: 25809313, Equation A.4, page 2075)
    # Old in 'Old vs Placebo'
    poE.OP[i, ] <- max(0, min(1, (-AE.OP[i, ] - sqrt(AE.OP[i, ]*AE.OP[i, ] - 2*pE.OP[i, ]*BE.OP[i, ]))/BE.OP[i, ] ))  
    # New in 'New vs Placebo'  
    poE.NP[i, ] <- max(0, min(1, (-AE.NP[i, ] - sqrt(AE.NP[i, ]*AE.NP[i, ] - 2*pE.NP[i, ]*BE.NP[i, ]))/BE.NP[i, ] ))  
    # New in 'New vs Old'
    poE.NO[i, ] <- max(0, min(1, (-AE.NO[i, ] - sqrt(AE.NO[i, ]*AE.NO[i, ] - 2*pE.NO[i, ]*BE.NO[i, ]))/BE.NO[i, ] ))        

    ## Obtain risk of observed events in Control-arm 
    # Placebo in 'Old vs Placebo'
    poC.OP[i, ] <- max(0, min(1, (-AC.OP[i, ] - sqrt(AC.OP[i, ]*AC.OP[i, ] - 2*pC.OP[i, ]*BC.OP[i, ]))/BC.OP[i, ] ))  
    # Placebo in 'New vs Placebo'  
    poC.NP[i, ] <- max(0, min(1, (-AC.NP[i, ] - sqrt(AC.NP[i, ]*AC.NP[i, ] - 2*pC.NP[i, ]*BC.NP[i, ]))/BC.NP[i, ] ))
    # Old in 'New vs Old'    
    poC.NO[i, ] <- max(0, min(1, (-AC.NO[i, ] - sqrt(AC.NO[i, ]*AC.NO[i, ] - 2*pC.NO[i, ]*BC.NO[i, ]))/BC.NO[i, ] ))    
    
    
    # Generate number of observed events in Experimental-arm 
    eventE.OP[i, ] <- rbinom(4, nE.OP[i, ] - mE.OP[i, ], poE.OP[i, ])    # Old in 'Old vs Placebo'
    eventE.NP[i, ] <- rbinom(3, nE.NP[i, ] - mE.NP[i, ], poE.NP[i, ])    # New in 'New vs Placebo'
    eventE.NO[i, ] <- rbinom(1, nE.NO[i, ] - mE.NO[i, ], poE.NO[i, ])    # New in 'New vs Old'
    # Generate number of observed events in Control-arm 
    eventC.OP[i, ] <- rbinom(4, nC.OP[i, ] - mC.OP[i, ], poC.OP[i, ])    # Placebo in 'Old vs Placebo'
    eventC.NP[i, ] <- rbinom(3, nC.NP[i, ] - mC.NP[i, ], poC.NP[i, ])    # Placebo in 'New vs Placebo'
    eventC.NO[i, ] <- rbinom(1, nC.NO[i, ] - mC.NO[i, ], poC.NO[i, ])    # Old in 'New vs Old'
    
    
    # List of separate datasets: one for each direct comparison of the network
    dataset.OP[,, i] <- cbind(id.OP[i, ], eventE.OP[i, ], eventC.OP[i, ], mE.OP[i, ], mC.OP[i, ], nE.OP[i, ], nC.OP[i, ], 
                              treatO.OP[i, ], treatP.OP[i, ], na.OP[i, ])    # 'Old vs Placebo'
    dataset.NP[,, i] <- cbind(id.NP[i, ], eventE.NP[i, ], eventC.NP[i, ], mE.NP[i, ], mC.NP[i, ], nE.NP[i, ], nC.NP[i, ], 
                              treatN.NP[i, ], treatP.NP[i, ], na.NP[i, ])    # 'New vs Placebo'
    dataset.NO[,, i] <- cbind(id.NO[i, ], eventE.NO[i, ], eventC.NO[i, ], mE.NO[i, ], mC.NO[i, ], nE.NO[i, ], nC.NO[i, ], 
                              treatN.NO[i, ], treatO.NO[i, ], na.NO[i, ])    # 'New vs Old'
    
    
    # Array with ALL datasets 
    arm.level.dataset[,, i] <- rbind(dataset.OP[,, i], dataset.NP[,, i], dataset.NO[,, i])
  }
  return(arm.level.dataset)
} 

# Collect the simulated triangles for all scenarios and necessary number of simulations:

mat <- list()
set.seed(123)
for(l in 1:length(compile[, 1])){  
# The index 'l' corresponds to a specific scenario as a list and it contains all simulations as indicated by 'sim'.
  mat[[l]] <- simulation.binmod.network(mOP, mNP, compile[l, 'size'], compile[l, 'event'], compile[l, 'tau2'], 
                                        compile[l, 'pm'], sim)
}

# Save the generated triangles as .RData (here, it is saved as '1000 simulated two-arm networks'):

save(mat, file = "./One-stage vs two-stage pattern-mixture models/R scripts/Simulation study/1000 simulated two-arm networks.RData")

# Then, the generated triangles can be analyzed using the one-stage and two-stage pattern-mixture models (see, Part B)

