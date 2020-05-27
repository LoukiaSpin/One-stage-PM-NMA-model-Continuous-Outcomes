#####################################################################################################################
#                                                                                                                   #
#               A function to perform Bayesian network meta-analysis for aggregate continuous outcomes              # 
#                                 <Normal likelihood, identity link, Random Effects>                                # 
#                                  <Normal likelihood, identity link, Fixed Effect>                                 #    
#                                 (Dias et al., 2013 in Appendix  - PMID: 23104435)                                 #    
#           One-stage pattern-mixture model with relative missingness parameters under Missing at Random            #
#               Informative Missingness Difference of Means & Informative Missingness Ratio of Means                #
#                                     (Mavridis et al., 2015 - PMID: 25393541)                                      #
#                      <Various structures of prior distribution on the missingness parameter>                      #
#                                                                                                                   #
#####################################################################################################################



## Necessary arguments (with explanations)
# data: an arm-level, wide-format dataset where each line is a trial. The file should be a delimited text-file (.txt). 
#       tX, intervention identifier in arm X; y: the observed mean value in arm X; seX, the observed standard error in arm X; 
#       mX, number of missing outcome data in arm X; cX: the number of completers in arm X; nX, the number randomised in arm X;
#       NMA, identifier for each network (necessary even when a single network is analysed);
# meta.model: a character to indicate whether random-effects or fixed-effect network meta-analysis will be applied.
#             Possible values are "RE", and "FE" to indicate a random-effects and a fixed-effect network meta-analysis, respectively.
# measure: a character to indicate the effect measure for continuous outcomes. Possible values are "MD", "SMD", and "ROM"to indicate the mean difference,
#          standardised mean difference, and ratio of means (analysed in the logarithmic scale), respectively.
# model: a character to indicate the prior structure for the missingness parameter.
#        Possible values are "identical", "hierarchical", "independent". 
# var.misspar: a numerical value for the prior variance of the informative missingness difference of means (IMDOM) and the informative missingness ratio of means (IMROM; in the logarithim value). 
#              The corresponding Plausible values for IMDOM and logIMROM include 1 and 0.01 for liberal belief, 9 and 0.04 for conservative belief about the missing at random assumption.
# dir: a binary value that indicates the direction of the  outcome (1: positive outcome; 0: negative outcome) 
# n.chains: number of Markov chains (default: 3 - https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).
# n.iter: number of total iterations per chain (including burn in; default: 2000 - https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).
# n.burnin: length of burn in, i.e. number of iterations to discard at the beginning (https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).
# n.thin: thinning rate. Must be a positive integer (https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).



NMA.IMDOM.IMROM.PM.model <- function(data, meta.model, measure, model, var.misspar, dir, n.chains, n.iter, n.burnin, n.thin){

  
  ## Necessary function to return the posterior distribution of core NMA parameters from each model
  source("./Pattern-mixture IMDOM & IMROM models/R scripts/Functions/Full NMA model/collect jags results Full NMA PM Continuous_function.R")
  
  library(dplyr)    # For data-management
  library(R2jags)   # To implement Bayesian analysis in JAGS 
  
  
  ## Arm-level, wide-format dataset  
  # START
  (y0 <- data %>% dplyr::select(starts_with("y")))            # Observed mean value in each arm of every trial
  (se0 <- data %>% dplyr::select(starts_with("se")))          # Observed standard error in each arm of every trial
  (m <- data %>% dplyr::select(starts_with("m")))             # Number of missing participants in each arm of every trial
  (c <- data %>% dplyr::select(starts_with("c")))             # Number completers in each arm of every trial
  (sd0 <- round(se0*sqrt(c), 2))                              # Observed standard deviation in each arm of every trial
  (N <- m + c)                                                # Number of randomised participants in each arm of every trial
  (t <- data %>% dplyr::select(starts_with("t")))             # Intervention studied in each arm of every trial 
  (J <- data[, c("NMA")])                                     # Identify unique networks (e.g., 1, 2, etc)
  (L <- length(unique(J)))                                    # Number of unique networks in the dataset
  narms <- apply(t, 1, function(x) length(which(!is.na(x))))
  ref <- ns <- nt <- rep(0, L) 
  na <- sigma <- list()
  for(i in 1:L){
    nt[i] <- length(table(as.matrix(t[J == i, ])))            # Total number of interventions per network
    ns[i] <- length(y0[J == i, 1])                            # Total number of included trials per network
    ref[i] <- which.max(table(as.matrix(t[J == i, ])))        # Reference intervention per network: the most frequently appeared intervention in the network
    na[[i]] <- narms[J == i]                                  # Number of interventions investigated in every trial per network
    # Trial-specific observed pooled standard deviation
    (sigma[[i]] <- sqrt(apply((sd0[J == i, ]^2)*(c[J == i, ] - 1), 1, sum, na.rm = T)/(apply(c[J == i, ], 1, sum, na.rm = T) - unlist(na[[i]])))) 
  }
  #END


  ## Information for the prior distribution on the missingness parameter (IMDOM or logIMROM) 
  # START
  M <- list()
  for(i in 1:L){
    M[[i]] <- ifelse(!is.na(y0[J == i, ]), 0, NA)             # Vector of as many 0s (missing at random assumption) as the number of arms in trial i (independent structure)
  }
  prec.misspar <- 1/var.misspar
  psi.misspar <- sqrt(var.misspar)                            # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
  cov.misspar <- 0.5*var.misspar                              # covariance of pair of missingness parameters in a trial (independent structure)
  # END

  
  if(meta.model == "RE"){
    
    ## Define the pattern-mixture models for various prior structures of the missingness parameter under the random-effects assumption
    # START
    if(model == "identical"){
      
      data.jag.ide <- jagsfit.ide.common <- jagsfit.ide.trial <- jagsfit.ide.arm <- list()
      
      if(measure == "MD"){
        
        param.jags.ide <- c("theta", "MD", "tau2", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ide[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ide.common[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture IDE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.trial[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture IDE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.arm[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture IDE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ide.common, meta.model = "RE", model = "identical", assumption = "common", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.trial, meta.model = "RE", model = "identical", assumption = "trial", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.arm, meta.model = "RE", model = "identical", assumption = "arm", measure = "MD", nt, ns, na, L))
        
      } else if(measure == "SMD"){
        
        param.jags.ide <- c("theta", "SMD", "tau2", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ide[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ide.common[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture IDE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.trial[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture IDE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.arm[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture IDE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ide.common, meta.model = "RE", model = "identical", assumption = "common", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.trial, meta.model = "RE", model = "identical", assumption = "trial", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.arm, meta.model = "RE", model = "identical", assumption = "arm", measure = "SMD", nt, ns, na, L))
        
      } else {
        
        param.jags.ide <- c("theta", "LRoM", "tau2", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ide[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ide.common[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture IDE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.trial[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture IDE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.arm[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture IDE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ide.common, meta.model = "RE", model = "identical", assumption = "common", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.trial, meta.model = "RE", model = "identical", assumption = "trial", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.arm, meta.model = "RE", model = "identical", assumption = "arm", measure = "ROM", nt, ns, na, L))
        
      }
      
    } else if(model == "hierarchical"){
      
      data.jag.hie <- jagsfit.hie.common <- jagsfit.hie.trial <- jagsfit.hie.arm <- list()
      
      if(measure == "MD"){
        
        param.jags.hie <- c("theta", "MD", "tau2", "SUCRA", "order", "mean.phi", "sd.phi")
        
        for(i in 1:L){
          
          data.jag.hie[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "precd.phi" = prec.misspar, "psi.phi" = psi.misspar, "D" = dir[i])
          
          jagsfit.hie.common[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture HIE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.trial[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture HIE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.arm[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture HIE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.hie.common, meta.model = "RE", model = "hierarchical", assumption = "common", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.trial, meta.model = "RE", model = "hierarchical", assumption = "trial", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.arm, meta.model = "RE", model = "hierarchical", assumption = "arm", measure = "MD", nt, ns, na, L))
        
      } else if(measure == "SMD"){
        
        param.jags.hie <- c("theta", "SMD", "tau2", "SUCRA", "order", "mean.phi", "sd.phi")
        
        for(i in 1:L){
          
          data.jag.hie[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "precd.phi" = prec.misspar, "psi.phi" = psi.misspar, "D" = dir[i])
          
          jagsfit.hie.common[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture HIE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.trial[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture HIE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.arm[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture HIE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.hie.common, meta.model = "RE", model = "hierarchical", assumption = "common", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.trial, meta.model = "RE", model = "hierarchical", assumption = "trial", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.arm, meta.model = "RE", model = "hierarchical", assumption = "arm", measure = "SMD", nt, ns, na, L))
        
      } else {
        
        param.jags.hie <- c("theta", "LRoM", "tau2", "SUCRA", "order", "mean.phi", "sd.phi")
        
        for(i in 1:L){
          
          data.jag.hie[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "precd.phi" = prec.misspar, "psi.phi" = psi.misspar, "D" = dir[i])
          
          jagsfit.hie.common[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture HIE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.trial[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture HIE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.arm[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture HIE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.hie.common, meta.model = "RE", model = "hierarchical", assumption = "common", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.trial, meta.model = "RE", model = "hierarchical", assumption = "trial", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.arm, meta.model = "RE", model = "hierarchical", assumption = "arm", measure = "ROM", nt, ns, na, L))
      } 
      
    } else {
      
      data.jag.ind.corr <- data.jag.ind.unc <- jagsfit.ind.corr <- jagsfit.ind.unc <- list()
      
      if(measure == "MD"){
        
        param.jags.ind <- c("theta", "MD", "tau2", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ind.corr[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "M" = M[[i]], "var.phi" = var.misspar, "cov.phi" = cov.misspar, "D" = dir[i])
          
          data.jag.ind.unc[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ind.corr[[i]] <- jags(data = data.jag.ind.corr[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture IND-Corr.txt",
                                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ind.unc[[i]] <- jags(data = data.jag.ind.unc[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Mean Difference/Full RE-NMA MD-IMDOM Pattern-mixture IND-Uncorr.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ind.corr, meta.model = "RE", model = "independent", assumption = "correlated", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ind.unc, meta.model = "RE", model = "independent", assumption = "uncorrelated", measure = "MD", nt, ns, na, L))
        
      } else if(measure == "SMD"){
        
        param.jags.ind <- c("theta", "SMD", "tau2", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ind.corr[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "M" = M[[i]], "var.phi" = var.misspar, "cov.phi" = cov.misspar, "D" = dir[i])
          
          data.jag.ind.unc[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ind.corr[[i]] <- jags(data = data.jag.ind.corr[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture IND-Corr.txt",
                                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ind.unc[[i]] <- jags(data = data.jag.ind.unc[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Standardised Mean Difference/Full RE-NMA SMD-IMDOM Pattern-mixture IND-Uncorr.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ind.corr, meta.model = "RE", model = "independent", assumption = "correlated", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ind.unc, meta.model = "RE", model = "independent", assumption = "uncorrelated", measure = "SMD", nt, ns, na, L))
        
      } else {
        
        param.jags.ind <- c("theta", "LRoM", "tau2", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ind.corr[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "M" = M[[i]], "var.phi" = var.misspar, "cov.phi" = cov.misspar, "D" = dir[i])
          
          data.jag.ind.unc[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ind.corr[[i]] <- jags(data = data.jag.ind.corr[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture IND-Corr.txt",
                                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ind.unc[[i]] <- jags(data = data.jag.ind.unc[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Random-effects model/Ratio of Means/Full RE-NMA ROM-IMROM Pattern-mixture IND-Uncorr.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ind.corr, meta.model = "RE", model = "independent", assumption = "correlated", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ind.unc, meta.model = "RE", model = "independent", assumption = "uncorrelated", measure = "ROM", nt, ns, na, L))
        
      }
      
    } 
    
  } else {  
    
    ## Define the pattern-mixture models for various prior structures of the missingness parameter under the fixed-effect assumption
    # START
    if(model == "identical"){
      
      data.jag.ide <- jagsfit.ide.common <- jagsfit.ide.trial <- jagsfit.ide.arm <- list()
      
      if(measure == "MD"){
        
        param.jags.ide <- c("MD", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ide[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ide.common[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture IDE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.trial[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture IDE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.arm[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture IDE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ide.common, meta.model = "FE", model = "identical", assumption = "common", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.trial, meta.model = "FE", model = "identical", assumption = "trial", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.arm, meta.model = "FE", model = "identical", assumption = "arm", measure = "MD", nt, ns, na, L))
        
      } else if(measure == "SMD"){
        
        param.jags.ide <- c("SMD", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ide[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ide.common[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture IDE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.trial[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture IDE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.arm[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture IDE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ide.common, meta.model = "FE", model = "identical", assumption = "common", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.trial, meta.model = "FE", model = "identical", assumption = "trial", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.arm, meta.model = "FE", model = "identical", assumption = "arm", measure = "SMD", nt, ns, na, L))
        
      } else {
        
        param.jags.ide <- c("LRoM", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ide[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ide.common[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture IDE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.trial[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture IDE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ide.arm[[i]] <- jags(data = data.jag.ide[[i]], parameters.to.save = param.jags.ide, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture IDE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ide.common, meta.model = "FE", model = "identical", assumption = "common", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.trial, meta.model = "FE", model = "identical", assumption = "trial", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ide.arm, meta.model = "FE", model = "identical", assumption = "arm", measure = "ROM", nt, ns, na, L))
        
      }
      
    } else if(model == "hierarchical"){
      
      data.jag.hie <- jagsfit.hie.common <- jagsfit.hie.trial <- jagsfit.hie.arm <- list()
      
      if(measure == "MD"){
        
        param.jags.hie <- c("MD", "SUCRA", "order", "mean.phi", "sd.phi")
        
        for(i in 1:L){
          
          data.jag.hie[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "precd.phi" = prec.misspar, "psi.phi" = psi.misspar, "D" = dir[i])
          
          jagsfit.hie.common[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture HIE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.trial[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture HIE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.arm[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture HIE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.hie.common, meta.model = "FE", model = "hierarchical", assumption = "common", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.trial, meta.model = "FE", model = "hierarchical", assumption = "trial", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.arm, meta.model = "FE", model = "hierarchical", assumption = "arm", measure = "MD", nt, ns, na, L))
        
      } else if(measure == "SMD"){
        
        param.jags.hie <- c("SMD", "SUCRA", "order", "mean.phi", "sd.phi")
        
        for(i in 1:L){
          
          data.jag.hie[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "precd.phi" = prec.misspar, "psi.phi" = psi.misspar, "D" = dir[i])
          
          jagsfit.hie.common[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture HIE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.trial[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture HIE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.arm[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture HIE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.hie.common, meta.model = "FE", model = "hierarchical", assumption = "common", measure = "sMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.trial, meta.model = "FE", model = "hierarchical", assumption = "trial", measure = "sMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.arm, meta.model = "FE", model = "hierarchical", assumption = "arm", measure = "sMD", nt, ns, na, L))
        
      } else {
        
        param.jags.hie <- c("LRoM", "SUCRA", "order", "mean.phi", "sd.phi")
        
        for(i in 1:L){
          
          data.jag.hie[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "precd.phi" = prec.misspar, "psi.phi" = psi.misspar, "D" = dir[i])
          
          jagsfit.hie.common[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture HIE-Common.txt",
                                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.trial[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture HIE-Trial.txt",
                                         n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.hie.arm[[i]] <- jags(data = data.jag.hie[[i]], parameters.to.save = param.jags.hie, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture HIE-Arm.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.hie.common, meta.model = "FE", model = "hierarchical", assumption = "common", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.trial, meta.model = "FE", model = "hierarchical", assumption = "trial", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.hie.arm, meta.model = "FE", model = "hierarchical", assumption = "arm", measure = "ROM", nt, ns, na, L))
      } 
      
    } else {
      
      data.jag.ind.corr <- data.jag.ind.unc <- jagsfit.ind.corr <- jagsfit.ind.unc <- list()
      
      if(measure == "MD"){
        
        param.jags.ind <- c("MD", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ind.corr[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "M" = M[[i]], "var.phi" = var.misspar, "cov.phi" = cov.misspar, "D" = dir[i])
          
          data.jag.ind.unc[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ind.corr[[i]] <- jags(data = data.jag.ind.corr[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture IND-Corr.txt",
                                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ind.unc[[i]] <- jags(data = data.jag.ind.unc[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Mean Difference/Full FE-NMA MD-IMDOM Pattern-mixture IND-Uncorr.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ind.corr, meta.model = "FE", model = "independent", assumption = "correlated", measure = "MD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ind.unc, meta.model = "FE", model = "independent", assumption = "uncorrelated", measure = "MD", nt, ns, na, L))
        
      } else if(measure == "SMD"){
        
        param.jags.ind <- c("SMD", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ind.corr[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "M" = M[[i]], "var.phi" = var.misspar, "cov.phi" = cov.misspar, "D" = dir[i])
          
          data.jag.ind.unc[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "sigma" = sigma[[i]], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ind.corr[[i]] <- jags(data = data.jag.ind.corr[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture IND-Corr.txt",
                                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ind.unc[[i]] <- jags(data = data.jag.ind.unc[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Standardised Mean Difference/Full FE-NMA SMD-IMDOM Pattern-mixture IND-Uncorr.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ind.corr, meta.model = "FE", model = "independent", assumption = "correlated", measure = "SMD", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ind.unc, meta.model = "FE", model = "independent", assumption = "uncorrelated", measure = "SMD", nt, ns, na, L))
        
      } else {
        
        param.jags.ind <- c("LRoM", "SUCRA", "order", "phi")
        
        for(i in 1:L){
          
          data.jag.ind.corr[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "M" = M[[i]], "var.phi" = var.misspar, "cov.phi" = cov.misspar, "D" = dir[i])
          
          data.jag.ind.unc[[i]] <- list("y.o" = y0[J == i, ], "se.o" = se0[J == i, ], "m" = m[J == i, ], "N" = N[J == i, ], "t" = t[J == i, ], "na" = na[[i]], "nt" = nt[i], "ns" = ns[i], "ref" = ref[i], "prec.phi" = prec.misspar, "D" = dir[i])
          
          jagsfit.ind.corr[[i]] <- jags(data = data.jag.ind.corr[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture IND-Corr.txt",
                                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
          jagsfit.ind.unc[[i]] <- jags(data = data.jag.ind.unc[[i]], parameters.to.save = param.jags.ind, model.file = "./Pattern-mixture IMDOM & IMROM models/Model scripts/Full NMA models/Fixed-effect model/Ratio of Means/Full FE-NMA ROM-IMROM Pattern-mixture IND-Uncorr.txt",
                                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F, working.directory = getwd())
          
        }
        
        results <- list(collect.jags.results.Full.PM.continuous(jagsfit.ind.corr, meta.model = "FE", model = "independent", assumption = "correlated", measure = "ROM", nt, ns, na, L), 
                        collect.jags.results.Full.PM.continuous(jagsfit.ind.unc, meta.model = "FE", model = "independent", assumption = "uncorrelated", measure = "ROM", nt, ns, na, L))
        
      }
      
    } # End of models for FE
    
  } # End of RE/FE
    
return(results)
}

