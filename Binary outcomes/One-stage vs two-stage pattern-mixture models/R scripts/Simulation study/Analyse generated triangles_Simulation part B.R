########################################################################
#                                                                      #
#           Code to analyse the simulated triangles - Part B           #
#                                                                      #
########################################################################



# Initially, load the following necessary packages for the simulation: R2jags, MCMCpack, and mcmcplots:

list.of.packages <- c("R2jags", "MCMCpack", "mcmcplots")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)  

# Then, load the function to apply the Taylor series for the first-stage of the two-stage pattern-mixture model:

source("./One-stage vs two-stage pattern-mixture models/R scripts/Simulation study/Two-stage PM model_Function for First-stage.R") 
       
# Load the 1,000 generated triangle networks with missing outcome data:

load("./One-stage vs two-stage pattern-mixture models/R scripts/Simulation study/1000 simulated two-arm networks.RData")

# Description of true parameters with their pre-specified scenarios.

sim <- 1000                         # Number of simulations
mOP <- 1.5                          # True odds ratio on 'Old intervention versus Placebo' 
mNP <- 2                            # True odds ratio on 'New intervention versus Placebo'
size <- c("<50", ">100")            # Size of the trials
event <- c("low", "frequent")       # Event frequency
tau2 <- c("small", "substantial")   # Extent of between-trial variance (tau2)
pm <- c("moderate", "large")        # Extent of missing outcome data 

# The scenarios considered in the present simulation study:

(compile <- expand.grid(list(size = size, event = event, tau2 = tau2, pm = pm)))

# Then, the 1,000 generated arm-level networks are analyzed using one-stage and two-stage pattern-mixture network meta-analysis (NMA). 
# The NMA estimates of interest are the NMA log odds ratio of an event for all comparisons, the common between-trial variance, and 
# the SUCRA value of new and old intervention and placebo. For both models, we applied 80,000 MCMC iterations (n.iter) with 20,000 burn-in (n.burnin), 
# 3 parallel chains (n.chains) and thinning (n.thin) equal to 10.

n.iter <- 80000                     # Number of MCMC iterations
n.burnin <- 20000                   # Number of burn-in samples
n.chains <- 3                       # Number of parallel chains
n.thin <- 10                        # Thinning

# The following code indicates the number of MCMC samples to be saved. We need this information later in the code to obtain the posterior distribution of each parameter 
# for every simulation and scenario.

(plus <- ((n.iter - n.burnin)/n.thin)*n.chains)   

# The function binmod.network.analysis performs the necessary analysis for each model and each scenario. 
# The parameters lor32, lor31 and lor21 refer to the NMA log odds ratio of event for new versus old intervention, for new intervention versus palcebo, and 
# for old intervention versus placebo, respectively, whereas the parameters sucraN, sucraO, and sucraP refer to the SUCRA value for the new intervention, old intervention, 
# and palcebo, respectively.

binmod.network.analysis <- function(dataset, size, event, tau2, pm, sim){ 
  
  counter1 <- 0                        # Returning the % of completed simulations
  index <- 1
  
  
  
  #################################################################################################################
  ###  CREATE VECTOR FOR EACH STATISTIC TO STORE POSTERIOR DISTRIBUTION FOR ALL SIMULATIONS AND FOR EACH MODEL  ###
  #################################################################################################################
  
  lor32.tempmat.one.pm <- lor31.tempmat.one.pm <- lor21.tempmat.one.pm <- rep(NA, sim*plus) 
  lor32.tempmat.two.pm <- lor31.tempmat.two.pm <- lor21.tempmat.two.pm <- rep(NA, sim*plus) 
  tau2.tempmat.one.pm <- sucraN.tempmat.one.pm <- sucraO.tempmat.one.pm <- sucraP.tempmat.one.pm <- rep(NA, sim*plus) 
  tau2.tempmat.two.pm <- sucraN.tempmat.two.pm <- sucraO.tempmat.two.pm <- sucraP.tempmat.two.pm <- rep(NA, sim*plus) 
  converged.lor32.tempmat.one.pm <- converged.lor31.tempmat.one.pm <- converged.lor21.tempmat.one.pm <- rep(NA, sim)
  converged.lor32.tempmat.two.pm <- converged.lor31.tempmat.two.pm <- converged.lor21.tempmat.two.pm <- rep(NA, sim)
  converged.tau2.tempmat.one.pm <- converged.tau2.tempmat.two.pm <- rep(NA, sim)  

  
    
  #######################################################################################
  ###  CREATE DATAFRAME FOR EACH MODEL TO STORE SUMMARY STATISTICS FOR EACH SCENARIO  ###
  #######################################################################################
  
  results.one.pm <- results.two.pm <- matrix(NA, length(compile[, 1]), 22)    
  colnames(results.one.pm) <- colnames(results.two.pm) <- c("median_LOR32","lower_LOR32","upper_LOR32","median_LOR31",
                                                            "lower_LOR31","upper_LOR31","median_LOR21","lower_LOR21",
                                                            "upper_LOR21","median_tau2","lower_tau2","upper_tau2","mean_sucraN",
                                                            "sd_sucraN","mean_sucraO","sd_sucraO","mean_sucraP","sd_sucraP",
                                                            "converged_LOR32","converged_LOR31","converged_LOR21",
                                                            "converged_tau2")
  
  
  
  for(i in 1:length(compile[, 1])){  # LOOP For SCENARIOS
    
    index.sim <- 1
    index.convergence<-1
    
    for(l in 1:sim){    # LOOP FOR SIMULATED DATASET
      sim.data <- mat[[i]][,, l]
      r <- data.matrix(sim.data[, c("r2", "r1")])        # Number of events per arm (the first arm is the baseline in the code)
      m <- data.matrix(sim.data[, c("m2", "m1")])        # Number of missing outcome data per arm
      n <- data.matrix(sim.data[, c("n2", "n1")])        # Number of randomized participants per arm
      t <- data.matrix(sim.data[, c("t2", "t1")])        # Intervention studied per arm (1: Placebo; 2: Old; 3: New)
      nt <- max(t, na.rm = TRUE)                         # Number of interventions (three, since we are studying a triangle)
      ns <- nrow(r)                                      # Number of included trials
      na <- data.matrix(sim.data[,"na"])                 # Number of arms (two, since we are focusing on two-arm trials)
      ref <- 1                                           # Placebo as the reference intervention
      
      # Obtain within-trial logOR and SElogOR - For each trial, it calculates the first versus second treatment arm
      dataset.averageMAR <- Taylor.IMOR(data = sim.data[, 1:9], delta1 = 0, delta2 = 0, var.delta1 = 1, var.delta2 = 1, rho = 0)
      y <- cbind(rep(NA, length(dataset.averageMAR$logOR)), dataset.averageMAR$logOR)       # Within-trial LORs 
      se <- cbind(rep(NA, length(dataset.averageMAR$SElogOR)), dataset.averageMAR$SElogOR)  # Within-trial standard error of LOR 
      
      
      
      #########################################################     
      ###  DEFINE ALL NECESSARY DATA TO RUN FULL NMA MODEL  ###
      #########################################################
      
      Data.one   <- list("r" = r, "m" = m, "n" = n, "t" = t, "na" = na[, 1], "nt" = nt, "ns" = ns, "ref" = ref, "mean.tau2" = -2.06, "prec.tau2" = 0.44) # Symptoms/signs reflecting end of condition (Pharma vs PBO) [PMID: 25475839]
      Data.two <- list("y" = y, "se" = se, "t" = t, "na" = na[, 1], "nt" = nt, "ns" = ns, "ref" = ref, "mean.tau2" = -2.06, "prec.tau2" = 0.44)          # Symptoms/signs reflecting end of condition (Pharma vs PBO) [PMID: 25475839]
      
      
      #########################################
      ###  RUN FULL RE-NMA MODEL UNDER MAR  ###
      #########################################
      
      jagsfit.one.pm <- jags(data = Data.one, parameters.to.save = c("tausq", "LOR", "SUCRA"), n.chains = n.chains, 
                             n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, 
                             model.file = "./One-stage vs two-stage pattern-mixture models/Model scripts/Full RE-NMA One-stage IMOR Pattern-mixture model.txt", DIC = F, working.directory = getwd())
      jagsfit.two.pm <- jags(data = Data.two, parameters.to.save = c("tausq", "LOR", "SUCRA"), n.chains = n.chains, 
                             n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, 
                             model.file = "./One-stage vs two-stage pattern-mixture models/Model scripts/Full RE-NMA Two-stage IMOR Pattern-mixture model_Two-arm trials.txt", DIC = F, working.directory = getwd())
      
      
      #########################################################################################################
      ###  Collect results on (i) log OR & 95% CrI, (ii) tau & 95% CrI, and (iii) SUCRA for each treatment  ###
      #########################################################################################################
      


      ## Collect posterior distribution of NMA LOR for all pairwise comparisons, tau2 and SUCRAs
      # One-stage pattern-mixture model                                                                                             
      lor32.tempmat.one.pm[c(index.sim:(index.sim + plus - 1))] <- c(as.mcmc(jagsfit.one.pm)[[1]][, "LOR[3,2]"],   
                                                                     as.mcmc(jagsfit.one.pm)[[2]][, "LOR[3,2]"], 
                                                                     as.mcmc(jagsfit.one.pm)[[3]][, "LOR[3,2]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)


      lor31.tempmat.one.pm[c(index.sim:(index.sim + plus - 1))] <- c(as.mcmc(jagsfit.one.pm)[[1]][, "LOR[3,1]"],  
                                                                     as.mcmc(jagsfit.one.pm)[[2]][, "LOR[3,1]"], 
                                                                     as.mcmc(jagsfit.one.pm)[[3]][, "LOR[3,1]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      lor21.tempmat.one.pm[c(index.sim:(index.sim + plus - 1))] <- c(as.mcmc(jagsfit.one.pm)[[1]][, "LOR[2,1]"],  
                                                                     as.mcmc(jagsfit.one.pm)[[2]][, "LOR[2,1]"], 
                                                                     as.mcmc(jagsfit.one.pm)[[3]][, "LOR[2,1]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      tau2.tempmat.one.pm[c(index.sim:(index.sim + plus - 1))] <- c(as.mcmc(jagsfit.one.pm)[[1]][, "tausq"],
                                                                    as.mcmc(jagsfit.one.pm)[[2]][, "tausq"],
                                                                    as.mcmc(jagsfit.one.pm)[[3]][, "tausq"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      sucraN.tempmat.one.pm[c(index.sim:(index.sim + plus - 1))] <- c(as.mcmc(jagsfit.one.pm)[[1]][, "SUCRA[3]"],
                                                                      as.mcmc(jagsfit.one.pm)[[2]][, "SUCRA[3]"],
                                                                      as.mcmc(jagsfit.one.pm)[[3]][, "SUCRA[3]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      sucraO.tempmat.one.pm[c(index.sim:(index.sim + plus - 1))] <- c(as.mcmc(jagsfit.one.pm)[[1]][, "SUCRA[2]"],
                                                                      as.mcmc(jagsfit.one.pm)[[2]][, "SUCRA[2]"],
                                                                      as.mcmc(jagsfit.one.pm)[[3]][, "SUCRA[2]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      sucraP.tempmat.one.pm[c(index.sim:(index.sim + plus - 1))] <- c(as.mcmc(jagsfit.one.pm)[[1]][, "SUCRA[1]"],
                                                                      as.mcmc(jagsfit.one.pm)[[2]][, "SUCRA[1]"],
                                                                      as.mcmc(jagsfit.one.pm)[[3]][, "SUCRA[1]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      # If R-hat > 1, no convergence for logOR (for 3 vs 2 is found in the 3nd line of the summary results)
      converged.lor32.tempmat.one.pm[index.convergence] <- ifelse(jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1, NA, 1) 

      # If R-hat > 1, no convergence for logOR (for 3 vs 1 is found in the 2nd line of the summary results)   
      converged.lor31.tempmat.one.pm[index.convergence] <- ifelse(jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1, NA, 1) 

      # If R-hat > 1, no convergence for logOR (for 2 vs 1 is found in the 1rst line of the summary results)
      converged.lor21.tempmat.one.pm[index.convergence] <- ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1, NA, 1) 

      # If R-hat > 1, no convergence for tau2 (it is found in the 7th line of the summary results)
      converged.tau2.tempmat.one.pm[index.convergence] <- ifelse(jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, NA, 1)     
      
      
      # Two-stage pattern-mixture model                                                                                             
      lor32.tempmat.two.pm[c(index.sim:(index.sim + plus))] <- c(as.mcmc(jagsfit.two.pm)[[1]][, "LOR[3,2]"],
                                                                 as.mcmc(jagsfit.two.pm)[[2]][, "LOR[3,2]"],
                                                                 as.mcmc(jagsfit.two.pm)[[3]][, "LOR[3,2]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      lor31.tempmat.two.pm[c(index.sim:(index.sim + plus))] <- c(as.mcmc(jagsfit.two.pm)[[1]][, "LOR[3,1]"],
                                                                 as.mcmc(jagsfit.two.pm)[[2]][, "LOR[3,1]"],
                                                                 as.mcmc(jagsfit.two.pm)[[3]][, "LOR[3,1]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      lor21.tempmat.two.pm[c(index.sim:(index.sim + plus))] <- c(as.mcmc(jagsfit.two.pm)[[1]][, "LOR[2,1]"],
                                                                 as.mcmc(jagsfit.two.pm)[[2]][, "LOR[2,1]"],
                                                                 as.mcmc(jagsfit.two.pm)[[3]][, "LOR[2,1]"])*
     rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      tau2.tempmat.two.pm[c(index.sim:(index.sim + plus))] <- c(as.mcmc(jagsfit.two.pm)[[1]][, "tausq"],
                                                                as.mcmc(jagsfit.two.pm)[[2]][, "tausq"],
                                                                as.mcmc(jagsfit.two.pm)[[3]][, "tausq"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      sucraN.tempmat.two.pm[c(index.sim:(index.sim + plus))] <- c(as.mcmc(jagsfit.two.pm)[[1]][, "SUCRA[3]"],
                                                                  as.mcmc(jagsfit.two.pm)[[2]][, "SUCRA[3]"],
                                                                  as.mcmc(jagsfit.two.pm)[[3]][, "SUCRA[3]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      sucraO.tempmat.two.pm[c(index.sim:(index.sim + plus))] <- c(as.mcmc(jagsfit.two.pm)[[1]][, "SUCRA[2]"],
                                                                  as.mcmc(jagsfit.two.pm)[[2]][, "SUCRA[2]"],
                                                                  as.mcmc(jagsfit.two.pm)[[3]][, "SUCRA[2]"])*
      rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      sucraP.tempmat.two.pm[c(index.sim:(index.sim + plus))] <- c(as.mcmc(jagsfit.two.pm)[[1]][, "SUCRA[1]"],
                                                                  as.mcmc(jagsfit.two.pm)[[2]][, "SUCRA[1]"],
                                                                  as.mcmc(jagsfit.two.pm)[[3]][, "SUCRA[1]"])*
     rep(ifelse(jagsfit.one.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1 || 
                 jagsfit.one.pm$BUGSoutput$summary[7, "Rhat"] > 1.1 || jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, 
                 NA, 1), plus)

      # If R-hat > 1, no convergence for logOR (for 3 vs 2 is found in the 3nd line of the summary results)
      converged.lor32.tempmat.two.pm[index.convergence] <- ifelse(jagsfit.two.pm$BUGSoutput$summary[3, "Rhat"] > 1.1, NA, 1) 

      # If R-hat > 1, no convergence for logOR (for 3 vs 1 is found in the 2nd line of the summary results)  
      converged.lor31.tempmat.two.pm[index.convergence] <- ifelse(jagsfit.two.pm$BUGSoutput$summary[2, "Rhat"] > 1.1, NA, 1) 

      # If R-hat > 1, no convergence for logOR (for 2 vs 1 is found in the 1rst line of the summary results) 
      converged.lor21.tempmat.two.pm[index.convergence] <- ifelse(jagsfit.two.pm$BUGSoutput$summary[1, "Rhat"] > 1.1, NA, 1) 

      # If R-hat > 1, no convergence for tau2 (it is found in the 7th line of the summary results)  
      converged.tau2.tempmat.two.pm[index.convergence] <- ifelse(jagsfit.two.pm$BUGSoutput$summary[7, "Rhat"] > 1.1, NA, 1)    
      
      index.sim <- index.sim + plus
      index.convergence <- index.convergence + 1
      counter1 <- counter1 + 1 
      
      print(paste("Percentage completed: ", round((counter1/(sim*length(size)*length(event)*length(tau2)*length(pm)))*100, 2),"%")) 
      
    } # END of SIMULATIONS LOOP
    
    
    ############################################################################################################################
    ###  Collect results on (i) mean and 95% CrI of LOR, (ii) median and 95% CrI of tau (iii)  mean SUCRA per intervention   ###
    ############################################################################################################################
    
    ## Collect mean NMA-logOR for all possible comparisons and lower-upper bound of 95% CrIs 
    # One-stage pattern-mixture model
    results.one.pm[index, 1] <- mean(lor32.tempmat.one.pm, na.rm = T)    # New vs Old intervention 
    results.one.pm[index, 2] <- quantile(lor32.tempmat.one.pm, prob = 0.025, na.rm = T)
    results.one.pm[index, 3] <- quantile(lor32.tempmat.one.pm, prob = 0.975, na.rm = T)
    results.one.pm[index, 4] <- mean(lor31.tempmat.one.pm, na.rm = T)    # New vs PBO   
    results.one.pm[index, 5] <- quantile(lor31.tempmat.one.pm, prob = 0.025, na.rm = T)
    results.one.pm[index, 6] <- quantile(lor31.tempmat.one.pm, prob = 0.975, na.rm = T)
    results.one.pm[index, 7] <- mean(lor21.tempmat.one.pm, na.rm = T)    # Old vs PBO   
    results.one.pm[index, 8] <- quantile(lor21.tempmat.one.pm, prob = 0.025, na.rm = T)
    results.one.pm[index, 9] <- quantile(lor21.tempmat.one.pm, prob = 0.975, na.rm = T)
    # Two-stage pattern-mixture model
    results.two.pm[index, 1] <- mean(lor32.tempmat.two.pm, na.rm = T)     # New vs Old intervention   
    results.two.pm[index, 2] <- quantile(lor32.tempmat.two.pm, prob = 0.025, na.rm = T)
    results.two.pm[index, 3] <- quantile(lor32.tempmat.two.pm, prob = 0.975, na.rm = T)
    results.two.pm[index, 4] <- mean(lor31.tempmat.two.pm, na.rm = T)     # New vs PBO   
    results.two.pm[index, 5] <- quantile(lor31.tempmat.two.pm, prob = 0.025, na.rm = T)
    results.two.pm[index, 6] <- quantile(lor31.tempmat.two.pm, prob = 0.975, na.rm = T)
    results.two.pm[index, 7] <- mean(lor21.tempmat.two.pm, na.rm = T)     # Old vs PBO   
    results.two.pm[index, 8] <- quantile(lor21.tempmat.two.pm, prob = 0.025, na.rm = T)
    results.two.pm[index, 9] <- quantile(lor21.tempmat.two.pm, prob = 0.975, na.rm = T)
    
    
    ## Collect median between-study variance and lower-upper bound of 95% CrI
    # One-stage pattern-mixture model
    results.one.pm[index, 10] <- median(tau2.tempmat.one.pm, na.rm = T)
    results.one.pm[index, 11] <- quantile(tau2.tempmat.one.pm, prob = 0.02, na.rm = T)
    results.one.pm[index, 12] <- quantile(tau2.tempmat.one.pm, prob = 0.975, na.rm = T)
    # Two-stage pattern-mixture model
    results.two.pm[index, 10] <- median(tau2.tempmat.two.pm, na.rm = T)
    results.two.pm[index, 11] <- quantile(tau2.tempmat.two.pm, prob = 0.025, na.rm = T)
    results.two.pm[index, 12] <- quantile(tau2.tempmat.two.pm, prob = 0.975, na.rm = T)
    
    
    ## Collect mean SUCRA value and standard deviation for NEW 
    # One-stage pattern-mixture model
    results.one.pm[index, 13] <- mean(sucraN.tempmat.one.pm, na.rm = T)
    results.one.pm[index, 14] <- sd(sucraN.tempmat.one.pm, na.rm = T)
    # Two-stage pattern-mixture model
    results.two.pm[index, 13] <- mean(sucraN.tempmat.two.pm, na.rm = T)  
    results.two.pm[index, 14] <- sd(sucraN.tempmat.two.pm, na.rm = T)
    
    
    ## Collect mean SUCRA value and standard deviation for OLD 
    # One-stage pattern-mixture model
    results.one.pm[index, 15] <- mean(sucraO.tempmat.one.pm, na.rm = T)  
    results.one.pm[index, 16] <- sd(sucraO.tempmat.one.pm, na.rm = T)
    # Two-stage pattern-mixture model
    results.two.pm[index, 15] <- mean(sucraO.tempmat.two.pm, na.rm = T)  
    results.two.pm[index, 16] <- sd(sucraO.tempmat.two.pm, na.rm = T)
    
    
    ## Collect mean SUCRA value and standard deviation for PBO 
    # One-stage pattern-mixture model
    results.one.pm[index, 17] <- mean(sucraP.tempmat.one.pm, na.rm = T)  
    results.one.pm[index, 18] <- sd(sucraP.tempmat.one.pm, na.rm = T)
    # Two-stage pattern-mixture model
    results.two.pm[index, 17] <- mean(sucraP.tempmat.two.pm, na.rm = T)  
    results.two.pm[index, 18] <- sd(sucraP.tempmat.two.pm, na.rm = T)
    
    
    ## Collect Gelman-Rubin's R-hat (R-hat < 1.1 -> convergence) 
    # One-stage pattern-mixture model
    results.one.pm[index, 19] <- sum(converged.lor32.tempmat.one.pm, na.rm = T)/sim    # % of simulations with converged LOR32 
    results.one.pm[index, 20] <- sum(converged.lor31.tempmat.one.pm, na.rm = T)/sim    # % of simulations with converged LOR31
    results.one.pm[index, 21] <- sum(converged.lor21.tempmat.one.pm, na.rm = T)/sim    # % of simulations with converged LOR21
    results.one.pm[index, 22] <- sum(converged.tau2.tempmat.one.pm, na.rm = T)/sim     # % of simulations with converged tau2

    ## Collect Gelman-Rubin's R-hat (R-hat < 1.1 -> convergence) 
    # Two-stage pattern-mixture model
    results.two.pm[index, 19] <- sum(converged.lor32.tempmat.two.pm, na.rm = T)/sim    # % of simulations with converged LOR32
    results.two.pm[index, 20] <- sum(converged.lor31.tempmat.two.pm, na.rm = T)/sim    # % of simulations with converged LOR31
    results.two.pm[index, 21] <- sum(converged.lor21.tempmat.two.pm, na.rm = T)/sim    # % of simulations with converged LOR21
    results.two.pm[index, 22] <- sum(converged.tau2.tempmat.two.pm, na.rm = T)/sim     # % of simulations with converged tau2
    
    index <- index + 1
    
  }   # END of SCENARIOS LOOP
  final.results <- list(results.one.pm, results.two.pm)
  
  return(final.results)
} 

# Calculate the time needed for the whole simulation.

start.time <- Sys.time()

set.seed(123)
memory.limit(size = 16000)
mtc.results <- binmod.network.analysis(mat, size, event, tau2, pm, sim)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



