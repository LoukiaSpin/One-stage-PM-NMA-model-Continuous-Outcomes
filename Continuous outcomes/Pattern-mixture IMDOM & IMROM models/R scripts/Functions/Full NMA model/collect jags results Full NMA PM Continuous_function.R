########################################################################################################
#                                                                                                      #
#               A function to obtain the posterior distribution of core NMA parameters                 # 
#                One-stage pattern-mixture model for aggregate continuous outcome data                 #
#         Informative Missingness Difference of Means & Informative Missingness Ratio of Means         #
#                               (Mavridis et al., 2015 - PMID: 25393541)                               #
#                                                                                                      #
########################################################################################################



## Necessary arguments (with explanations)
# name: the name assigned on the 'jags' function.
# meta.model: a character to indicate whether random-effects or fixed-effect network meta-analysis will be applied.
#             Possible values are "RE", and "FE" to indicate a random-effects and a fixed-effect network meta-analysis, respectively.
# model: a character to indicate the prior structure for the informative missingness odd ratio (IMOR).
#        Possible values are "identical", "hierarchical", "independent".
# assumption: a character to indicate the structural assumption about the missingness parameter.
#             Possible values are "common", "arm", "trial", and "correlated".
# measure: a character to indicate the effect measure for continuous outcomes. Possible values are "MD", "SMD", and "ROM".
# nt: total number of interventions per network (necessary even when a single network is analysed).
# ns: total number of included trials per network (necessary even when a single network is analysed).
# na: number of interventions investigated in every trial per network (necessary even when a single network is analysed).
# L: number of analysed networks (necessary even when a single network is analysed).

## Other abbreviations
# EM: effect measure, namely, mean difference, standardised mean difference, and ratio of means (in the logarithmic scale)
# MP: missingness parameter, namely, informative missingness difference of means, and informative missingness ratio of means (log scale)

collect.jags.results.Full.PM.continuous <- function(name, meta.model, model, assumption, measure, nt, ns, na, L){
  
  if(meta.model == "RE"){
    
    EM <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[1:(nt[i]*(nt[i] - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
    tau2 <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary["tau2", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
    SUCRA <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("SUCRA[", seq(1:nt[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
    order <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("order[", seq(1:nt[i]), "]"), c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
    
    if(model == "identical" & assumption == "common"){
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 3):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 3 + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary["phi", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, phi)
      names(results) <- c("IDE Common - EM", "IDE Common - theta", "IDE Common - tau2", "IDE Common - SUCRA", "IDE Common - order", "IDE Common - MP")
      
      if(measure == "MD"){

        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Common/MD IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Common/theta IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Common/tau2 IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Common/SUCRA IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Common/order IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Common/IMDOM IDE-Common.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Common/SMD IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Common/theta IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Common/tau2 IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Common/SUCRA IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Common/order IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Common/IMDOM IDE-Common.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Common/LROM IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Common/theta IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Common/tau2 IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Common/SUCRA IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Common/order IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Common/logIMROM IDE-Common.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "identical" & assumption == "arm"){
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 3*nt[i] + 2):(nt[i]*(nt[i] - 1)*0.5 + 3*nt[i] + 2 + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("phi[", seq(1:nt[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, phi)
      names(results) <- c("IDE Arm - EM", "IDE Arm - theta", "IDE Arm - tau2", "IDE Arm - SUCRA", "IDE Arm - order", "IDE Arm - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Arm/MD IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Arm/theta IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Arm/tau2 IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Arm/SUCRA IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Arm/order IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Arm/IMDOM IDE-Arm.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Arm/SMD IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Arm/theta IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Arm/tau2 IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Arm/SUCRA IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Arm/order IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Arm/IMDOM IDE-Arm.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Arm/LROM IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Arm/theta IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Arm/tau2 IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Arm/SUCRA IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Arm/order IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Arm/logIMROM IDE-Arm.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "identical" & assumption == "trial"){
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + ns[i] + 2):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + ns[i] + 2 + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("phi[", seq(1:ns[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, phi)
      names(results) <- c("IDE Trial - EM", "IDE Trial - theta", "IDE Trial - tau2", "IDE Trial - SUCRA", "IDE Trial - order", "IDE Trial - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Trial/MD IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Trial/theta IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Trial/tau2 IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Trial/SUCRA IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Trial/order IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Identical Trial/IMDOM IDE-Trial.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Trial/SMD IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Trial/theta IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Trial/tau2 IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Trial/SUCRA IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Trial/order IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Identical Trial/IMDOM IDE-Trial.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Trial/LROM IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Trial/theta IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Trial/tau2 IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Trial/SUCRA IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Trial/order IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Identical Trial/logIMROM IDE-Trial.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "hierarchical" & assumption == "common"){
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 4):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 4 + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      mean.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary["mean.phi", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      sd.phi<- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary["sd.phi", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, mean.phi, sd.phi)
      names(results) <- c("HIE Common - EM", "HIE Common - theta", "HIE Common - tau2", "HIE Common - SUCRA", "HIE Common - order", "HIE Common - mean MP", "HIE Common - sd MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common/MD HIE Common.txt", sep = "\t", quote = F)
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common/theta HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common/tau2 HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common/SUCRA HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common/order HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common/mean IMDOM HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Common/sd IMDOM HIE Common.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common/SMD HIE Common.txt", sep = "\t", quote = F)
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common/theta HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common/tau2 HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common/SUCRA HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common/order HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common/mean IMDOM HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Common/sd IMDOM HIE Common.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common/LROM HIE Common.txt", sep = "\t", quote = F)
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common/theta HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common/tau2 HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common/SUCRA HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common/order HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common/mean logIMROM HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Common/sd logIMROM HIE Common.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "hierarchical" & assumption == "arm"){
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 4*nt[i] + 2):(nt[i]*(nt[i] - 1)*0.5 + 4*nt[i] + 2 + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      mean.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("mean.phi[", seq(1:nt[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      sd.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("sd.phi[", seq(1:nt[i]), "]"), c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, mean.phi, sd.phi)
      names(results) <- c("HIE Arm - EM", "HIE Arm - theta", "HIE Arm - tau2", "HIE Arm - SUCRA", "HIE Arm - order", "HIE Arm - mean MP", "HIE Arm - sd MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm/MD HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm/theta HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm/tau2 HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm/SUCRA HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm/order HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm/mean IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Arm/sd IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm/SMD HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm/theta HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm/tau2 HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm/SUCRA HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm/order HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm/mean IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Arm/sd IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm/LROM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm/theta HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm/tau2 HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm/SUCRA HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm/order HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm/mean logIMROM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Arm/sd logIMROM HIE Arm.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "hierarchical" & assumption == "trial"){
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 2*ns[i] + 2):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 2*ns[i] + 2 + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      mean.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("mean.phi[", seq(1:ns[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      sd.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("sd.phi[", seq(1:ns[i]), "]"), c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, mean.phi, sd.phi)
      names(results) <- c("HIE Trial - EM", "HIE Trial - theta", "HIE Trial - tau2", "HIE Trial - SUCRA", "HIE Trial - order", "HIE Trial - mean MP", "HIE Trial - sd MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial/MD HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial/theta HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial/tau2 HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial/SUCRA HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial/order HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial/mean IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Hierarchical Trial/sd IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial/SMD HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial/theta HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial/tau2 HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial/SUCRA HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial/order HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial/mean IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Hierarchical Trial/sd IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial/LROM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial/theta HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial/tau2 HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial/SUCRA HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial/order HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial/mean logIMROM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Hierarchical Trial/sd logIMROM HIE Trial.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "independent" & assumption == "correlated"){
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 2 + sum(na[[i]])):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 2 + sum(na[[i]]) + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      phi <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1 + sum(na[[i]]) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, phi)
      names(results) <- c("IND Correlated - EM", "IND Correlated - theta", "IND Correlated - tau2", "IND Correlated - SUCRA", "IND Correlated - order", "IND Correlated - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Correlated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Correlated/MD IND-Corr.txt", sep = "\t", quote = F)
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Correlated/theta IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Correlated/tau2 IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Correlated/SUCRA IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Correlated/order IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Correlated/IMDOM IND-Corr.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Correlated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Correlated/SMD IND-Corr.txt", sep = "\t", quote = F)
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Correlated/theta IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Correlated/tau2 IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Correlated/SUCRA IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Correlated/order IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Correlated/IMDOM IND-Corr.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Correlated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Correlated/LROM IND-Corr.txt", sep = "\t", quote = F)
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Correlated/theta IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Correlated/tau2 IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Correlated/SUCRA IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Correlated/order IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Correlated/logIMROM IND-Corr.txt", sep = "\t", quote = F) 
        
      }
      
    } else { 
      
      theta <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 2 + sum(na[[i]])):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 2 + sum(na[[i]]) + sum(na[[i]] - 1) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      phi <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1 + sum(na[[i]]) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, theta, tau2, SUCRA, order, phi) 
      names(results) <- c("IND Uncorrelated - EM", "IND Uncorrelated - theta", "IND Uncorrelated - tau2", "IND Uncorrelated - SUCRA", "IND Uncorrelated - order", "IND Uncorrelated - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Uncorrelated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Uncorrelated/MD IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Uncorrelated/theta IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Uncorrelated/tau2 IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Uncorrelated/SUCRA IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Uncorrelated/order IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Mean Difference/Independent Uncorrelated/IMDOM IND-Unc.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Uncorrelated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Uncorrelated/SMD IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Uncorrelated/theta IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Uncorrelated/tau2 IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Uncorrelated/SUCRA IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Uncorrelated/order IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Standardised Mean Difference/Independent Uncorrelated/IMDOM IND-Unc.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Uncorrelated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Uncorrelated/LROM IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(theta, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Uncorrelated/theta IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(tau2, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Uncorrelated/tau2 IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Uncorrelated/SUCRA IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Uncorrelated/order IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Random-effects model/Ratio of Means/Independent Uncorrelated/logIMROM IND-Unc.txt", sep = "\t", quote = F) 
        
      }
      
    }
    
    ## Pattern-mixture models for various prior structures of the missingness parameter under the fixed-effect assumption
  } else {
    
    EM <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[1:(nt[i]*(nt[i] - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
    SUCRA <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("SUCRA[", seq(1:nt[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
    order <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("order[", seq(1:nt[i]), "]"), c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
    
    if(model == "identical" & assumption == "common"){
      
      phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary["phi", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, phi)
      names(results) <- c("IDE Common - EM", "IDE Common - SUCRA", "IDE Common - order", "IDE Common - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Common/MD IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Common/SUCRA IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Common/order IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Common/IMDOM IDE-Common.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Common/SMD IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Common/SUCRA IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Common/order IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Common/IMDOM IDE-Common.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Common/LROM IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Common/SUCRA IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Common/order IDE-Common.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Common/logIMROM IDE-Common.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "identical" & assumption == "arm"){
      
      phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("phi[", seq(1:nt[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, phi)
      names(results) <- c("IDE Arm - EM", "IDE Arm - SUCRA", "IDE Arm - order", "IDE Arm - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Arm/MD IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Arm/SUCRA IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Arm/order IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Arm/IMDOM IDE-Arm.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Arm/SMD IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Arm/SUCRA IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Arm/order IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Arm/IMDOM IDE-Arm.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Arm/LROM IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Arm/SUCRA IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Arm/order IDE-Arm.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Arm/logIMROM IDE-Arm.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "identical" & assumption == "trial"){
      
      phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("phi[", seq(1:ns[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, phi)
      names(results) <- c("IDE Trial - EM",  "IDE Trial - SUCRA", "IDE Trial - order", "IDE Trial - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Trial/MD IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Trial/SUCRA IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Trial/order IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Identical Trial/IMDOM IDE-Trial.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Trial/SMD IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Trial/SUCRA IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Trial/order IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Identical Trial/IMDOM IDE-Trial.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Trial/LROM IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Trial/SUCRA IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Trial/order IDE-Trial.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Identical Trial/logIMROM IDE-Trial.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "hierarchical" & assumption == "common"){
      
      mean.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary["mean.phi", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      sd.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary["sd.phi", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, mean.phi, sd.phi)
      names(results) <- c("HIE Common - EM", "HIE Common - SUCRA", "HIE Common - order", "HIE Common - mean MP", "HIE Common - sd MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Common/MD HIE Common.txt", sep = "\t", quote = F)
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Common/SUCRA HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Common/order HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Common/mean IMDOM HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Common/sd IMDOM HIE Common.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Common/SMD HIE Common.txt", sep = "\t", quote = F)
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Common/SUCRA HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Common/order HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Common/mean IMDOM HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Common/sd IMDOM HIE Common.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Common")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Common/LROM HIE Common.txt", sep = "\t", quote = F)
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Common/SUCRA HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Common/order HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Common/mean logIMROM HIE Common.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Common/sd logIMROM HIE Common.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "hierarchical" & assumption == "arm"){
      
      mean.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("mean.phi[", seq(1:nt[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      sd.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("sd.phi[", seq(1:nt[i]), "]"), c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, mean.phi, sd.phi)
      names(results) <- c("HIE Arm - EM", "HIE Arm - SUCRA", "HIE Arm - order", "HIE Arm - mean MP", "HIE Arm - sd MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Arm/MD HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Arm/SUCRA HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Arm/order HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Arm/mean IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Arm/sd IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Arm/SMD HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Arm/SUCRA HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Arm/order HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Arm/mean IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Arm/sd IMDOM HIE Arm.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Arm")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Arm/LROM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Arm/SUCRA HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Arm/order HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Arm/mean logIMROM HIE Arm.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Arm/sd logIMROM HIE Arm.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "hierarchical" & assumption == "trial"){
      
      mean.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("mean.phi[", seq(1:ns[i]), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      sd.phi <- do.call(rbind, lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[paste0("sd.phi[", seq(1:ns[i]), "]"), c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, mean.phi, sd.phi)
      names(results) <- c("HIE Trial - EM", "HIE Trial - SUCRA", "HIE Trial - order", "HIE Trial - mean MP", "HIE Trial - sd MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Trial/MD HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Trial/SUCRA HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Trial/order HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Trial/mean IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Hierarchical Trial/sd IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Trial/SMD HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Trial/SUCRA HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Trial/order HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Trial/mean IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Hierarchical Trial/sd IMDOM HIE Trial.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Trial")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Trial/LROM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Trial/SUCRA HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Trial/order HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(mean.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Trial/mean logIMROM HIE Trial.txt", sep = "\t", quote = F) 
        write.table(round(sd.phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Hierarchical Trial/sd logIMROM HIE Trial.txt", sep = "\t", quote = F) 
        
      }
      
    } else if(model == "independent" & assumption == "correlated"){
      
      phi <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1 + sum(na[[i]]) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, phi)
      names(results) <- c("IND Correlated - EM", "IND Correlated - SUCRA", "IND Correlated - order", "IND Correlated - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Correlated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Correlated/MD IND-Corr.txt", sep = "\t", quote = F)
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Correlated/SUCRA IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Correlated/order IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Correlated/IMDOM IND-Corr.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Correlated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Correlated/SMD IND-Corr.txt", sep = "\t", quote = F)
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Correlated/SUCRA IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Correlated/order IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Correlated/IMDOM IND-Corr.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Correlated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Correlated/LROM IND-Corr.txt", sep = "\t", quote = F)
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Correlated/SUCRA IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Correlated/order IND-Corr.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Correlated/logIMROM IND-Corr.txt", sep = "\t", quote = F)
        
      }
      
    } else {
      
      phi <- do.call(rbind,lapply(1:L, function(i) name[[i]]$BUGSoutput$summary[(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1):(nt[i]*(nt[i] - 1)*0.5 + 2*nt[i] + 1 + sum(na[[i]]) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
      
      results <- list(EM, SUCRA, order, phi) 
      names(results) <- c("IND Uncorrelated - EM", "IND Uncorrelated - SUCRA", "IND Uncorrelated - order", "IND Uncorrelated - MP")
      
      if(measure == "MD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Uncorrelated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Uncorrelated/MD IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Uncorrelated/SUCRA IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Uncorrelated/order IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Mean Difference/Independent Uncorrelated/IMDOM IND-Unc.txt", sep = "\t", quote = F) 
        
      } else if(measure == "SMD"){
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Uncorrelated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Uncorrelated/SMD IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Uncorrelated/SUCRA IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Uncorrelated/order IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Standardised Mean Difference/Independent Uncorrelated/IMDOM IND-Unc.txt", sep = "\t", quote = F) 
        
      } else {
        
        dir.create("./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Uncorrelated")
        write.table(round(EM, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Uncorrelated/LROM IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(SUCRA, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Uncorrelated/SUCRA IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(order, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Uncorrelated/order IND-Unc.txt", sep = "\t", quote = F) 
        write.table(round(phi, 3), file = "./Pattern-mixture IMDOM & IMROM models/Output/Data/Full NMA models/Fixed-effect model/Ratio of Means/Independent Uncorrelated/logIMROM IND-Unc.txt", sep = "\t", quote = F) 
        
      }
      
    }
    
  }
  
  return(list(results))
}