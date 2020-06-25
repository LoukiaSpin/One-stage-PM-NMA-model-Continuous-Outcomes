#####################################################################################################################
#                                                                                                                   #
#               Two-stage Bayesian random-effects network meta-analysis for aggregate binary outcomes               #                      
#                     <Normal likelihood, identity link, treatment differences, Random Effects>                     #           
#                          (Dias et al., 2013 in Appendix, Example 7(a)  - PMID: 23104435)                          #        
#               Pattern-mixture model with Informative Missingness Odds Ratio under Missing at Random               #
#                                       (White et al., 2008 - PMID: 17703496)                                       #
#                                                                                                                   #
#####################################################################################################################



## Load library
library(R2jags); library(runjags); library(MCMCpack)


## Load datasets
load("./Dataset/Two-stage model_NMA Dataset.RData")  # For Full NMA


## Run 'full NMA with consistency'
L <- 29
nt <- rep(0, 29)
y <- list()
se <- t <- V <- na <- y
ns2 <- ns3 <- ns4 <- ref <- ns20 <- ns30 <- ns40 <- nt
 
for(i in 1:L){
  y[[i]] <- two.stage.dataset.NMA[[i]][, c("y0", "y1", "y2", "y3")] 
  se[[i]] <- two.stage.dataset.NMA[[i]][, c("se0", "se1", "se2", "se3")]
  t[[i]] <- two.stage.dataset.NMA[[i]][, c("t1", "t2", "t3", "t4")] 
  V[[i]] <- two.stage.dataset.NMA[[i]][, 13] 
  na[[i]] <- two.stage.dataset.NMA[[i]][, 14]
  nt[i] <- length(table(two.stage.dataset.NMA[[i]][, 1:4]))   
  ref[i] <- which.max(table(two.stage.dataset.NMA[[i]][, 1:4]))
  ns20[i] <- table(two.stage.dataset.NMA[[i]][, 14] == 2)[2] 
  ns30[i] <- table(two.stage.dataset.NMA[[i]][, 14] == 3)[2]
  ns40[i] <- table(two.stage.dataset.NMA[[i]][, 14] == 4)[2]
  ns2[i] <- ifelse(is.na(table(two.stage.dataset.NMA[[i]][, 14] == 2)[2]), table(two.stage.dataset.NMA[[i]][, 14] == 2)[1], table(two.stage.dataset.NMA[[i]][, 14] == 2)[2])
  ns3[i] <- ifelse(is.na(ns30[i]), 0, ns30[i])
  ns4[i] <- ifelse(is.na(ns40[i]), 0, ns40[i])
}


## Predictive prior distributions for between-trial variance
(tau2.priors <- read.table("./Dataset/Empirical priors for between-trial variance.txt", head = T)) 
mean.tausq <- tau2.priors[, 2]
sd.tausq <- tau2.priors[, 3]
prec.tausq <- 1/sd.tausq^2


## Prepare parameters for JAGS
data.jag <- list()
jagsfit.NMA <- list()


## calculate time needed for all models
start.time <- Sys.time()
set.seed(123)
memory.limit(size = 40000)
for(i in 1:L){
  data.jag[[i]] <- list("y" = y[[i]], "se" = se[[i]], "t" = t[[i]], "V" = V[[i]], "nt" = nt[i], "ns2" = ns2[i], "ns3" = ns3[i], "ns4" = ns4[i],
                        "ref" = ref[i], "mean.tau2" = mean.tausq[i], "prec.tau2" = prec.tausq[i], "na" = na[[i]])
  param.jags <- c("delta", "LOR.ref", "tausq")
  
  jagsfit.NMA[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = "./One-stage vs two-stage pattern-mixture models/Model scripts/Full RE-NMA Two-stage IMOR Pattern-mixture model.txt",
                           n.chains = 3, n.iter = 100000, n.burnin = 10000, n.thin = 10, DIC = F, working.directory = getwd())
}
print(jagsfit.NMA)

end.time <- Sys.time()
time.taken <- end.time - start.time  
time.taken



## Results on NMA parameters  
(LOR <- do.call(rbind,lapply(1:L, function(i) jagsfit.NMA[[i]]$BUGSoutput$summary[1:(nt[i] - 1), c("mean", "sd", "2.5%", "97.5%")])))
(theta <- do.call(rbind,lapply(1:L, function(i) jagsfit.NMA[[i]]$BUGSoutput$summary[(nt[i] - 1 + 1):(nt[i] - 1 + sum(na[[i]] - 1)), c("mean", "sd", "2.5%", "97.5%")])))
(tausq <- do.call(rbind,lapply(1:L, function(i) jagsfit.NMA[[i]]$BUGSoutput$summary["tausq", c("50%", "sd")])))
write.table(round(LOR, 4), file = "./One-stage vs two-stage pattern-mixture models/Output/Empirical study/LOR_Two-stage_Results.txt", sep = "\t", quote = F) # LOR.ams[-134, ]
write.table(round(theta, 4), file = "./One-stage vs two-stage pattern-mixture models/Output/Empirical study/theta_Two-stage_Results.txt", sep = "\t", quote = F) # LOR.ams[-134, ]
write.table(round(tausq, 4), file = "./One-stage vs two-stage pattern-mixture models/Output/Empirical study/tausq_Two-stage_Results.txt", sep = "\t", quote = F)

