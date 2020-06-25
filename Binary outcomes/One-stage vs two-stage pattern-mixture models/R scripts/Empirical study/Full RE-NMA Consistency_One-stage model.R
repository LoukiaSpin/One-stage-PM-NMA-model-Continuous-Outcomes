#####################################################################################################################
#                                                                                                                   #
#               One-stage Bayesian random-effects network meta-analysis for aggregate binary outcomes               #                       # 
#                                 <Binomial likelihood, logit link, Random Effects>                                 # 
#                                 (Dias et al., 2013 in Appendix  - PMID: 23104435)                                 #    
#               Pattern-mixture model with Informative Missingness Odds Ratio under Missing at Random               #
#                                       (Turner et al., 2015 - PMID: 25809313)                                      #
#                                                                                                                   #
#####################################################################################################################


## Load libraries
library(R2jags); library(mcmcplots)


## Load datasets
load("./Dataset/One-stage model_NMA Dataset.RData")  # For Full NMA



## Run 'full NMA with consistency'
L <- 29
r <- list()
na <- t <- n <- m <- r
for(i in 1:L){
  r[[i]] <- one.stage.dataset.NMA[[i]][, c("r1", "r2", "r3", "r4")] 
  m[[i]] <- one.stage.dataset.NMA[[i]][, c("m1", "m2", "m3", "m4")]
  n[[i]] <- one.stage.dataset.NMA[[i]][, c("n1", "n2", "n3", "n4")] 
  t[[i]] <- one.stage.dataset.NMA[[i]][, c("t1", "t2", "t3", "t4")]  
  na[[i]] <- one.stage.dataset.NMA[[i]][, "na"]  
}
nt <- rep(0, L)
ref <- ns <- nt 
for(i in 1:L){
  nt[i] <- length(table(t[[i]]))     
  ref[i] <- which.max(table(t[[i]])) # Consider as reference the most frequently appeared internvetion
  ns[i] <- length(r[[i]][, 1])
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
  data.jag[[i]] <- list("r" = r[[i]], "m" = m[[i]], "n" = n[[i]], "t" = t[[i]], "na" = na[[i]], "nt" = nt[i], 
                        "ns" = ns[i], "ref" = ref[i], "mean.tau2" = mean.tausq[i], "prec.tau2" = prec.tausq[i])
  param.jags <- c("theta", "LOR.ref", "tausq")
  
  jagsfit.NMA[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = "./One-stage vs two-stage pattern-mixture models/Model scripts/Full RE-NMA One-stage IMOR Pattern-mixture model.txt",
                           n.chains = 3, n.iter = 100000, n.burnin = 10000, n.thin = 10, DIC = F, working.directory = getwd())
}
print(jagsfit.NMA)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken  # Time difference of 40.82558 mins



## Results for NMA parameters
(LOR <- do.call(rbind,lapply(1:L, function(i) jagsfit.NMA[[i]]$BUGSoutput$summary[1:(nt[i] - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
(theta <- do.call(rbind,lapply(1:L, function(i) jagsfit.NMA[[i]]$BUGSoutput$summary[(nt[i] - 1 + 2):(nt[i] - 1 + 1 + sum(na[[i]] - 1)), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
(tausq <- do.call(rbind,lapply(1:L, function(i) jagsfit.NMA[[i]]$BUGSoutput$summary["tausq", c("50%", "sd", "Rhat", "n.eff")])))
write.table(round(LOR, 4), file = "./One-stage vs two-stage pattern-mixture models/Output/Empirical study/LOR_One-stage_Results.txt", sep = "\t", quote = F) # LOR.ams[-134, ]
write.table(round(theta, 4), file = "./One-stage vs two-stage pattern-mixture models/Output/Empirical study/theta_One-stage_Results.txt", sep = "\t", quote = F) # theta[-430, ]
write.table(round(tausq, 4), file = "./One-stage vs two-stage pattern-mixture models/Output/Empirical study/tausq_One-stage_Results.txt", sep = "\t", quote = F)


