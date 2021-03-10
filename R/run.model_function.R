# Perform network meta-analysis for an aggregate continuous outcome with missing participant data
#
# data: A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
# measure: Character string indicating the effect measure with values "MD", "SMD", or "ROM".
# assumption: Character string indicating the structure of the informative missingness parameter. Set assumption equal to one of the following: "HIE-COMMON", "HIE-TRIAL", "HIE-ARM", "IDE-COMMON", "IDE-TRIAL", "IDE-ARM", "IND-CORR", or "IND-UNCORR".
# mean.misspar: A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
# var.misspar: A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
# D: A binary number for the direction of the outcome. Set D = 1 for a positive outcome and D = 0 for a negative outcome.
# n.chains: Integer specifying the number of chains for the MCMC sampling; an argument of the jags function in R2jags.
# n.iter: Integer specifying the number of Markov chains for the MCMC sampling; an argument of the jags function in R2jags.
# n.burnin: Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the jags function in R2jags.
# n.thin: Integer specifying the thinning rate for the MCMC sampling; an argument of the jags function in R2jags.
#
# OUTPUT
# It returns sn R2jags output on the summaries of the posterior distribution, and the Gelman–Rubin convergence diagnostic of the following parameters:
# EM, the effect estimate of all possible comparisons of interventions.
# SUCRA, the surface under the cumulative ranking curve for each intervention.
# phi, the informative missingness parameter.
# delta, the underlying trial-specific effect estimate. For a multi-arm trial, we estimate T-1 trial-specific effect estimates, where T is the number of interventions in the trial.
# tau, The between-trial standard deviation assumed to be common for all observed comparisons.
# 
# FORMAT
# The columns of the data (a data frame) refer to the following ordered elements for a continuous outcome:
# t, an intervention identifier.
# y, the observed mean value of the outcome.
# se, the observed standard error of the outcome.
# m, the number of missing outcome data.
# c, the number of participants completing the assigned intervention.
# na, the number of compared interventions.
#
# Apart from na, all other elements appear in data as many times as the maximum number of interventions compared in a trial.

run.model <- function(data, measure, assumption, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin){

  ## Continuous: arm-level, wide-format dataset
  (y.obs <- data %>% dplyr::select(starts_with("y")))             # Observed mean value in each arm of every trial
  (se.obs <- data %>% dplyr::select(starts_with("se")))           # Observed standard error in each arm of every trial
  (mod <- data %>% dplyr::select(starts_with("m")))               # Number of missing participants in each arm of every trial
  (c <- data %>% dplyr::select(starts_with("c")))                 # Number of completers in each arm of every trial
  (sd.obs <- se.obs*c)                                            # Observed standard deviation in each arm of every trial
  (rand <- mod + c)                                               # Number of randomised participants in each arm of every trial
  (treat <- data %>% dplyr::select(starts_with("t")))             # Intervention studied in each arm of every trial
  na <- apply(treat, 1, function(x) length(which(!is.na(x))))     # Number of interventions investigated in every trial per network
  nt <- length(table(as.matrix(treat)))                           # Total number of interventions per network
  ns <- length(y.obs[, 1])                                        # Total number of included trials per network
  ref <- ifelse(nt > 2, which.max(table(as.matrix(treat))), 1)    # Reference intervention per network: the most frequently appeared intervention in the network
  # Trial-specific observed pooled standard deviation
  (sigma <- sqrt(apply((sd.obs^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - na)))

              
  ## Order by 'id of t1' < 'id of t1'
  y0 <- se0 <- m <- N <- t <- t0 <- treat
  for(i in 1:ns){
    t0[i, ] <- order(treat[i, ], na.last = T)
    y0[i, ] <- y.obs[i, order(t0[i, ], na.last = T)]
    se0[i, ] <- se.obs[i, order(t0[i, ], na.last = T)]
    m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
    N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
    t[i, ] <- sort(treat[i, ], na.last = T)
  }


  if((assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & !is.null(dim(mean.misspar))) {

    mean.misspar <- as.vector(mean.misspar)

  } else if((assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & is.null(dim(mean.misspar))) {

    mean.misspar <- rep(mean.misspar, 2)

  } else {

    mean.misspar <- mean.misspar
  }


  ## Information for the prior distribution on the missingness parameter (IMDOM or logIMROM)
  M <- ifelse(!is.na(y0), mean.misspar, NA)  # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
  prec.misspar <- 1/var.misspar
  psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
  cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)


  # Under the Independent structure with or without SMD as effect measure
  if (measure == "SMD" & assumption != "IND-CORR") {

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D)

  } else if (measure == "SMD" & assumption == "IND-CORR"){

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D)

  } else if (measure != "SMD" & assumption == "IND-CORR") {

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D)

  } else {

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D)

  }


  ## Condition for the hierarchical structure of the missingness parameter
  if (assumption == "HIE-COMMON" || assumption == "HIE-TRIAL" || assumption == "HIE-ARM") {

    param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "SUCRA", "mean.phi", "effectiveness")

  } else {

    param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "SUCRA", "phi", "effectiveness")

  }



  ## Run the Bayesian analysis
  jagsfit <- jags(data = data.jag, parameters.to.save = param.jags, model.file = paste0("./model/Full RE-NMA/Full RE-NMA_", measure, "_Pattern-mixture_", assumption, ".txt"),
                    n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = T)



  ## Obtain the posterior distribution of the necessary model paramters
  EM <- jagsfit$BUGSoutput$summary[1:(nt*(nt - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  EM.pred <- jagsfit$BUGSoutput$summary[(nt*(nt - 1)*0.5 + 1):(2*nt*(nt - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  EM.ref <- jagsfit$BUGSoutput$summary[paste0("EM.ref[", seq(1:nt)[-ref], "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  SUCRA <- jagsfit$BUGSoutput$summary[paste0("SUCRA[", seq(1:nt), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  delta <- jagsfit$BUGSoutput$summary[(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + 1):(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + sum(na - 1)), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  effectiveness <- jagsfit$BUGSoutput$summary[(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + sum(na - 1) + 2):(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + sum(na - 1) + 1 + nt*nt), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  pred.ref <- jagsfit$BUGSoutput$summary[paste0("pred.ref[", seq(1:nt)[-ref], "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  tau <- jagsfit$BUGSoutput$summary["tau", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  DIC <- jagsfit$BUGSoutput$DIC
  pD <- jagsfit$BUGSoutput$pD # pD = var(deviance)/2
  dev <- jagsfit$BUGSoutput$summary["deviance", "50%"]
  model.assessment <- data.frame(DIC, pD, dev)


  ## Conditions to obtain the posterior distribution of the missingness parameter
  if (assumption == "IDE-COMMON") {

    phi <- jagsfit$BUGSoutput$summary["phi", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

  } else if (assumption == "HIE-COMMON"){

    phi <- jagsfit$BUGSoutput$summary["mean.phi", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

  } else if (assumption == "IDE-TRIAL") {

    phi <- jagsfit$BUGSoutput$summary[paste0("phi[", seq(1:ns), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

  } else if (assumption == "HIE-TRIAL") {

    phi <- jagsfit$BUGSoutput$summary[paste0("mean.phi[", seq(1:ns), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

  } else if (assumption == "IDE-ARM") {

    phi <- jagsfit$BUGSoutput$summary[paste0("phi[", seq(1:nt), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

  } else if (assumption == "HIE-ARM") {

    phi <- jagsfit$BUGSoutput$summary[paste0("mean.phi[", seq(1:nt), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

  } else {

    phi <- jagsfit$BUGSoutput$summary[(nt*(nt - 1)*0.5 + nt + nt*nt + 1):(nt*(nt - 1)*0.5 + nt + nt*nt + 1 + sum(na) - 1), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

  }


  return(list(EM = EM, EM.ref = EM.ref, EM.pred = EM.pred, pred.ref = pred.ref, tau = tau, SUCRA = SUCRA, delta = delta, effectiveness = effectiveness, phi = phi, model.assessment = model.assessment, ref = ref))

}

