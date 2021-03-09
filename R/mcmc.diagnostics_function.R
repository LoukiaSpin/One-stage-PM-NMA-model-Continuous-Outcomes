#' Markov Chain Monte Carlo Diagnostics
#'
#' @param par A vector of the names of the parameters to be monitored. Up to three parameters can be named. See, section \strong{return} in function \code{nma.continuous.full.model}.
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#' @param mean.misspar A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a positive outcome and \code{D = 0} for a negative outcome.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#'
#' @return A panel of autocorrelation plots where the rows correspond to the chains and the columns correspond to the monitor parameters (maximum three).
#' Additionally, it uses the \code{\link[mcmcplots]{mcmcplot}} function to create an HTML file with a panel of diagnostic plots (trace, density, and autocorrelation) for each monitored parameter.
#'
#' @format See, function \code{nma.continuous.full.model}.
#'
#' @seealso \code{\link{mcmcplots}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. Stat Sci. 1992;7:457–472.
#'
#' \dontshow{load("netmodr/data/One-stage model_NMA Dataset.RData")}
#' @examples
#' ### Obtain the diagnostic plots and check convergence for all monitored parameters using the R.hat
#' mcmc.diagnostics(par = c("tau2", "EM[3,1]", "EM[3,2]"), data = data, measure = "ROM", assumption = "HIE-TRIAL", mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
mcmc.diagnostics <- function(par, data, measure, assumption, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin){



  if(measure == "MD" || measure == "SMD"|| measure == "ROM"){


    ## Continuous: arm-level, wide-format dataset
    (y.obs <- data %>% dplyr::select(starts_with("y")))             # Observed mean value in each arm of every trial
    (sd.obs <- data %>% dplyr::select(starts_with("sd")))           # Observed standard deviation in each arm of every trial
    (mod <- data %>% dplyr::select(starts_with("m")))               # Number of missing participants in each arm of every trial
    (c <- data %>% dplyr::select(starts_with("c")))                 # Number of completers in each arm of every trial
    (se.obs <- sd.obs/sqrt(c))                                      # Observed standard error in each arm of every trial
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



  } else {



    ## Binary: arm-level, wide-format dataset
    (event <- data %>% dplyr::select(starts_with("r")))             # Number of observed events in each arm of every trial
    (mod <- data %>% dplyr::select(starts_with("m")))               # Number of missing participants in each arm of every trial
    (rand <- data %>% dplyr::select(starts_with("n")))              # Number randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))             # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                           # Total number of interventions per network
    ns <- length(event[, 1])                                        # Total number of included trials per network
    ref <- ifelse(nt > 2, which.max(table(as.matrix(treat))), 1)    # Reference intervention per network: the most frequently appeared intervention in the network



    ## Order by 'id of t1' < 'id of t1'
    r <- m <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      r[i, ] <- event[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }


    ## Information for the prior distribution on log IMOR
    if((assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & !is.null(dim(mean.misspar))) {

      mean.misspar <- as.vector(mean.misspar)
      mean.misspar[1] <- ifelse(mean.misspar[1] == 0, 0.0001, mean.misspar[1])
      mean.misspar[2] <- ifelse(mean.misspar[2] == 0, 0.0001, mean.misspar[2])

    } else if((assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & is.null(dim(mean.misspar))) {

      mean.misspar <- rep(ifelse(mean.misspar == 0, 0.0001, mean.misspar), 2)

    } else if(assumption != "HIE-ARM" || assumption != "IDE-ARM" ) {

      mean.misspar <- ifelse(mean.misspar == 0, 0.0001, mean.misspar)
    }


    M <- ifelse(!is.na(r), mean.misspar, NA)   # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
    prec.misspar <- 1/var.misspar
    psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
    cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)



    ## Condition for the Independent structure
    if (assumption != "IND-CORR") {

      data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D)

    } else {

      data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D)

    }


  }



  ## Condition for the hierarchical structure of the missingness parameter
  if (assumption == "HIE-COMMON" || assumption == "HIE-TRIAL" || assumption == "HIE-ARM") {

    param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "SUCRA", "mean.phi", "effectiveness")

  } else {

    param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "SUCRA", "phi", "effectiveness")

  }



  ## Run the Bayesian analysis
  jagsfit <- jags(data = data.jag, parameters.to.save = param.jags, model.file = paste0("./model/Full RE-NMA/Full RE-NMA_", measure, "_Pattern-mixture_", assumption, ".txt"),
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F)



  ## Obtain the posterior distribution of the necessary model paramters
  EM <- jagsfit$BUGSoutput$summary[1:(nt*(nt - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  EM.pred <- jagsfit$BUGSoutput$summary[(nt*(nt - 1)*0.5 + 1):(2*nt*(nt - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  EM.ref <- jagsfit$BUGSoutput$summary[paste0("EM.ref[", seq(1:nt)[-ref], "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  SUCRA <- jagsfit$BUGSoutput$summary[paste0("SUCRA[", seq(1:nt), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  delta <- jagsfit$BUGSoutput$summary[(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + 1):(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + sum(na - 1)), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  effectiveness <- jagsfit$BUGSoutput$summary[(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + sum(na - 1) + 1):(2*nt*(nt - 1)*0.5 + (nt - 1) + nt + sum(na - 1) + nt*nt), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  pred.ref <- jagsfit$BUGSoutput$summary[paste0("pred.ref[", seq(1:nt)[-ref], "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  tau <- jagsfit$BUGSoutput$summary["tau", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]


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



  jagsfit.mcmc <- as.mcmc(jagsfit)
  ## A panel of autocorrelation plots for each chain and every monitored parameter
  autocorrelation <- par(mfrow = c(3, n.chains))
  for(i in 1:n.chains){

    autplot1(jagsfit.mcmc[, par[1]], chain = i, main = paste(par[1], "-", "chain", i))
    autplot1(jagsfit.mcmc[, par[2]], chain = i, main = paste(par[2], "-","chain", i))
    autplot1(jagsfit.mcmc[, par[3]], chain = i, main = paste(par[3], "-","chain", i))

  }

  ## An HTML file with a panel of diagnostic plots per monitored paraemter
  mcmcplot <- mcmcplot(jagsfit.mcmc, parms = par)


  if(assumption == "IDE-COMMON" || assumption == "HIE-COMMON"){

    R.hat.max <- c(max(EM[, 5]), max(EM.pred[, 5]), max(delta[, 5]), max(tau[5]), max(SUCRA[, 5]), max(effectiveness[, 5]), phi[5])

  } else {

    R.hat.max <- c(max(EM[, 5]), max(EM.pred[, 5]), max(delta[, 5]), max(tau[5]), max(SUCRA[, 5]), max(effectiveness[, 5]), max(phi[, 5]))

  }

  conv <- rep(NA, length(R.hat.max))
  for(i in 1:length(R.hat.max)) {

    conv[i] <- ifelse(R.hat.max[i] < 1.1, "achieved", "failed")

  }

  # Check convergence for all monitored parameters using the Rhat
  convergence <- data.frame(R.hat.max, conv)
  rownames(convergence) <- c("EM", "Pred", "delta", "tau", "SUCRA", "effectiveness", "phi")
  colnames(convergence) <- c("R.hat max", "convergence status")

  return(list(convergence = convergence))
}

