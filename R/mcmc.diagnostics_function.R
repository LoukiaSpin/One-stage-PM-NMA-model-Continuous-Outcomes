#' Markov Chain Monte Carlo Diagnostics
#'
#' @param par A vector of the names of the parameters to be monitored. Up to three parameters can be named. See, section \strong{return} in function \code{nma.continuous.full.model}.
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#' @param mean.misspar A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
#' @param D A binary number for the direction of the outcome. Set \code{direction = 0} for a positive outcome and \code{direction = 1} for a negative outcome.
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
#' \dontshow{load("./data/One-stage model_NMA Dataset.RData")}
#' @examples
#' ### Obtain the diagnostic plots and check convergence for all monitored parameters using the R.hat
#' mcmc.diagnostics(par = c("tau2", "EM[3,1]", "EM[3,2]"), data = data, measure = "ROM", assumption = "HIE-TRIAL", mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
mcmc.diagnostics <- function(par, data, measure, assumption, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin){

  ## Arm-level, wide-format dataset
  (y0 <- data %>% dplyr::select(starts_with("y")))            # Observed mean value in each arm of every trial
  (sd0 <- data %>% dplyr::select(starts_with("sd")))          # Observed standard deviation in each arm of every trial
  (m <- data %>% dplyr::select(starts_with("m")))             # Number of missing participants in each arm of every trial
  (c <- data %>% dplyr::select(starts_with("c")))             # Number completers in each arm of every trial
  (se0 <- sd0/sqrt(c))                                        # Observed standard error in each arm of every trial
  (N <- m + c)                                                # Number of randomised participants in each arm of every trial
  (t <- data %>% dplyr::select(starts_with("t")))             # Intervention studied in each arm of every trial
  na <- apply(t, 1, function(x) length(which(!is.na(x))))     # Number of interventions investigated in every trial per network
  nt <- length(table(as.matrix(t)))                           # Total number of interventions per network
  ns <- length(y0[, 1])                                       # Total number of included trials per network
  ref <- which.max(table(as.matrix(t)))                       # Reference intervention per network: the most frequently appeared intervention in the network
  # Trial-specific observed pooled standard deviation
  (sigma <- sqrt(apply((sd0^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - na)))



  ## Information for the prior distribution on the missingness parameter (IMDOM or logIMROM)
  M <- ifelse(!is.na(y0), mean.misspar, NA)  # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
  prec.misspar <- 1/var.misspar
  psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
  cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)


  ## Define the pattern-mixture models for various prior structures of the missingness parameter under the random-effects assumption
  if (assumption == "HIE-COMMON" || assumption == "HIE-TRIAL" || assumption == "HIE-ARM") {

    param.jags <- c("theta", "EM", "tau2", "SUCRA", "order", "mean.phi", "sd.phi", "effectiveness")

  } else {

    param.jags <- c("theta", "EM", "tau2", "SUCRA", "phi", "effectiveness")

  }

  if (measure == "SMD" & assumption != "IND-CORR") {

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D)

  } else if (measure == "SMD" & assumption == "IND-CORR"){

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D)

  } else if (measure != "SMD" & assumption == "IND-CORR") {

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D)

  } else {

    data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D)

  }


  jagsfit <- jags(data = data.jag, parameters.to.save = param.jags, model.file = paste0("./model/Full RE-NMA_", measure, "_Pattern-mixture_", assumption, ".txt"),
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F)

  EM <- jagsfit$BUGSoutput$summary[1:(nt*(nt - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  tausq <- jagsfit$BUGSoutput$summary["tau2", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  SUCRA <- jagsfit$BUGSoutput$summary[paste0("SUCRA[", seq(1:nt), "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]
  effectiveness <- jagsfit$BUGSoutput$summary[(nt*(nt - 1)*0.5 + nt + 1):(nt*(nt - 1)*0.5 + nt + nt*nt), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]

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

  # An HTML file with a panel of diagnostic plots per monitored paraemter
  mcmcplot <- mcmcplot(jagsfit.mcmc, parms = par)

  R.hat.max <- c(max(EM[, 5]), max(tausq[5]), max(SUCRA[, 5]), max(effectiveness[, 5]), max(phi[, 5]))
  conv <- rep(NA, length(R.hat.max))
  for(i in 1:length(R.hat.max)) {

    conv[i] <- ifelse(R.hat.max[i] < 1.1, "convergence achieved", "convergence issue")

  }

  # Check convergence for all monitored parameters using the Rhat
  convergence <- data.frame(R.hat.max, conv)
  rownames(convergence) <- c("EM", "tausq", "SUCRA", "effectiveness", "phi")
  colnames(convergence) <- c("R.hat max", "convergence status")

  return(list(convergence = convergence))
}

