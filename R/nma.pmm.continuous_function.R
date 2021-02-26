#' Perform network meta-analysis for an aggregate continuous outcome with missing participant data
#'
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
#' @return An R2jags output on the summaries of the posterior distribution, and the Gelman–Rubin convergence diagnostic of the following parameters:
#' \describe{
#'  \item{\code{EM}}{The effect estimate of all possible comparisons of interventions.}
#'  \item{\code{SUCRA}}{The surface under the cumulative ranking curve for each intervention.}
#'  \item{\code{phi}}{The informative missingness parameter.}
#'  \item{\code{theta}}{The trial-specific effect estimate. For a multi-arm trial, we estimate \emph{T-1} trial-specific effect estimates, where \emph{T} is the number of interventions in the trial.}
#'  \item{\code{tausq}}{The between-trial variance assumed to be common for all observed comparisons.}
#' }
#'
#' @format The columns of the data frame \code{data} refer to the following ordered elements:
#' \describe{
#'  \item{\strong{t}}{An intervention identifier.}
#'  \item{\strong{y}}{The observed mean value of the outcome.}
#'  \item{\strong{sd}}{The observed standard deviation of the outcome.}
#'  \item{\strong{m}}{The number of missing outcome data.}
#'  \item{\strong{c}}{The number of participants completing the assigned intervention.}
#'  \item{\strong{na}}{The number of compared interventions.}
#' }
#' Apart from \strong{na}, all other elements appear in \code{data} as many times as the maximum number of interventions compared in a trial. See, 'Example'.
#'
#' @seealso \code{\link{R2jags}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. Stat Sci. 1992;7:457–472.
#'
#' \dontshow{load("netmodr/data/One-stage model_NMA Dataset.RData")}
#' @examples
#' ### Show the data (one-trial-per-row format)
#' (data <- as.data.frame(one.stage.dataset.NMA[[3]]))
#'
#' ### Run a random-effects network meta-analysis with consistency equations for the standardised mean difference
#' ### assuming missing at random for identical, common informative missingness difference of means.
#' run.model(data = data, measure = "SMD", assumption = "IDE-COMMON", mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
run.model <- function(data, measure, assumption, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin){

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


  jagsfit <- jags(data = data.jag, parameters.to.save = param.jags, model.file = paste0("./model/Full RE-NMA/Full RE-NMA_", measure, "_Pattern-mixture_", assumption, ".txt"),
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

return(list(EM = EM, tausq = tausq, SUCRA = SUCRA, effectiveness = effectiveness, phi = phi))
}

