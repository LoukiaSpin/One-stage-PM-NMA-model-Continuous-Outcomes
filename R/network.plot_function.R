#' The network plot
#'
#' @description
#' A function to create the network plot and illustrate the risk of bias due to missing outcome data in each node and link using different colours.
#' The node refers to the intervention and the link refers to the observed pairwise comparison. The risk of bias due to missing outcome data
#' can be low (the proportion of missing outcome data is up to 5%), moderate (the proportion of missing outcome data is more than 5% and up to 20%),
#' or high (the proportion of missing outcome data is more than 20%). Green, orange and red represent low, moderate, and high risk of bias.
#' For each node, the risk of bias is determined by the total proportion of missing outcome data in the corresponding intervention, namely, the ratio of
#' the sum of missing outcome data to the sum of randomised in that intervention.
#' For each link, the risk of bias is determined by the total proportion of missing outcome data across the corresponding trials, namely, the ratio of
#' the sum of missing outcome data to the sum of randomised in these trials.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#' See 'Format' in \code{nma.continuous.full.model} for the specification of the columns.
#' @param drug.names A vector of characteristics with name of the interventions as appear in the function \code{nma.continuous.full.model}.
#' @param outcome Character string indicating the type of outcome with values \code{"binary"}, or \code{"continuous"}.
#'
#' @return A network plot with coloured nodes and links to indicate the risk of bias due to missing outcome data.
#'
#' \dontshow{load("./data/NMA Dataset Continuous.RData")}
#' @examples
#' ### Show the data (one-trial-per-row format)
#' (data <- as.data.frame(one.stage.dataset.NMA[[3]]))
#' drug.names <- sapply(1:8, function(x) letters[x])
#'
#' netplot(data = res, outcome = "binary, drug.names = drug.names)
#'
#' @export
netplot <- function(data, outcome, drug.names){


  if(outcome == "binary") {


    ## Obtain dataset
    (r <- data %>% dplyr::select(starts_with("r")))
    (m <- data %>% dplyr::select(starts_with("m")))
    (n <- data %>% dplyr::select(starts_with("n")))
    (t <- data %>% dplyr::select(starts_with("t")))
    (nt <- length(table(as.matrix(t))))
    (ns <- length(r[, 1]))
    na..  <- rep(0, length(r[, 1]))
    for(i in 1:length(r[, 1])){
      na..[i] <- table(!is.na(t[i, ]))["TRUE"]
    }


    ## Prepare data for gemtc
    data.gemtc <- cbind(t, r, n, na..)


    ## Rename columns to agree with gemtc
    names(r) <- paste0("r..",1:length(r[1, ]),".")
    names(n) <- paste0("n..",1:length(n[1, ]),".")
    names(t) <- paste0("t..",1:length(t[1, ]),".")


    ## one row per study arm
    transform0 <- mtc.data.studyrow(cbind(t, r, n, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
    for(i in 1:length(unique(transform0$study))){
      transform0[transform0$treatment == i, 2] <- drug.names[i]
    }


  } else {

    ## Obtain dataset
    (y <- data %>% dplyr::select(starts_with("y")))
    (sd <- data %>% dplyr::select(starts_with("sd")))
    (c <- data %>% dplyr::select(starts_with("c")))
    (se <- sd/sqrt(c))
    (m <- data %>% dplyr::select(starts_with("m")))
    (n <- data %>% dplyr::select(starts_with("n")))
    (t <- data %>% dplyr::select(starts_with("t")))
    (nt <- length(table(as.matrix(t))))
    (ns <- length(y[, 1]))
    na..  <- rep(0, length(y[, 1]))
    for(i in 1:length(y[, 1])){
      na..[i] <- table(!is.na(t[i, ]))["TRUE"]
    }


    ## Prepare data for gemtc
    data.gemtc <- cbind(t, y, se, n, na..)


    ## Rename columns to agree with gemtc
    names(y) <- paste0("y..",1:length(y[1, ]),".")
    names(se) <- paste0("se..",1:length(se[1, ]),".")
    names(n) <- paste0("n..",1:length(n[1, ]),".")
    names(t) <- paste0("t..",1:length(t[1, ]),".")


    ## one row per study arm
    transform0 <- mtc.data.studyrow(cbind(t, y, se, n, na..), armVars = c('treatment'= 't', 'mean'='y', 'std.err'='se', 'sampleSize'='n'), nArmsVar='na')
    for(i in 1:length(unique(transform0$study))){
      transform0[transform0$treatment == i, 2] <- drug.names[i]
    }

  }


  ## Turn arm-level to contrast-level dataset
  (pairwise <- pairwise(as.list(t), event = as.list(m), n = as.list(n), data = data, studlab = 1:ns)[, c(3:6, 8, 7, 9)])
  colnames(pairwise) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")


  ## Calculate summary of %MOD in each intervention
  risk.drug <- round(aggregate(unlist(m), by = list(unlist(t)), sum)[, 2]/aggregate(unlist(n), by = list(unlist(t)), sum)[, 2], 2)


  ## Calculate %total MOD per observed comparison
  trial.mod <- apply(pairwise[, 4:5], 1, sum)
  trial.size <- apply(pairwise[, 6:7], 1, sum)
  comp <- paste0(pairwise[, 3], "vs", pairwise[, 2])
  risk.comp <- round(aggregate(trial.mod, by = list(comp), sum)[, 2]/aggregate(trial.size, by = list(comp), sum)[, 2], 2)


  ## Replace coded 'treatment' with the original names
  transform <- data.prep(arm.data = transform0, varname.t = "treatment", varname.s = "study")


  ## Obtain the network plot
  net.plot(transform, node.lab.cex = 1.5,
           node.colour = ifelse(risk.drug <= 0.05, "green4", ifelse(risk.drug > 0.20, "brown1", "orange")),
           edge.colour = ifelse(risk.comp <= 0.05, "green4", ifelse(risk.comp > 0.20, "brown1", "orange")))

#  legend("bottom", legend = c("low", "moderate", "high"), horiz = T, seg.len = 0.2, x.intersp = 0.1,
#         inset = c(0, -0.1), text.width = c(0, 0.08, 0.05),
#         bty = "n", xpd = TRUE, cex = 1.2, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("green4", "orange", "brown1"))

}



