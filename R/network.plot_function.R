##############################################################################################
#                                                                                            #
#                                                                                            #
#             Create network plots that illustrate MOD extent in nodes and links             #
#                                                                                            # 
#                                                                                            # 
##############################################################################################



network.plot.mod.extension <- function(dataset, names){
  

  ## Load packages
  list.of.packages <- c("dplyr", "gemtc", "BUGSnet")
  lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages) 
  
  
  ## Load functions
  source("./32_Functions/Contrast-level binary data_Function.R")
  source("./32_Functions/Arm-level Long dataset_Arm re-arrangement_Function.R")
  

  ## Obtain dataset
  (r <- dataset %>% dplyr::select(starts_with("r")))
  (m <- dataset %>% dplyr::select(starts_with("m")))
  (n <- dataset %>% dplyr::select(starts_with("n")))
  (t <- dataset %>% dplyr::select(starts_with("t")))
  (nt <- length(table(as.matrix(t))))
  (ns <- length(r[, 1]))

  
  ## Turn arm-level to contrast-level dataset
  pairwise <- contrast.level.dataset(as.matrix(r), as.matrix(m), as.matrix(n), as.matrix(t))

  
  ## Convert one-study-per-row data to one-arm-per-row as required in GeMTC
  data.arranged.arms <- as.data.frame(arm.level.long.dataset.arm.arrange(pairwise))
  (r.new <- data.arranged.arms %>% dplyr::select(starts_with("r")))
  (n.new <- data.arranged.arms %>% dplyr::select(starts_with("n") & !ends_with("a")))
  (t.new <- data.arranged.arms %>% dplyr::select(starts_with("t")))
  na..  <- rep(0, length(r.new[, 1]))
  for(i in 1:length(r.new[, 1])){
    na..[i] <- table(!is.na(t.new[i, ]))["TRUE"] 
  }
  
  
  ## Rename columns to agree with gemtc  
  names(r.new) <- paste0("r..",1:length(r.new[1, ]),".")
  names(n.new) <- paste0("n..",1:length(n.new[1, ]),".")
  names(t.new) <- paste0("t..",1:length(t.new[1, ]),".")

  
  ## Calculate summary of %MOD in each intervention
  risk.drug <- round(aggregate(unlist(m), by = list(unlist(t)), sum)[, 2]/aggregate(unlist(n), by = list(unlist(t)), sum)[, 2], 2)
  
  
  ## Calculate %total MOD per observed comparison
  trial.mod <- apply(pairwise[, 4:5], 1, sum)   
  trial.size <- apply(pairwise[, 6:7], 1, sum) 
  comp <- paste0(pairwise[, 8], "vs", pairwise[, 9])
  risk.comp <- round(aggregate(trial.mod, by = list(comp), sum)[, 2]/aggregate(trial.size, by = list(comp), sum)[, 2], 2)
  
  
  ## one row per study arm
  transform0 <- mtc.data.studyrow(cbind(t.new , r.new , n.new , na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
  for(i in 1:length(unique(transform0$study))){
    transform0[transform0$treatment == i, 2] <- names[i]
  }
  
  
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



