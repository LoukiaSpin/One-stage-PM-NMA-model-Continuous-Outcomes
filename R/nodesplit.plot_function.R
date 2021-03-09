#' Plot the results from the node-splitting approach
#'
#' @export
nodesplit.plot <- function(net, drug.names) {

  drug.names <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N")
  net <- node1

  ## Keep results on 'direct evidence', 'indirect evidence', 'inconsistency factor', 'between-trial standard deviation',
  ## and model assessment measures (i.e., DIC, posterior mean of refisual deviance, and pD)
  direct <-net$direct; indirect <-net$EM; IF <-net$diff; tau <- net$tau; model.assess <- net$model.assessment


  ## Interventions' name: Replace code with original names
  # For treat1 (non-baseline arm)
  for(i in sort(unique(unlist(direct[, 1])))) {
    direct[direct$treat1 == i, 1] <- drug.names[i]
  }

  # For treat2 (baseline arm)
  for(i in sort(unique(unlist(direct[, 2])))) {
    direct[direct$treat2 == i, 2] <- drug.names[i]
  }


  ## Prepare the dataset to create the panel of forest-plot on the 'direct evidence', 'indirect evidence', and 'inconsistency factor' for each split node
  comp <- paste(direct[, 1], "vs", direct[, 2])
  prepare <- data.frame(rep(comp, 3), rbind(direct[, c(3, 5:6)], indirect[, c(3, 5:6)], IF[, c(3, 5:6)]), rep(c("direct", "indirect", "IF"), each = length(direct[, 1])))
  colnames(prepare) <- c("node", "mean", "lower", "upper", "evidence")
  prepare$stat.signif <- ifelse(prepare$lower > 0 | prepare$upper < 0  , "statistically significant", "statistically non-significant")
  prepare$stat.signif <- ifelse(prepare$evidence != "IF", NA, prepare$stat.signif)



  ## Create the panel of forest-plots
  if(length(unique(comp)) <= 16) {

    ggplot(data = prepare[prepare$node == "10vs1", ], aes(x = factor(evidence, levels = c("IF", "indirect", "direct")), y = mean, ymin = lower, ymax = upper, colour = stat.signif) ) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, lty = 2, size = 1.5, col = "grey") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(evidence), y = round(mean, 2), label = round(mean, 2)), color = "black", hjust = -0.2, vjust = -0.3, size = 4.0,
                check_overlap = F, parse = F, position = position_dodge(width = 0.8),  inherit.aes = T) +
      facet_grid(. ~ node) +
      coord_flip() +
      labs(x = "", y = "", colour = "") +
      scale_color_manual(breaks = c("statistically significant", "statistically non-significant"), values = c("green4", "red"), na.value = "black") +
      theme_classic() +
      theme(axis.title.y = element_text(color = "black", size = 14, face = "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
            axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14), legend.position = "none")


  } else {

    #filter(prepare, node %in% unique(comp)[1:16])
 p1 <- lapply(1:length(direct[, 1]), function(.x)
   ggplot(data = prepare, aes(x = factor(evidence, levels = c("IF", "indirect", "direct")), y = mean, ymin = lower, ymax = upper, colour = stat.signif) ) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, lty = 2, size = 1.5, col = "grey") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(evidence), y = round(mean, 2), label = round(mean, 2)), color = "black", hjust = -0.2, vjust = -0.3, size = 4.0,
                check_overlap = F, parse = F, position = position_dodge(width = 0.8),  inherit.aes = T) +
      facet_wrap(vars(node), scales = "free_x") +
      coord_flip() +
      labs(x = "", y = "", colour = "") +
      scale_color_manual(breaks = c("statistically significant", "statistically non-significant"), values = c("green4", "red"), na.value = "black") +
      theme_classic() +
      theme(axis.title.y = element_text(color = "black", size = 14, face = "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
            axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14), legend.position = "none"))

     marrangeGrob(p1, nrow = 4, ncol = 4)
  }





}
