#functions for plotting beta parameters from HMSC models





plotBeta2 = function(data, post, supportLevel = 0.9, covlabels = data$covNames) {
  mbeta = post$mean
  betaP = post$support
  toPlot = sign(mbeta)
  toPlot = toPlot * ((betaP > supportLevel) + (betaP < 
                                                 (1 - supportLevel)) > 0)
  toPlot = as.data.frame(toPlot)
  toPlot$covnames = factor(data$covNames, levels = data$covNames,
                           labels = covlabels)
  
  toPlot2 = pivot_longer(toPlot, cols = 1:last_col(1), names_to = "Species", values_to = "support") %>%
    mutate(support = as.factor(support), Species2 = factor(Species, levels = Specieslevels, labels = Specieslabels)) 
  ggplot(toPlot2) + geom_tile(aes(x = covnames, y = Species2, fill = support), color = "grey") +
    scale_fill_manual(values = c("blue", "white", "red"), 
                      labels = c("Negative", "Neutral","Positive"),
                      name = "Relationship")+
    scale_x_discrete(name = NULL) + 
    theme(legend.position = "bottom")+
    scale_y_discrete(name = NULL)
  
}


plotBeta3 = function(data, post, supportLevel, covlabels = data$covNames) {
  mbeta = post$mean
  betaP = post$support
  toPlot = mbeta * ((betaP > supportLevel) + (betaP < 
                                                 (1 - supportLevel)) > 0)
  toPlot = as.data.frame(toPlot)
  toPlot$covnames = factor(data$covNames, levels = data$covNames,
                           labels = covlabels)
  
  toPlot2 = pivot_longer(toPlot, cols = 1:last_col(1), names_to = "Species", values_to = "support") %>%
    mutate( Species2 = factor(Species, levels = Specieslevels, labels = Specieslabels)) 
  ggplot(toPlot2) + geom_tile(aes(x = covnames, y = Species2, fill = support), color = "gray") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
    scale_x_discrete(name = NULL)
  
}


