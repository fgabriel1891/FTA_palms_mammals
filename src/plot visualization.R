# visualization.R
# This script is responsible for creating visualizations of the analysis results.

# Load required libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(sf)

# Load pruned assemblages data
all_assemblages_prunned <- readRDS('00_Data/02_species_interactions/Metaweb.RDS')

# Plot the output of the latent network models
plot_latent_network_models <- function(latent_network_models, nRLQ) {
  par(mfrow = c(2, 2))
  
  heatmap(latent_network_models$SBM_ProbsMat, main = 'Stochastic Block Model')
  heatmap(latent_network_models$B_ProbsMat, main = 'Matching-centrality Model')
  heatmap(latent_network_models$C_ProbsMatrix, main = 'Connectance Model')
  heatmap(latent_network_models$M_ProbsMatrix, main = 'Trait Matching Model')
  heatmap(nRLQ, main = 'Observed', col = c('white', 'black'), Rowv = NA, Colv = NA)
}

# Plot the ROC curves for all models
plot_roc_curves <- function(roc_data_combined, title) {
  ggplot(roc_data_combined, aes(x = 1 - specificity, y = sensitivity, color = Model)) +
    geom_smooth(aes(fill = Model), alpha = 0.2, size = 2) +
    geom_abline(aes(intercept = 0, slope = 1), size = 3) +
    labs(x = "1 - Specificity", y = "Sensitivity", title = title) +
    theme_minimal()
}

# Plot variable importance
plot_variable_importance <- function(var_imp, title) {
  barplot(var_imp, horiz = TRUE, las = 1, main = title)
}

# Plot species groupings
plot_species_groupings <- function(PalmNet, MammNet, SBMs) {
  palm_groupings_p <- PalmNet %>% 
    ggplot(aes(factor(SBMs.SB_H), AverageFruitLength_cm)) + 
    geom_boxplot(aes(fill = factor(SBMs.SB_H)), col = 'black', alpha = 0.6) +
    theme_minimal() + 
    xlab('') + 
    ylab('Palm Fruit Length (log)') +
    theme(legend.position = "none")
  
  mamma_groupings_p <- MammNet %>% 
    ggplot(aes(factor(SBMs.SB_W), BodyMass.Value)) + 
    geom_boxplot(aes(fill = factor(SBMs.SB_W)), col = 'black', alpha = 0.6) + 
    theme_minimal() + 
    xlab('SBM group') + 
    ylab('Mammal body mass (log)') +
    theme(legend.position = "none")
  
  grid_sbm <- SBMs$SBM1$Omega_rs %>% 
    reshape2::melt() %>% 
    ggplot() + 
    geom_tile(aes(Var1, Var2, fill = value), col = 'black', size = 1) + 
    theme_minimal() + 
    xlab('SBM group') + 
    ylab('SBM group') + 
    scale_fill_gradient(low = 'white', high = 'firebrick') + 
    theme(legend.position = "none")
  
  panel1 <- gridExtra::grid.arrange(palm_groupings_p, mamma_groupings_p)
  full_panel <- gridExtra::grid.arrange(panel1, grid_sbm, ncol = 2)
  full_panel
}

# Plot Youden's J statistic comparison
plot_youden_j <- function(YJtestin) {
  ggplot(YJtestin, aes(1 - speci, sens, color = factor(id))) +
    geom_point(size = 3, alpha = 0.4) + 
    geom_line(size = 2, alpha = 0.5) + 
    theme_bw() + 
    geom_abline(aes(slope = 1, intercept = 0), size = 1.5) + 
    xlab('1 - Specificity') + 
    ylab('Sensitivity') + 
    theme_classic()
}

# Example usage:
# plot_latent_network_models(latent_network_models, nRLQ)
# plot_roc_curves(roc_data_combined, "ROC Curves for All Models")
# plot_variable_importance(var_imp_mam, "Mammal Variable Importance")
# plot_variable_importance(var_im_palm, "Palm Variable Importance")
# plot_species_groupings(PalmNet, MammNet, SBMs)
# plot_youden_j(YJtestin)
