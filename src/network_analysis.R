# network_analysis.R
# This script is responsible for performing network modeling and statistical analysis.

# Load required libraries
library(cassandRa)
library(vegan)

# Prepare object to fit latent trait models
Ng <- cassandRa::CreateListObject(nRLQ)

# Fit all models
latent_network_models <- cassandRa::FitAllModels(Ng)

# Alternatively, load pre-fitted models if they exist
# latent_network_models <- readRDS('path_to_saved_models/latent_net_mod.RDS')

# Function to compute Youden's J statistic
TestYJ <- function(probNet, obs, n){
  sq <- seq(range(probNet)[1], range(probNet)[2], diff(range(probNet))/n)
  sens <- c()
  speci <- c()
  YJ <- c()
  
  for (i in 1:n) {
    prob10 <- ifelse(probNet > sq[i], 1, 0)
    Ttab <- prop.table(table(obs, prob10))
    sens[i] <- Ttab[4] / (Ttab[4] + Ttab[2])
    speci[i] <- Ttab[1] / (Ttab[1] + Ttab[3])
    YJ[i] <- sens[i] + speci[i] - 1
  }
  
  data.frame(sens, speci, YJ)
}

# Apply TestYJ function to the models
probNet <- list(latent_network_models$SBM_ProbsMat, 
                latent_network_models$C_ProbsMatrix, 
                latent_network_models$M_ProbsMatrix, 
                latent_network_models$B_ProbsMat)

YJtestin <- lapply(1:4, function(i) TestYJ(probNet[[i]], latent_network_models$obs, 100))

# Rearrange the resulting dataset and name variables appropriately
YJtestin <- YJtestin %>%
  set_names(c('SBM', 'Cent', 'Matc', 'Match_Cent')) %>%
  imap(~{
    .x %>%
      mutate(id = .y)
  }) %>%
  bind_rows()

# Refit SBM model to the continental network data
SBMs <- cassandRa::FitSBM(Ng)

# Calculate summary statistics of the fitted model
mean(diag(SBMs$SBM1$Omega_rs)) / mean(SBMs$SBM1$Omega_rs[upper.tri(SBMs$SBM1$Omega_rs)])

# Extract fitted matrix of species group associations
PalmNet <- data.frame(Ng$HostNames, "SBMs.SB_H" = SBMs$SBM1$SB_H)
MammNet <- data.frame(Ng$WaspNames, "SBMs.SB_W" = SBMs$SBM1$SB_W)

# Join synthetic data with observed trait data
PalmNet <- data.frame(PalmNet, palm_traits[match(PalmNet$Ng.HostNames, palm_traits$SpecName),])
MammNet <- data.frame(MammNet, mammal_traits[match(MammNet$Ng.WaspNames, mammal_traits$Scientific),])

# Fit multinomial logistic regression models
library(tidymodels)

# Prepare data for mammals
MammNet <- MammNet %>% mutate(SBMs.SB_W = as.factor(SBMs.SB_W))
data_split <- initial_split(MammNet, prop = 0.80)
data_train <- training(data_split)
data_test <- testing(data_split)

rec_list <- list(
  "recipe1" = recipe(SBMs.SB_W ~ BodyMass.Value + Diet.Fruit, data = MammNet) %>%
    step_log(BodyMass.Value, base = 10),
  "recipe2" = recipe(SBMs.SB_W ~ BodyMass.Value + Diet.Fruit, data = MammNet),
  "recipe3" = recipe(SBMs.SB_W ~ BodyMass.Value + Diet.Fruit, data = MammNet) %>%
    step_log(BodyMass.Value, base = 10),
  "recipe4" = recipe(SBMs.SB_W ~ BodyMass.Value + Diet.Fruit, data = MammNet)
)

model_spec <- multinom_reg() %>%
  set_engine("nnet") %>%
  set_mode("classification")

workflows <- map(rec_list, ~workflow() %>%
                   add_recipe(.x) %>%
                   add_model(model_spec) %>%
                   fit(data_train))

# Calculate AUC for each model
aucs_mam <- map_df(workflows, ~.x %>%
                     augment(data_test) %>%
                     roc_auc(truth = SBMs.SB_W, .pred_1:.pred_7))

# Plot ROC curves
roc_data <- map(workflows, ~.x %>%
                  augment(data_test) %>%
                  roc_curve(truth = SBMs.SB_W, .pred_1:.pred_7))

roc_data_combined <- bind_rows(roc_data, .id = "Model")

roc_mammals <- roc_data_combined %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = Model)) +
  geom_smooth(aes(fill = Model), alpha = 0.2, size = 2) +
  geom_abline(aes(intercept = 0, slope = 1), size = 3) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves for All Models") +
  theme_minimal()

# Prepare data for palms
PalmNet <- PalmNet %>% mutate(SBMs.SB_H = as.factor(SBMs.SB_H))
data_split <- initial_split(PalmNet, prop = 0.80)
data_train <- training(data_split)
data_test <- testing(data_split)

rec_list <- list(
  "recipe1" = recipe(SBMs.SB_H ~ MaxStemHeight_m + AverageFruitLength_cm, data = data_train) %>%
    step_log(AverageFruitLength_cm, base = 10, offset = 1),
  "recipe2" = recipe(SBMs.SB_H ~ MaxStemHeight_m + AverageFruitLength_cm, data = data_train),
  "recipe3" = recipe(SBMs.SB_H ~ MaxStemHeight_m, data = data_train) %>%
    step_log(MaxStemHeight_m, base = 10, offset = 1),
  "recipe4" = recipe(SBMs.SB_H ~ AverageFruitLength_cm, data = data_train) %>%
    step_log(AverageFruitLength_cm, base = 10, offset = 1),
  "recipe5" = recipe(SBMs.SB_H ~ AverageFruitLength_cm, data = data_train)
)

model_spec <- multinom_reg() %>%
  set_engine("nnet") %>%
  set_mode("classification")

workflows <- map(rec_list, ~workflow() %>%
                   add_recipe(.x) %>%
                   add_model(model_spec) %>%
                   fit(data_train))

# Calculate AUC for each model
aucs <- map_df(workflows, ~.x %>%
                 augment(data_test) %>%
                 roc_auc(truth = SBMs.SB_H, .pred_1:.pred_7))

# Plot ROC curves
roc_data <- map(workflows, ~.x %>%
                  augment(data_test) %>%
                  roc_curve(truth = SBMs.SB_H, .pred_1:.pred_7))

roc_data_combined <- bind_rows(roc_data, .id = "Model")

roc_palm_plot <- roc_data_combined %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = Model)) +
  geom_point() +
  geom_smooth() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves for All Models")

# Refit using nnet and calculate variable importance for mammals
refit_mammal <- nnet::multinom(SBMs.SB_W ~ log(BodyMass.Value) + Diet.Fruit, data = MammNet)
var_imp_mam <- t((caret::varImp(refit_mammal)))
colnames(var_imp_mam) <- str_remove(colnames(var_imp_mam), 'MSWFamilyLatin')
colnames(var_imp_mam) <- str_remove(colnames(var_imp_mam), '.Value')
sjPlot::tab_model(refit_mammal)

# Plot variable importance for mammals
var_imp_plot_mam <- barplot(var_imp_mam, horiz = TRUE, las = 1)

# Refit using nnet and calculate variable importance for palms
refit_palm <- nnet::multinom(SBMs.SB_H ~ MaxStemHeight_m + AverageFruitLength_cm + Acaulescent + Erect, data = PalmNet)
var_im_palm <- t((caret::varImp(refit_palm)))
colnames(var_im_palm) <- str_remove(colnames(var_im_palm), 'PalmTribe')
sjPlot::tab_model(refit_palm)

# Plot variable importance for palms
var_imp_plot_palm <- barplot(var_im_palm, horiz = TRUE, las = 1)

# Predict species group associations and make assemblages
PalmPreds <- data.frame("spNamePalm" = palm_traits$SpecName, "group" = predict(refit_palm, palm_traits, allow.new.levels = TRUE))
mammPreds <- data.frame("spNameMam" = MammNet$Scientific, "group" = predict(refit_mammal, MammNet))
palm_grids <- palm_grids %>% set_names(str_replace(str_remove(basename(palm_shp_files), '.shp'), '_', " "))
palm_grids <- keep(palm_grids, ~ !is.null(.x$result))

# Obtain centroids for grid features
sf::sf_use_s2(FALSE)
centroids_mammals <- mammal_grids %>% imap(~st_centroid(.x) %>% st_coordinates() %>% data.frame() %>% mutate(id = .y, area = st_area(.x))) %>% bind_rows()
get_palm_centroids <- function(palm_grids) {
  palm_grids %>% imap(~st_centroid(.x$result) %>% st_coordinates() %>% data.frame() %>% mutate(id = .y, area = st_area(.x$result))) %>% bind_rows()
}
safe_get_palm_centroids <- safely(get_palm_centroids)
centroids_palms <- safe_get_palm_centroids(palm_grids)
centroids_palms$result %>% head()

# Group by matching coordinates and make assemblages for species co-occurring in the same grid cell
all_assemblages <- centroids_mammals %>% rbind(centroids_palms$result %>% dplyr::select(!X2))
all_assemblages <- all_assemblages %>% mutate(taxa = case_when(id %in% palm_traits$SpecName ~ 'palm', id %in% mammal_traits$Scientific ~ 'mammals', TRUE ~ NA_character_), grid_id = paste0(X, '_', Y))
all_assemblages <- st_as_sf(all_assemblages, coords = c('X', 'Y'), crs = st_crs(grid))
all_assemblages <- st_set_crs(all_assemblages, value = st_crs(grid))
int <- st_intersects(all_assemblages$geometry, grid)
all_assemblages$grid <- unlist(int)

# Join trait data and filter grids with at least 5 species
all_preds_sbm <- rbind(PalmPreds %>% setNames(c('id', 'SBM_G')), mammPreds %>% setNames(c('id', 'SBM_G')))
all_assemblages <- all_assemblages %>% left_join(all_preds_sbm, c('id'))
table_taxa_grid <- all_assemblages %>% split(.$grid) %>% imap(~(table(.$taxa)) %>% data.frame() %>% mutate(id = .y))
table_taxa_grid <- table_taxa_grid %>% bind_rows()
richtab <- xtabs(Freq ~ id + Var1, table_taxa_grid) %>% matrix(., ncol = 2) %>% rownames()
all_assemblages_prunned <- all_assemblages %>% filter(grid %in% richtab)
saveRDS(all_assemblages_prunned, '00_Data/02_species_interactions/Metaweb.RDS')
