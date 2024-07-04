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

# Load prunned assemblages data
all_assemblages_prunned <- readRDS('00_Data/02_species_interactions/Metaweb.RDS')

# Filter for specific criteria and extract IDs
idtst <- all_assemblages_prunned %>%
  filter(SBM_G == '5', taxa == 'palm') %>%
  pull(id)

# Calculate mean fruit length for specific palm traits
mean(palm_traits %>%
  filter(SpecName %in% idtst) %>%
  pull(AverageFruitLength_cm))

# Function to calculate network metrics for a given grid
calc_net_metric <- function(grid_test, SBMs) {
  # Distinct species by taxa
  a <- grid_test %>% 
    distinct(id, taxa, SBM_G) %>% 
    split(.$taxa)
  
  # Frequency of species within SBM groups
  fr_palm <- table(a[[1]]$SBM_G)
  fr_mammals <- table(a[[2]]$SBM_G)
  
  # Normalize frequency
  fr_norm_palm <- fr_palm / sum(fr_palm)
  fr_norm_mammals <- fr_mammals / sum(fr_mammals)
  
  # Compute FTA
  fta <- abs(fr_norm_palm - fr_norm_mammals)
  
  # Create combinations of species pairs
  n <- expand.grid(pluck(a, 'mammals', 'id'), pluck(a, 'palm', 'id')) 
  
  # Compute area sums for each species
  areas <- grid_test %>% 
    split(.$taxa) %>% 
    map(~ .x %>% 
      mutate(area = as.numeric(area)) %>% 
      group_by(id) %>% 
      summarize(area_sum = sum(area))) %>% 
    bind_rows()
  
  # Join data to combinations
  n <- n %>% 
    left_join(pluck(a, 'mammals'), by = c('Var1' = 'id')) %>% 
    left_join(pluck(a, 'palm'), by = c('Var2' = 'id')) %>% 
    left_join(areas, by = c('Var1' = 'id')) %>% 
    left_join(areas, by = c('Var2' = 'id'))
  
  # Compute interaction probabilities
  n$intPro <- sapply(1:length(n$Var1), function(i) 
    (SBMs$SBM1$Omega_rs[n$SBM_G.x[i], n$SBM_G.y[i]]))
  
  # Compute geographic distance and final interaction
  n <- n %>% 
    mutate(int_area = ((area_sum.x) / sum(area_sum.x, area_sum.y)) * ((area_sum.y) / sum(area_sum.x, area_sum.y)))
  n$n_geog_dist <- st_distance(x = n$geometry.x, y = n$geometry.y) %>% diag()
  n$int_final <- scales::rescale(n$int_area, c(0, 1)) * (scales::rescale(n$int_area, c(0, 1)) * scales::rescale(as.numeric(n$n_geog_dist), c(0, 1)))
  
  # Create interaction matrix and compute specialization metric
  netT <- xtabs(int_final ~ Var1 + Var2, n)
  h2 <- cassandRa::RarefyNetwork(netT, abs_sample_levels = 100, metrics = c("H2"))$H2 %>% median()
  
  return(list(fr_palm = fr_palm, fr_mammals = fr_mammals, fr_norm_palm = fr_norm_palm, fr_norm_mammals = fr_norm_mammals, fta = fta, netT = netT, h2 = h2))
}

# Safely wrap the calc_net_metric function
cal_net_metric_safe <- safely(calc_net_metric)

# Test the function on a single grid and measure execution time
system.time({
  test_res <- cal_net_metric_safe(all_assemblages_prunned %>% split(.$grid) %>% pluck(1), SBMs)
})

# Check test results
test_res$result

# Remove geometry and group by SBM group
all_assemblages_prunned2 <- all_assemblages_prunned
all_assemblages_prunned2$geometry <- NULL
all_assemblages_prunned2 %>% group_by(SBM_G)

# Count unique grids
length(unique(all_assemblages_prunned$grid))

# Load required libraries for parallel processing
library(furrr)
plan(multisession, workers = 10)

# Compute network metrics for each grid in parallel
my_net__output <- all_assemblages_prunned %>% 
  split(.$grid) %>% 
  furrr::future_map(function(grid) {
    cal_net_metric_safe(grid, SBMs)
  })

# Save computed networks to file
saveRDS(my_net__output, '00_Data/02_species_interactions/final-networks-grid.RDS')

# Read saved networks
my_networks <- readRDS('00_Data/02_species_interactions/final-networks-grid.RDS')

# Extract and combine results
fta_obs <- my_networks %>% 
  map(~ .x$result) %>% 
  map(~ .x$fta %>% unlist()) %>% 
  bind_rows() %>% 
  mutate('grid' = names(my_networks))

fr_palms <- my_networks %>% 
  map(~ .x$result) %>% 
  map(~ .x$fr_palm %>% unlist()) %>% 
  bind_rows() %>% 
  mutate('grid' = names(my_networks))

fr_mammals <- my_networks %>% 
  map(~ .x$result) %>% 
  map(~ .x$fr_mammals %>% unlist()) %>% 
  bind_rows() %>% 
  mutate('grid' = names(my_networks))

fr_norm_palms <- my_networks %>% 
  map(~ .x$result) %>% 
  map(~ .x$fr_norm_palm %>% unlist()) %>% 
  bind_rows() %>% 
  mutate('grid' = names(my_networks))

fr_norm_mammals <- my_networks %>% 
  map(~ .x$result) %>% 
  map(~ .x$fr_norm_mammals %>% unlist()) %>% 
  bind_rows() %>% 
  mutate('grid' = names(my_networks))

h2_grid <- my_networks %>% 
  map(~ .x$result) %>% 
  map(~ .x$h2 %>% unlist()) %>% 
  unlist() %>% 
  data.frame() %>% 
  setNames('h2') %>% 
  rownames_to_column('grid')

# Add biogeographic region based on coordinates
xy_sf <- st_as_sf(all_assemblages_prunned, coords = c("cord_x", "cord_y"), crs = st_crs(neotropics))
xy_sf <- st_set_crs(xy_sf, st_crs(neotropics))

all_assemblages_prunned_biog <- st_join(xy_sf, neotropics)
all_assemblages_prunned_biog <- all_assemblages_prunned_biog %>% 
  group_by(id, taxa, grid, Dominions) %>% 
  slice(1)

# Plot species distribution across SBM groups
all_assemblages_prunned_biog2 <- all_assemblages_prunned_biog
all_assemblages_prunned_biog2$geometry <- NULL
sp_per_sbm <- all_assemblages_prunned_biog2 %>%
  group_by(Dominions, SBM_G, taxa) %>%
  summarise(n = n_distinct(id)) %>%
  filter(!is.na(taxa), !is.na(SBM_G), !is.na(Dominions)) %>%
  ggplot() + 
  geom_bar(aes(x = SBM_G, y = n, fill = taxa), stat = "identity") +
  theme_minimal() + 
  facet_wrap(~ Dominions, scales = 'free_y') +
  scale_fill_manual(values = c('firebrick2', 'darkgreen')) 

# Filter for specific SBM group and extract unique IDs
all_assemblages_prunned_biog2 %>%
  filter(SBM_G == 5) %>%
  pull(id) %>%
  unique()

# Define function to compute expected values
get_expected_val <- function(all_assemblages_prunned_biog, grid_to_sample, SBMs) {
  biog_to_sample <- (all_assemblages_prunned_biog$Dominions[all_assemblages_prunned_biog$grid == grid_to_sample] %>% 
                       table() %>% sort(decreasing = TRUE))[1] %>% names()
  
  expected_comm <- all_assemblages_prunned_biog %>%
    filter(Dominions %in% biog_to_sample) %>%
    split(.$taxa) %>%
    map(~ .x %>% 
          group_by(id) %>% 
          slice(1) %>% 
          ungroup() %>% 
          slice_sample(n = 10)) %>%
    bind_rows()
  
  return(calc_net_metric2(expected_comm, SBMs))
}

# Safely wrap the get_expected_val function
safe_expected_values <- function(n_rep, all_assemblages_prunned_biog, grids_to_sample, SBMs) {
  replicate(n_rep, get_expected_val(all_assemblages_prunned_biog, grids_to_sample, SBMs))
}

safe_expected_values <- safely(safe_expected_values)

# Open a

 parallel cluster to run safe_expected_values in parallel
library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(10)
registerDoParallel(cl)

# Export variables and libraries to each cluster
clusterExport(cl, c('all_assemblages_prunned_biog', 'grids_to_sample', 'SBMs', 'calc_net_metric2', 'get_expected_val', 'safe_expected_values'))
clusterEvalQ(cl, {
  library(tidyverse)
  library(sf)
  library(cassandRa)
  library(vegan)
})

# Sample 100 grids and apply the function in parallel using foreach
gsample <- grids_to_sample
my_null_result_full <- foreach(grid = gsample, .packages = c('tidyverse', 'sf', 'cassandRa', 'vegan')) %dopar% {
  safe_expected_values(50, all_assemblages_prunned_biog, grid, SBMs)
}

# Save null networks to file
saveRDS(my_null_result_full, '00_Data/02_species_interactions/null-networks-grid_final_all2.RDS')

# Read null networks from file and extract results
my_null_result_full <- readRDS('00_Data/02_species_interactions/null-networks-grid_final_all2.RDS')
names(my_null_result_full) <- gsample
mres <- keep(my_null_result_full, ~ !is.null(.x$result))
mres <- mres %>% map(~ .x$result)

# Compute z-scores for FTA and H2
fta_expected_mean <- mres %>% 
  map(~ .x['fta', ] %>% bind_rows() %>% colMeans()) %>%
  bind_rows() %>%
  mutate(grid = names(mres))

fta_expected_sd <- mres %>% 
  map(~ .x['fta', ] %>% bind_rows() %>% apply(2, sd)) %>%
  bind_rows() %>%
  mutate(grid = names(mres))

fr_palm_mean <- mres %>% 
  map(~ .x['fr_palm', ] %>% bind_rows() %>% colMeans()) %>%
  bind_rows() %>%
  mutate(grid = names(mres))

fr_mammals_mean <- mres %>% 
  map(~ .x['fr_mammals', ] %>% bind_rows() %>% colMeans()) %>%
  bind_rows() %>%
  mutate(grid = names(mres))

fr_norm_palm_mean <- mres %>% 
  map(~ .x['fr_norm_palm', ] %>% bind_rows() %>% colMeans()) %>%
  bind_rows() %>%
  mutate(grid = names(mres))

fr_norm_mammals_mean <- mres %>% 
  map(~ .x['fr_norm_mammals', ] %>% bind_rows() %>% colMeans()) %>%
  bind_rows() %>%
  mutate(grid = names(mres))

h2_mean <- mres %>% 
  map(~ .x['h2', ] %>% unlist() %>% mean(na.rm = TRUE)) %>%
  unlist() %>%
  data.frame() %>%
  setNames('h2_x') %>%
  mutate(grid = names(mres))

h2_sd <- mres %>% 
  map(~ .x['h2', ] %>% unlist() %>% sd(na.rm = TRUE)) %>%
  unlist() %>%
  data.frame() %>%
  setNames('h2sd') %>%
  mutate(grid = names(mres))

h2_obs <- mres %>% 
  imap(~ .x['h2', ] %>% unlist() %>% data.frame() %>% setNames('h2_obs') %>%
         mutate(grid = .y) %>%
         mutate(rep = 1:50)) %>%
  bind_rows()

# Split original vector for parallel processing
original_vector <- 1:length(mres)
size <- 100
split_vectors <- split(original_vector, gl(ceiling(length(original_vector) / size), size, length(original_vector)))

# Compute full FTA expected values
full_fta_expected <- split_vectors %>%
  map(function(split_vec) {
    split_vec %>%
      map(~ {
        1:50 %>%
          map(~ {
            expand.grid(
              mres[[1]]['fr_norm_palm', .x] %>% bind_rows() %>% as.matrix(),
              mres[[1]]['fr_norm_mammals', .x] %>% bind_rows() %>% as.matrix()
            ) %>%
              mutate(lab = bind_rows(
                replicate(50,
                          expand.grid(matrix(rep(1:7), ncol = 7, byrow = TRUE),
                                      matrix(rep(1:7), ncol = 7, byrow = TRUE)),
                          simplify = FALSE)) %>%
                mutate(label = paste0('p', Var1, 'm', Var2)) %>%
                pull(label)) %>%
              mutate(fta = abs(Var1 - Var2)) %>%
              mutate(grid = grids_to_sample[.x]) %>%
              group_by(lab) %>%
              mutate(h2_obs = mres[[.x]]['h2', ] %>% unlist()) %>%
              ungroup()
          }) %>%
          bind_rows()
      }) %>%
      bind_rows()
  }) %>%
  bind_rows()

# Summarize full FTA expected values
full_fta_expected_summ <- full_fta_expected %>%
  group_by(lab, grid) %>%
  summarise(fta_mean = mean(fta), fta_sd = sd(fta))

full_fta_expected_summ$grid <- as.character(full_fta_expected_summ$grid)

# Join summarized FTA expected values with observed values
full_fta_val <- full_fta_expected_summ %>%
  left_join(full_fta, by = c('lab', 'grid'))

full_fta_val$zscore <- (full_fta_val$fta - full_fta_val$fta_mean) / full_fta_val$fta_sd

full_fta_val <- full_fta_val %>%
  left_join(h2_zscore, 'grid')

full_fta_val %>%
  summarize(mn = mean(fta, na.rm = TRUE),
            sd_fta = sd(fta, na.rm = TRUE),
            mnz = median(zscore, na.rm = TRUE)) %>%
  arrange(desc(mnz))

# Plot interaction probabilities within SBM groups
int_per_smb <- SBMs$SBM1$Omega_rs %>%
  reshape2::melt() %>%
  ggplot() +
  geom_tile(aes(Var1, Var2, fill = value), col = 'black', size = 1) +
  theme_minimal() +
  xlab('SBM group (mammals)') +
  ylab('SBM group (palms)') +
  scale_fill_gradient(low = 'white', high = 'firebrick') +
  theme(legend.position = "none") +
  geom_text(aes(Var1, Var2, label = round(value, 2))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Plot average asymmetry heatmap
average_asymmetry_heatmap <- full_fta_val %>%
  group_by(lab) %>%
  summarize(zscore = median(zscore, na.rm = TRUE)) %>%
  xtabs(zscore ~ lab, .) %>%
  matrix(., nrow = 7, byrow = TRUE) %>%
  as.data.frame() %>%
  mutate(across(everything(), ~case_when(is.infinite(.) ~ 0, TRUE ~ .))) %>%
  as.matrix() %>%
  reshape2::melt() %>%
  ggplot() +
  geom_tile(aes(Var1, Var2, fill = value), col = 'black', size = 1) +
  theme_minimal() +
  xlab('SBM group (mammals)') +
  ylab('SBM group (palms)') +
  scale_fill_gradient(low = 'skyblue', high = '#FF6066') +
  theme(legend.position = "none") +
  geom_text(aes(Var1, Var2, label = round(value, 2))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Assemble and display combined plot
library(gridExtra)
library(grid)
gridExtra::grid.arrange(int_per_smb, average_asymmetry_heatmap, ncol = 2)

# Load climatic variables data
library(sf)
library(raster)
neotropics <- st_read('00_Data/03_Landscape/Morrone_Neotropics/Lowenberg_Neto_2014.shp')
grid <- st_make_grid(neotropics, cellsize = c(1, 1), what = "polygons", crs = sf::st_crs(st_read(palm_shp_files[1])))
grid <- st_sf(grid)
WCLim <- raster::getData("worldclim", var = "bio", res = 10)
cropMask <- function(raster, prov) {
  r2 <- crop(raster, extent(prov))
  r3 <- mask(r2, prov)
  return(r3)
}
WCLim <- cropMask(WCL

im, neotropics)
Temp <- WCLim[[1]]
Prec <- WCLim[[12]]
PrecSe <- WCLim[[15]]
IsoTer <- WCLim[[3]]
TempSeaso <- WCLim[[14]]
Temp <- aggregate(Temp, 1 / 0.17)
Prec <- aggregate(Prec, 1 / 0.17)
PrecSe <- aggregate(PrecSe, 1 / 0.17)
TempSeaso <- aggregate(TempSeaso, 1 / 0.17)

# Extract climatic data for grid
gridTemp <- raster::extract(Temp, Assemblages)
gridPrec <- raster::extract(Prec, Assemblages)
gridTS <- raster::extract(TempSeaso, Assemblages)
gridPS <- raster::extract(PrecSe, Assemblages)

# Function to join climatic variables to biological data
add_clim_data <- function(z_score_table, coordinates_grid) {
  grid_coords <- as.data.frame(coordinates_grid[as.numeric(z_score_table$grid), ])
  z_score_table <- cbind(z_score_table, grid_coords)
  grid_coords <- st_as_sf(grid_coords, coords = c("X", "Y"), crs = st_crs(neotropics))
  clim_var <- data.frame(
    'Temp' = raster::extract(Temp, sf::st_as_sf(grid_coords)),
    'Prec' = raster::extract(Prec, sf::st_as_sf(grid_coords)),
    'TS' = raster::extract(TempSeaso, sf::st_as_sf(grid_coords)),
    'PS' = raster::extract(PrecSe, sf::st_as_sf(grid_coords))
  )
  ggsf <- data.frame(z_score_table, clim_var)
  ggsf <- na.omit(ggsf)
  ggsf <- ggsf %>%
    reshape2::melt(id.vars = c('grid', 'Temp', 'Prec', 'TS', 'PS'), value.name = c('obs_ab'), variable.name = 'SBM_G') %>%
    reshape2::melt(id.vars = c('grid', 'SBM_G', 'obs_ab'), value.name = c('clim_val'), variable.name = 'clim_var') %>%
    filter(!SBM_G %in% c('X', 'Y'))
  return(ggsf)
}

add_clim_data2 <- function(z_score_table, coordinates_grid) {
  grid_coords <- as.data.frame(coordinates_grid[as.numeric(z_score_table$grid), ])
  z_score_table <- cbind(z_score_table, grid_coords)
  grid_coords <- st_as_sf(grid_coords, coords = c("X", "Y"), crs = st_crs(neotropics))
  clim_var <- data.frame(
    'Temp' = raster::extract(Temp, sf::st_as_sf(grid_coords)),
    'Prec' = raster::extract(Prec, sf::st_as_sf(grid_coords)),
    'TS' = raster::extract(TempSeaso, sf::st_as_sf(grid_coords)),
    'PS' = raster::extract(PrecSe, sf::st_as_sf(grid_coords))
  )
  ggsf <- data.frame(z_score_table, clim_var)
  ggsf <- na.omit(ggsf)
  return(ggsf)
}

# Plot palms relationship with climate
fr_palms_with_climate <- add_clim_data(fr_norm_palms, coordinates_grid)
fr_norm_palm_mean_with_climate <- add_clim_data(fr_norm_palm_mean, coordinates_grid)
fr_palms_with_climate %>%
  ggplot(aes(scale(clim_val), obs_ab)) +
  geom_point(alpha = 0.05, col = 'firebrick2', size = 0.3) + 
  facet_wrap(~ clim_var + SBM_G, ncol = 7, nrow = 4, scales = 'free') +
  theme_minimal() +
  labs(title = 'Mammal Functional Richness ~ Climate', 
       y = 'Proportional abundance', x = 'Climate variable') + 
  geom_smooth(method = 'lm', col = 'firebrick2') +
  theme(
    strip.text = element_text(size = 8, face = "plain"),
    strip.background = element_blank()
  ) +
  theme(strip.text = element_text(margin = margin(0, 0, 0, 0))) + 
  geom_smooth(aes(scale(clim_val), obs_ab), 
              method = 'lm',
              col = 'gray',
              data = fr_norm_palm_mean_with_climate)

# Plot mammals relationship with climate
fr_mammals_with_climate <- add_clim_data(fr_norm_mammals, coordinates_grid)
fr_norm_mammals_mean_with_climate <- add_clim_data(fr_norm_mammals_mean, coordinates_grid)
fr_mammals_with_climate %>%
  ggplot(aes(scale(clim_val), obs_ab)) +
  geom_point(alpha = 0.05, col = 'darkgreen', size = 0.3) + 
  facet_wrap(~ clim_var + SBM_G, ncol = 7, nrow = 4, scales = 'free') +
  theme_minimal() +
  labs(title = 'Palm Functional Richness ~ Climate', 
       y = 'Proportional abundance', x = 'Climate variable') + 
  geom_smooth(method = 'lm', col = 'darkgreen') +
  theme(
    strip.text = element_text(size = 8, face = "plain"),
    strip.background = element_blank()
  ) +
  theme(strip.text = element_text(margin = margin(0, 0, 0, 0))) + 
  geom_smooth(aes(scale(clim_val), obs_ab), 
              method = 'lm',
              col = 'gray',
              data = fr_norm_mammals_mean_with_climate)

# Compute the influence of climate on asymmetry
full_fta_val_wt_clim <- add_clim_data2(full_fta_val, coordinates_grid)
full_fta_val_wt_clim <- full_fta_val_wt_clim %>%
  filter(!is.infinite(zscore), !is.na(zscore))
lm_asym_all <- lm(zscore ~ scale(Temp) * lab + scale(Prec) * lab + scale(TS) * lab + scale(PS) * lab, data = full_fta_val_wt_clim)
summary(lm_asym_all)
sjPlot::tab_model(lm_asym_all)

# Define function to make plots as a function of the climate variable term
make_climate_data_plot <- function(model, var_term) {
  predDat <- sjPlot::plot_model(model, type = 'pred', terms = c(var_term, 'lab'), return.data = TRUE)
  plot_labels <- c('Temp' = 'Mean Annual Temperature', 'Prec' = 'Total Annual Precipitation', 'TS' = 'Temperature Seasonality', 'PS' = 'Precipitation Seasonality')
  predicted_data_df <- predDat$data
  predicted_data_df <- as.data.frame(predicted_data_df)
  model_p <- broom::tidy(model)
  p_signif <- model_p$term[model_p$p.value < 0.05] %>%
    data.frame() %>%
    setNames('term') %>%
    filter(str_detect(term, var_term)) %>%
    filter(str_detect(term, "p[0-9]m[0-9]")) %>%
    mutate(term = str_extract(term, 'p[0-9]m[0-9]')) %>%
    pull('term')
  predicted_data_df <- predicted_data_df %>%
    mutate(group = str_extract(group_col, 'p[0-9]')) %>%
    mutate(group2 = str_extract(group_col, 'm[0-9]')) %>%
    mutate(signif = ifelse(group_col %in% p_signif, "p < 0.05", "p > 0.05"))
  signif_levels <- unique(predicted_data_df$signif)
  if (all(length(signif_levels) == 1 & signif_levels == 'p > 0.05')) {
    ggplot(predicted_data_df, aes(x = scale(x), y = predicted, color = group2)) +
      geom_line(size = 1, linetype = 2) +
      labs(title = 'Predicted z-score of FTA',
           x = plot_labels[var_term],
           y = 'FTA z-score',
           color = 'Interaction Terms') +
      theme_minimal() +
      theme(legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8)) +
      facet_wrap(~ group, ncol = 7)
  } else {
    ggplot(predicted_data_df, aes(x = scale(x), y = predicted, color = group2)) +
      geom_line(aes(linetype = signif), size = 1) +
      labs(title = 'Predicted z-score of FTA',
           x = plot_labels[var_term],
           y = 'FTA z-score',
           color = 'Interaction Terms',
           linetype = 'p-values') +
      theme_minimal() +
      theme(legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8)) +
      facet_wrap(~ group, ncol = 7)
  }
}

# Generate climate data plots
T_fta_plot <- make_climate_data_plot(lm_asym_all, 'Temp')
Prec_fta_plot <- make_climate_data_plot(lm_asym_all, 'Prec')
PS

_fta_plot <- make_climate_data_plot(lm_asym_all, 'PS')
TS_fta_plot <- make_climate_data_plot(lm_asym_all, 'TS')
gridExtra::grid.arrange(T_fta_plot, Prec_fta_plot, PS_fta_plot, TS_fta_plot, ncol = 2)

# Plot precipitation seasonality map
library(viridis)
precip_df <- as.data.frame(scale(PrecSe), xy = TRUE)
ps_map <- ggplot() +
  geom_raster(data = precip_df, aes(x = x, y = y, fill = bio15)) +
  scale_fill_viridis(name = "Precipitation\nSeasonality", option = "D", na.value = "transparent") +
  geom_sf(data = neotropics, fill = NA, color = "black", size = 0.5) +
  coord_sf(xlim = c(-120, -30), ylim = c(-56, 33), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
gridExtra::grid.arrange(PS_fta_plot, ps_map, ncol = 2)

# Compute observed relationship between FTA and H2
full_fta <- 1:nrow(fr_norm_palms) %>%
  map(function(row) {
    expand.grid(
      fr_norm_palms[row, ] %>% select(!grid) %>% as.numeric(),
      fr_norm_mammals[row, ] %>% select(!grid) %>% as.numeric()) %>%
      mutate(lab = expand.grid(1:7, 1:7) %>%
               mutate(label = paste0('p', Var1, 'm', Var2)) %>%
               pull(label)) %>%
      mutate(fta = abs(Var1 - Var2)) %>%
      mutate(grid = fr_norm_palms$grid[row])
  }) %>%
  bind_rows()

# Join FTA and H2 grid data
my_fta_h2 <- full_fta %>%
  left_join(h2_grid, 'grid')

# Summarize FTA and H2 data
my_fta_h2_sum <- my_fta_h2 %>%
  group_by(grid) %>%
  summarise(mean_fta = mean(fta, na.rm = TRUE), 
            sd_fta = sd(fta, na.rm = TRUE),
            h2 = mean(h2))

# Fit observed model for FTA and H2
obs_model1 <- lm(h2 ~ mean_fta + sd_fta, data = my_fta_h2_sum)
sjPlot::tab_model(obs_model1)

# Save full FTA expected values to file
saveRDS(full_fta_expected, file = '00_Data/02_species_interactions/full_fta_expected.RDS')

# Compute distribution of expected coefficients between simulated FTA and H2
full_fta_expected <- full_fta_expected %>%
  group_by(grid, lab) %>%
  mutate(rep = rep(1:50))

h2_mod_fta <- full_fta_expected %>%
  group_by(grid, rep) %>%
  summarise(mean_fta = mean(fta, na.rm = TRUE), 
            sd_fta = sd(fta, na.rm = TRUE),
            h2_obs = mean(h2_obs)) %>%
  group_map(~ lm(h2_obs ~ mean_fta + sd_fta, data = .x))

h2_mod_fta_coef <- h2_mod_fta %>%
  map(~ coef(.x)) %>%
  bind_rows()

h2_mod_fta %>%
  map(~ vegan::RsquareAdj(.x)[[1]]) %>%
  unlist()

# Plot distribution of model estimates
vline_positions <- obs_model1 %>%
  coef() %>%
  data.frame()
vline_positions$variable <- rownames(vline_positions)
names(vline_positions) <- c('value', 'variable')

custom_labels <- c("(Intercept)" = "Model Intercept",
                   "mean_fta" = "FTA (mean)", 
                   "sd_fta" = "FTA (sd)")

h2_mod_fta_coef %>%
  reshape2::melt() %>%
  ggplot() +
  geom_histogram(aes(value)) + 
  facet_wrap(~ variable, labeller = as_labeller(custom_labels)) +
  geom_vline(data = vline_positions, aes(xintercept = value),
             color = "red", linetype = 1, size = 1) +
  theme_minimal() +
  ylab('Frequency count') +
  xlab('Model estimate')
