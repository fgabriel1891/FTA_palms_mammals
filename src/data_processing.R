# data_processing.R
# This script is responsible for processing the loaded data and preparing it for analysis.

# Load required libraries
library(tidyverse)
library(sf)

# Ensure palm_trait_data and mammal_trait_data are complete cases
palm_trait_data <- palm_traits[complete.cases(palm_traits),]
mammal_trait_data <- mammal_traits[complete.cases(mammal_traits),]

# Filter interaction data to match species between databases
int_data <- int_data %>% 
  filter(PALM %in% palm_trait_data$SpecName,
         FRUGIVORE %in% mammal_trait_data$Scientific)

# Create binary interaction matrix N
N <- int_data %>% 
  xtabs(~PALM + FRUGIVORE, .)
N[N > 1] <- 1

# Filter palm trait data to match interaction data
palm_trait_data <- palm_trait_data %>% 
  filter(SpecName %in% int_data$PALM)

# Select relevant columns for palm traits
pTRLQ <- palm_trait_data %>% 
  column_to_rownames("SpecName") %>% 
  select(Acaulescent, Erect, MaxStemHeight_m, AverageFruitLength_cm)

# Filter mammal trait data to match interaction data
mammal_trait_data <- mammal_trait_data %>% 
  filter(Scientific %in% int_data$FRUGIVORE)

# Select relevant columns for mammal traits
mTRLQ <- mammal_trait_data %>% 
  column_to_rownames("Scientific") %>%
  select(!MSWFamilyLatin)

# Normalize interaction matrix to binary presence-absence
nRLQ <- N[rownames(pTRLQ), rownames(mTRLQ)]
nRLQ <- decostand(nRLQ, "pa")
