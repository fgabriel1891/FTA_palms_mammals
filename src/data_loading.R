# Distinct functional responses of consumers and their producers to climate drive mutualistic network asymmetry
# Gabriel Munoz et al., in progress. 
# June 2024
# data_loading.R
# This script is responsible for loading the necessary data files. Paths might change and need to be adapted. 

# Load required libraries
library(tidyverse) # for data manipulation 
library(sf) # for simple feature manipulation 
library(vegan) # ecological analyses 
library(cassandRa) # latent network analyses 

# Load palm distribution shapefiles 
palm_all_files <- list.files("00_Data/00_species_distribution/Palm-distribution-ranges/Shapefiles/", full.names = TRUE)
palm_shp_files <- palm_all_files[str_detect(palm_all_files, ".shp")]
palm_shp_files <- palm_shp_files[!str_detect(palm_shp_files, ".xml")]

# Load Neotropics map shapefile
neotropics <- st_read('00_Data/03_Landscape/Morrone_Neotropics/Lowenberg_Neto_2014.shp')

# Create a grid over the Neotropics map
grid <- st_make_grid(neotropics, cellsize = c(1, 1), what = "polygons", 
                     crs = sf::st_crs(st_read(palm_shp_files[1])))
grid <- st_sf(grid) # Convert the grid to a simple feature collection

# Load gridded species data
palm_grids <- readRDS("00_Data/00_species_distribution/gridded_palm_data.RDS")
mammal_grids <- readRDS("00_Data/00_species_distribution/gridded_mammal_data.RDS")

# Load trait data for palms and mammals
palm_traits <- read.csv('00_Data/01_species_traits/final_palm_trait.csv')
mammal_traits <- read.csv('00_Data/01_species_traits/final_mammal_trait.csv')

# Load interaction data
int_data <- readRDS('00_Data/02_species_interactions/final_int_data.RDS')
