#==================================================================================================================================
#======================================= 1. FORMATION BINNER - DINOSAURS OF NORTH AMERICA =========================================
#==================================================================================================================================

# Christopher D. Dean, A. Alessandro Chiarenza & Susannah Maidment
# 2019

#================================================ DATA SETUP ======================================================================

# Set working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/FRM_BIN/Formation_Binner/Formation_Binner/") # Set your working directory

# Load in Functions
source("0_Functions_Form_Binner.R") # Import functions from other R file (must be in same working directory)

# Make vector of package names
packages <- c("pbapply", "dpylr", "tidyverse", "divDyn", "rowr", "matrixStats", "bleepr", "iNEXT", "reshape2", "RColorBrewer")

# Install packages
ipak(packages)

library(pbapply)
library(dplyr)
library(tidyverse)
library(divDyn)
library(rowr)
library(matrixStats)
library(beepr)
library(reshape2)
library(colorspace)

# Data input
formations <- read.csv (file = "Data/Formations_Final.csv")  #Read in formations
occs <- read.csv(file = "Data/Occurrences_Final.csv") # Read in occurrences

# Standard Bin setup - trim to fit relevant time frame. 
data(stages)
stages <- stages[75:81,] # Set stages to range from Albian to Maastrichtian

#===== Data Cleaning =====

# Make parts Numeric
formations$max_age <- as.numeric(as.character(formations$max_age)) # Make Numeric
formations$min_age <- as.numeric(as.character(formations$min_age)) # Make Numeric 

# Select appropriate formations and order
formations <- formations[which(formations$Location=='WI'),] # Only formations from Western Interior
myformations <- sort(as.vector(formations$Formation)) # Organise

# Select appropriate occurrences
testoccs <- occs[occs$formation %in% myformations,] # Only include occurrences from formation list
testoccs <- droplevels.data.frame(testoccs) # Remove old levels

# Create Formation/occurrences list
Form_list <- split(testoccs, testoccs$formation) # Makes inputted occ data into lists from Formations

# Reorganise formations
formations <- formations[order(formations$Formation),] # Reorganise formations
formations$forbinning <- 1:nrow(formations) # Number formations for easy plotting later

# Add mean occurrence latitude, raw diversity, umber of occurrences and length per formation
formations <- occs %>% 
  rename(Formation = formation) %>%
  group_by(Formation) %>%
  dplyr::summarize(Mean_Lat = mean(latdec, na.rm = TRUE),
                   Diversity = n_distinct(occurrence.genus_name),
                   Occurrences = n()) %>%
  inner_join(formations, by = "Formation")
formations$Range <- formations$max_age - formations$min_age

#=============================================== RUNNING TESTS ===================================================================

#===== Set up =====
Quorum <- c(0.4, 0.6, 0.8) # Sets quorums for running SQS
bin_limits <- c(4, max(formations$max_age), 66) # Set user defined bin size - change the first number to vary resolution in graphs.

#===== Bin generation and comparison =====
Scoring_Grid_1(formations) # Generates scoring grid. Currently set to default resolution (0.01 Ma intervals). Choose either Score_Grid_1 or 2 (find out more in Functions File)
Scoring_Grid_2(formations)
newBins(score_grid, formations, bin_limits, allbins, stages, smallamalg = TRUE) # Uses the scoring grid to generate new bins.
overlap_counter(score_grid) # Generates graph of the number of formations through time
FormationGraph(formations, form_bins, stages, score_grid_2 = TRUE, 
               draw_by = "Lat", Col = "Diversity", legend = TRUE, STAGE = FALSE) # Visualises the range of formations in comparison with stage level bins and new bins.

#===== Running diversity Methods =====
FormBin_M1(formations, binlist, Form_list, Quorum) # Generates formation binned plots of diversity, sampling proxies and SQS results using an inclusive model
FormBin_M2(formations, binlist, Form_list, Quorum) # Generates formation binned plots of diversity, sampling proxies and SQS results using an exclusive model
FormBin_M3(formations, binlist, Form_list, times = 100, Quorum, run_SQS = TRUE) # Generates plots of diversity, sampling proxies and SQS results using a representative model. 
# Must be run 2 times or more. At high times of times, might take a while!

