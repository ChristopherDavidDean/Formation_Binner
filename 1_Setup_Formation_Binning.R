#===============================================================================================================================================
#============================================== FORMATION BINNER - DINOSAURS OF NORTH AMERICA ==================================================
#===============================================================================================================================================

# Christopher D. Dean & Susannah Maidment
# 2019

#=========================================================== DATA SETUP ======================================================================

# Set working directory
setwd("C:/Users/deancd/Documents/RESEARCH/PROJECTS/FRM_BIN/Update_R/") # Set your working directory

# Load in Functions
source("0_Functions_Form_Binner.R") # Import functions from other R file (must be in same working directory)

# Make vector of package names
packages <- c("pbmcapply", "dpylr", "tidyverse", "divDyn", "rowr", "matrixStats", "bleepr")

# Install packages
ipak(packages)

library(pbmcapply)
library(dplyr)
library(tidyverse)
library(divDyn)
library(rowr)
library(matrixStats)
library(beepr)

# Data input
formations <- read.csv (file = "Data/Formations_test2.csv")  #Read in formations
occs <- read.csv(file = "Data/NADINOS-occs-edit.csv") # Read in occurrences

# Standard Bin setup
data(stages)

#==== Data Cleaning ====

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

#=============================================== RUNNING TESTS ===============================================

Quorum <- c(0.4, 0.6, 0.8)
bin_limits <- c(3, max(formations$max_age), 66) # Set user defined bin size
Scoring_Grid_1(formations) # Generates scoring grid. Currently set to default resolution (0.01 Ma intervals)
Scoring_Grid_2(formations)
newBins(score_grid, formations, bin_limits, allbins, stages) # Uses the scoring grid to generate new bins.

FormationGraph(formations, form_bins, stages) # Visualises the range of formations in comparison with stage level bins and new bins.

### Running Methods ###
FormBin_M1(formations, binlist, Form_list, Quorum) # Generates formation binned plots of diversity, sampling proxies and SQS results using an inclusive model
FormBin_M2(formations, binlist, Form_list, Quorum) # Generates formation binned plots of diversity, sampling proxies and SQS results using an exclusive model
FormBin_M3(formations, binlist, Form_list, times = 2, Quorum) # Generates plots of diversity, sampling proxies and SQS results using a representative model. 
# Must be run 2 times or more. At high times of times, might take a while!

# Testing Resolution

Scoring_Grid_1(formations, res=0.05)
newBins(score_grid, formations, bin_limits, allbins, stages)

Scoring_Grid_1(formations, res=0.1)
newBins(score_grid, formations, bin_limits, allbins, stages)

Scoring_Grid_1(formations, res=0.5)
newBins(score_grid, formations, bin_limits, allbins, stages)

Scoring_Grid_1(formations, res=1)
newBins(score_grid, formations, bin_limits, allbins, stages)

# Comparing Scoring Grid Methods
Scoring_Grid_1(formations, res=0.5)
newBins(score_grid, formations, bin_limits, bins, stages)
FormationGraph(formations, newbins_test, stages)

Scoring_Grid2(formations, res=0.5)
newBins(score_grid, formations, bin_limits, bins, stages)
FormationGraph(formations, newbins_test, stages)
