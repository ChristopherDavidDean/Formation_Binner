# FORMATION BINNER 2
# Written by Christopher Dean, 2016-2018
# 19/11/18
# General Public License 2
# R script for binning dinosaur occurrences by formation groups
#
# Data Required:
# Formations - file containing a list of all formations and their maximum and minimum ages. 
#              Required column headers: Formation, max_age, min_age
# Occurrences - file containing all fossil occurrences with the formation they were found in.
#              Required column headers: Taxa, Formation
# 
# NOTE: All formations in these files MUST be included in both files, and all must appear in alphabetical order.
#
# Functions:
# Scoring_Grid - function to create a score grid to identify best position to draw boundaries between formations.
#                The function looks through time at intervals specified by the user (res - automatically set at 0.001) and gives a suitability
#                score for each formation to draw a boundary at that point. If a formation does not cross that
#                boundary, the formation is given a score of 100. If it does cross that boundary, it works out 
#                how much of the formation sits each side of the boundary, and downweights accordingly - e.g. if 
#                a formation was found to have 10% of it's total range one side of a boundary and 90% the other,
#                it would receive a higher score than a formation where 50% of it's range sat either side. Mean
#                scores are generated for each time bin. 
# Bins -         Looks at the previously generated score_grid and graphically highlights all places where it would
#                be suitable to draw a boundary, above a user specified threshold (thresh). Boundaries are outputted 
#                as a list (binlist).
# FormBin_M1 -   Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and 
#                formations to all bins that they occur in. Produces graphs showing raw dinosaur occurrences,
#                rock outcrop area, and DPSK through time. 
# FormBin_M2 -   Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and 
#                formations to bins that the majority of the formation occurs in. Formations which have lengths 
#                greater than 2x the length of a bin are ignored. Produces graphs showing raw dinosaur occurrences, 
#                rock outcrop area, and DPSK through time.
# FormBin_M3 -   Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and 
#                formations to all bins that they occur in, based on the percentage of the formation that sits
#                within that bin. Dinosaur occs are selected at random from the formation list and not replaced.
#                The test is repeated according to user specified number of runs (times). Produces graphs showing
#                raw dinosaur occurrences, rock outcrop area, and DPSK through time. 

# Ideas
# One suggestion - user specified approx bin length. Choose a starting point that's
# related to the base of the bin length, draw up bins, then find highest score within
# that bin. Draw boundary there. Repeat this process for different starting point. 
# See if has influence on where bins are drawn.
#
#

# Necessary Packages #
install.packages("iNEXT")
install.packages("tidyverse")
install.packages("pbmcapply")
install.packages("divDyn")

library(pbmcapply)
library(iNEXT)
library(dplyr)
library(tidyverse)
library(divDyn)

select <- dplyr::select # ensure the select function is coming from tidyverse package

# Data input #

setwd(dir = "/Users/ChristopherDean/MEGA/PhD/PROJECTS/FRM_BIN/Data/")

formations <- read.csv (file = "Formations_test2.csv")  #Read in formations
occs <- read.csv(file = "NADINOS-occs-edit.csv") # Read in occurrences

formations$max_age <- as.numeric(as.character(formations$max_age))
formations$min_age <- as.numeric(as.character(formations$min_age))
formations <- formations[which(formations$Location=='WI'),]

myformations <- sort(as.vector(formations$Formation))

testoccs <- occs[occs$formation %in% myformations,]
testoccs <- droplevels.data.frame(testoccs)
Form_list <- split(testoccs, testoccs$formation) # Makes inputted occ data into lists from Formations
formations <- formations[order(formations$Formation),]


Berriasian <- c(145, 139.4)
Valanginian <- c(139.4, 134.7)
Hauterivian <- c(134.7, 130.8)
Barremian <- c(130.8, 126.3)
Aptian <- c(126.3, 113.1)
Albian <- c(113.1, 100.5)
Cenomanian <- c(100.5, 93.9)
Turonian <- c(93.9, 89.8)
Coniacian <- c(89.8, 86.3)
Santonian <- c(86.3, 83.6)
Campanian <- c(83.6, 72.1)
Maastrichtian <- c(72.1, 66)

stages <- c(145, 139.4, 134.7, 130.8, 126.3, 113.1, 100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66)
stages <- c(100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66)
stagelist <- list(Cenomanian, Turonian, Coniacian, Santonian, Campanian, Maastrichtian)

formations$forbinning <- 1:nrow(formations)

# FUNCTIONS #

GetInfo <- function(formations){
  mean <- rowMeans(subset(formations, select = c(3, 4)), na.rm = TRUE)
  range <- formations$max_age - formations$min_age
  par(mfrow=c(2,1))
  hist(mean, main = "Mean Age of Formations",xlab = "Mean Age (Ma)")
  hist(range, main = "Range of Formations (Ma)", xlab = "Range (Ma)")
  cat ('Mean range of formations: ', mean(range), fill = TRUE)
  cat ('Median range of formations: ', median(range),fill = TRUE)
}

Scoring_Grid1 <- function(formations, res=0.01) { #includes all formations
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  bins <- seq(min_age-1.0045, max_age+1.0045, res) # Makes 1ma bins in sequence based on max/min ages. 0.0001 added to ensure formation is never exactly equivalent to a bin.

  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(bins)) # makes a matrix for the scoring 
                               # of each time line in terms of how good it is to be a bin boundary
  colnames(score_grid) <- bins # All time lines
  rownames(score_grid) <- formations$Formation # All formations

  counter <- 0
  for(i in bins) { # Go through each time line
    counter <- sum(counter,1)
    for (f in 1:nrow(formations)){ # go through each formation 
      if (i <= formations$max_age[f] && i >= formations$min_age[f]){ # if timeline is between max/min age of a formation (i.e. formation crosses that line)
        # Need to translate i to position in grid - think this should work
        a <- formations$max_age[f] - i
        b <- i - formations$min_age[f]
        range <- formations$max_age[f] - formations$min_age[f]
        if (a > b){
          score_grid[,counter][f] <- (a/range)*100
        }
        else{ 
          score_grid[,counter][f] <- (b/range)*100 # Work out percentage that sits each side of line, reduce score by that amount.
        }
      }
      else {
        score_grid[,counter][f] = 100 # Otherwise, just score it 100. 
      }
    }
  }  
  means <- colMeans(score_grid) # Work out mean score for each time bin
  score_grid <- rbind(score_grid, means) # add to grid
  score_grid <<- score_grid
  bins <<- bins
}

Scoring_Grid2 <- function(formations, res=0.01) { # ignores formations longer than mean range
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  bins <- seq(min_age-1.0045, max_age+1.0045, res) # Makes 1ma bins in sequence based on max/min ages. 0.0001 added to ensure formation is never exactly equivalent to a bin.

  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(bins)) # makes a matrix for the scoring 
                               # of each time line in terms of how good it is to be a bin boundary
  colnames(score_grid) <- bins # All time lines
  rownames(score_grid) <- formations$Formation # All formations

  counter <- 0
  for(i in bins) { # Go through each time line
    counter <- sum(counter,1)
    for (f in 1:nrow(formations)){ # go through each formation 
      if (formations$max_age[f] - formations$min_age[f] < 
          mean(formations$max_age - formations$min_age)) {
         if (i <= formations$max_age[f] && i >= formations$min_age[f]){ # if timeline is between max/min age of a formation (i.e. formation crosses that line)
            # Need to translate i to position in grid - think this should work
            a <- formations$max_age[f] - i
            b <- i - formations$min_age[f]
            range <- formations$max_age[f] - formations$min_age[f]
            if (a > b){
              score_grid[,counter][f] <- (a/range)*100
            }
            else{ 
              score_grid[,counter][f] <- (b/range)*100 # Work out percentage that sits each side of line, reduce score by that amount.
            }
          }
          else {
            score_grid[,counter][f] = 100 # Otherwise, just score it 100. 
          }
      }
      else {
        next
      }
    }
  }  
  score_grid <- na.omit(score_grid)
  means <- colMeans(score_grid) # Work out mean score for each time bin
  score_grid <- rbind(score_grid, means) # add to grid
  score_grid <<- score_grid
  bins <<- bins
}

binsize <- 3

newBins <- function(score_grid, formations, bin_size, bins, stages){
  score_grid2<- as.data.frame(score_grid)
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- 66
  newbins_test <- c(min_age)
  testbin <- seq(min_age, max_age, binsize)
  for (i in 1:length((testbin)-1)){
    seqs <- seq(testbin[i],testbin[i]+binsize, 1)
  testmatch <- paste("^",seqs[1],"|", "^",seqs[2],"|","^",seqs[3],"|",
                     "^",seqs[4],"|",
                     "^",seqs[5], sep = "")
  a <- score_grid2[grep(testmatch, names(score_grid2))] 
  z <- apply(a,1,which.max) 
  newbins_test[i+1] <- names(a)[z][nrow(a)]
  }
  newbins_test <<- newbins_test
  newbinlist <- list()
  for (i in 1:(length(newbins_test)-1)) {
    newbinlist[[i]] <- c(as.numeric(newbins_test[i]), as.numeric(newbins_test[i+1]))
  }
  newbinlist <<- newbinlist
  par(mfrow=c(1,1))
  plot(bins, colMeans(score_grid2), type = 'l', main = 'Selected Bins', xlim=rev(range.default(as.numeric(newbins_test))), xlab="Time (Million Years)", ylab = "Bin Splitting Score") # plot, obviously. 
  abline(v=newbins_test, lty=2, col="blue")
  for(n in 1:length(stages)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(stages)){
        rect(stages[n], 90, stages[n]-1, 100.5, col = "#32323232", border = NA)
      }
      else{
      rect(stages[n], 90, stages[n+1], 100.5, col = "#32323232", border = NA)
      }
    }
  }
}

FormationGraph <- function(formations, newbins_test, stages){
  fp <- data.frame("x1" = formations$min_age, "y1" = formations$forbinning, 'x2' = formations$max_age, 'y2' = formations$forbinning)
  plot(range(fp$x1, fp$x2), range(fp$y1, fp$y2),  xlab="Age (Ma)", ylab="Individual formation",
  type = 'n', xlim=rev(range.default(66:max(fp$x2))))
  segments(fp$x1, fp$y1, fp$x2, fp$y2, lwd = 2)
  segments
  #abline(v=newbins_test, lty=2, col="blue")
  #abline(v=stages, lty=2, col="red")
  for(n in 1:length(stages)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(stages)){
        rect(stages[n], -2, stages[n]-1, max(fp$y2)+3, col = "#32323232", border = NA)
      }
      else{
      rect(stages[n], -2, stages[n+1], max(fp$y2)+3, col = "#32323232", border = NA)
      }
    }
  }
}

Scoring_Grid1(formations)
newBins(score_grid, formations, bin_size, bins, stages)
FormationGraph(formations, newbins_test, stages)

### CHECK FOR REPETITIONS ###
# Score1
newbinlist[[10]] <- NULL
newbinlist[[13]] <- NULL
# Score2
newbinlist[[6]] <- NULL
newbinlist[[9]] <- NULL

### METHOD 1: ASSIGN FORMATIONS AND SPECIES TO ALL BINS THEY APPEAR IN ###
FormBin_M1 <- function(formations, newbinlist, Form_list) {
  M1_Dino_List <- list()# make an empty list of the dinos in each bin
  
  for(b in 1:length(newbinlist)){ # for each new bin
    temp_dino_recs <- data.frame()
    for (f in 1:nrow(formations)){ 
      if (formations$max_age[f] > newbinlist[[b]][1] && formations$min_age[f] < newbinlist[[b]][2]){ # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
        temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
        print(formations$Formation[f])
      }
    }
    M1_Dino_List[[b]] <- temp_dino_recs
  }
  
  #Score 1
  newbins_test <- c(newbins_test[1:9], newbins_test[11:14])
  #Score 2
  newbins_test <- c(newbins_test[1:6], newbins_test[8:10], newbins_test[12:15])
  
  prefix <- "FB."
  suffix <- seq(1:(length(newbinlist)))
  my_names <- paste(prefix, suffix, sep = "")
  
  binlist <- data.frame(bin = my_names, 
                      bin_min = as.numeric(newbins_test[1:(length(newbins_test)-1)]), 
                      bin_max = as.numeric(newbins_test[2:(length(newbins_test))]))
  binlist$midpoint <- (binlist$bin_min + binlist$bin_max) / 2
  
  M1_incidence_data <- list()
  
  for (n in 1:length(M1_Dino_List)){
    incidence_raw <- M1_Dino_List[[n]] %>% .[,c("collection_no","occurrence.genus_name")] %>% distinct %>% table %>% t
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1]))
    M1_incidence_data[[n]] <- incidence_raw
  }
}

M1_incidence_data <- type.convert(M1_incidence_data)
names(M1_incidence_data) <- my_names

M1_raw_dinos <- unlist(lapply(1:length(M1_Dino_List), function(i) {
  nrow(unique(M1_Dino_List[[i]][6]))
}))

Ale_1 <- c(65, 95)

 par(mfrow=c(1,1))
  plot(binlist$midpoint, M1_raw_dinos, type = 'l', main = 'Selected Bins', xlim=rev(range.default(as.numeric(Ale_1))), xlab="Time (Ma)", ylab = "Raw Generic Diversity") # plot, obviously. 
  abline(v=newbins_test, lty=2, col="blue")
  for(n in 1:length(stages)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(stages)){
        rect(stages[n], 90, stages[n]-1, 100.5, col = "#32323232", border = NA)
      }
      else{
      rect(stages[n], 90, stages[n+1], 100.5, col = "#32323232", border = NA)
      }
    }
  }
  
## First make a list of quorum levels from 0.1-0.9
quorum_levels <- round(seq(from = 0.1, to = 0.9, by = 0.1), 1)

## Now, run the loop using the package pbmclapply
## If you have a large dataset, this part will put your computer under pressure...

  M1_estD_output3 <- estimateD(M1_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[3])
  M1_estD_output4 <- estimateD(M1_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[4])
  M1_estD_output5 <- estimateD(M1_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[5])
  M1_estD_output6 <- estimateD(M1_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[6])
  M1_estD_output7 <- estimateD(M1_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[7])
  
  M1_estD_output4 <- M1_estD_output4[M1_estD_output5$order == 0, ]
  M1_estD_output5 <- M1_estD_output5[M1_estD_output5$order == 0, ] #filter to just species richness
  M1_estD_output6 <- M1_estD_output6[M1_estD_output6$order == 0, ] #filter to just species richness
  M1_estD_output7 <- M1_estD_output7[M1_estD_output7$order == 0, ] #filter to just species richness
  
  M1_estD_output4$reference_t <- sapply(M1_incidence_data, sum)
  M1_estD_output5$reference_t <- sapply(M1_incidence_data, sum)
  M1_estD_output6$reference_t <- sapply(M1_incidence_data, sum)
  M1_estD_output7$reference_t <- sapply(M1_incidence_data, sum) #tally total occurrences in each bin

  M1_estD_output5[which(M1_estD_output4$t >= 2 * M1_estD_output5$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3)
  M1_estD_output5[which(M1_estD_output5$t >= 2 * M1_estD_output5$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) #no more than twice reference sample size
  M1_estD_output6[which(M1_estD_output6$t >= 2 * M1_estD_output6$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) 
  M1_estD_output7[which(M1_estD_output7$t >= 2 * M1_estD_output7$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) 
  
  M1_estD_output4$quorum_level <- quorum_levels[4]  
  M1_estD_output5$quorum_level <- quorum_levels[5] #create new column
  M1_estD_output6$quorum_level <- quorum_levels[6] #create new column
  M1_estD_output7$quorum_level <- quorum_levels[7] #create new column

plot(rev(M1_estD_output6$qD), type = "line")
## Take the output and make it more manageable:
M1_subsampled_richness <- bind_rows(M1_estD_output4, M1_estD_output5, M1_estD_output6, M1_estD_output7) #binds lists of dataframes

## add intervals names
M1_subsampled_richness <- M1_subsampled_richness %>% rename(bin = site) %>% full_join(., binlist, by = "bin") #NB: if this throws up a warning message its not a problem!

## ensuring the variables below are treated as factors is important for easy of plotting!
M1_subsampled_richness$stages <- as.factor(M1_subsampled_richness$bin) 
M1_subsampled_richness$quorum_level <- as.factor(M1_subsampled_richness$quorum_level) 

## Pullinng the quorum levels of interest (typically 0.3-0.7) out into separate lists can make plotting in ggplot much easier
M1_subsampled_richness_4 <- subset(M1_subsampled_richness, quorum_level == 0.4)
M1_subsampled_richness_5 <- subset(M1_subsampled_richness, quorum_level == 0.5)
M1_subsampled_richness_6 <- subset(M1_subsampled_richness, quorum_level == 0.6)

M1_subsampled_richnes_plot <- ggplot(filter(M1_subsampled_richness, quorum_level %in% quorum_levels[4:6]), aes(x = bin_min, y = qD, ymin = qD.LCL, ymax = qD.UCL, xmin = bin_min, xmax = bin_max)) + 
  geom_ribbon(data=M1_subsampled_richness_4, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "royalblue", alpha = 0.2) +
  geom_ribbon(data=M1_subsampled_richness_5, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "steelblue1", alpha = 0.2) +
  geom_ribbon(data=M1_subsampled_richness_6, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "skyblue", alpha = 0.2) +
  geom_line(data=M1_subsampled_richness_4, aes(x = midpoint, y = qD), colour = 'royalblue', size = 1.2) +
  geom_line(data=M1_subsampled_richness_5, aes(x = midpoint, y = qD), colour = 'steelblue1', size = 1.2) +
  geom_line(data=M1_subsampled_richness_6, aes(x = midpoint, y = qD), colour = 'skyblue', size = 1.2) +
  #geom_point(aes(pch = method), size = 4.5) +
  scale_x_reverse(expand=c(0,0)) +
  scale_y_continuous(trans = "log10", expand=c(0,0)) +
  labs(x = "Time (Ma)", y = "Coverage rarified richness (log scale)") +
  theme(panel.background = element_blank(),
        #legend.position="none",
        #plot.margin = margin(2, 2, 2, 2, "cm"),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=14, angle=0, hjust=0.5),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=12),
        aspect.ratio=1)
M1_subsampled_richnes_plot


### METHOD 2: ASSIGN FORMATIONS AND SPECIES TO 1 BIN, BASED ON % OF FORMATION IN BIN ###


FormBin_M2<- function(formations, binlist, Form_list) {
  M2_Dino_List <- list()# make an empty list of the dinos in each bin

  for(b in 1:length(newbinlist)){ # for each new bin
    temp_dino_recs <- data.frame()
    for (f in 1:nrow(formations)) { # for each formation
        if (formations$max_age[f] < newbinlist[[b]][1] |
            newbinlist[[b]][2] < formations$min_age[f]) { # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
        next # Just print, don't do anything else.
      }
      else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
          if (formations$max_age[f] < newbinlist[[b]][2] && 
              formations$min_age[f] > newbinlist[[b]][1]){# If formation sits within boundaries
          temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
          next
        }
          if (formations$max_age[f] > newbinlist[[b]][1] && 
              formations$min_age[f] < newbinlist[[b]][1] && 
              formations$max_age[f] < newbinlist[[b]][2]){ #If formation crosses upper age limit
          x <- as.numeric(formations$max_age[f] - as.numeric(newbinlist[[b]][1]))
          y <- as.numeric(as.numeric(newbinlist[[b]][1]) - formations$min_age[f])
          if (x > y){
            temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
            next
          }
        }
          if (formations$max_age[f] > newbinlist[[b]][2] && 
              formations$min_age[f] < newbinlist[[b]][2] &&
              formations$min_age[f] > newbinlist[[b]][1]){
          x <- as.numeric(formations$max_age[f] - as.numeric(newbinlist[[b]][2]))
          y <- as.numeric(as.numeric(newbinlist[[b]][2]) - formations$min_age[f])
          if (y > x){
            temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
            next
          }
        }
          if (formations$max_age[f] > newbinlist[[b]][2] &&
              formations$min_age[f] < newbinlist[[b]][1]){
              temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
              next
        }
      }
    }
    M2_Dino_List[[b]] <- temp_dino_recs
  }

}

M2_raw_dinos <- unlist(lapply(1:length(M2_Dino_List), function(i) {
  nrow(unique(M2_Dino_List[[i]][6]))
}))
 
  M2_Dino_List[[8]] <- NULL
  M2_incidence_data <- list()
  
  for (n in 1:length(M2_Dino_List)){
    if (nrow(M2_Dino_List[[n]]) == 0){
      next
    }
    else{
    incidence_raw <- M2_Dino_List[[n]] %>% .[,c("collection_no","occurrence.genus_name")] %>% distinct %>% table %>% t
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1]))
    M2_incidence_data[[n]] <- incidence_raw
  }
}

data(ciliates)

M2_incidence_data <- type.convert(M2_incidence_data)
names(M2_incidence_data) <- my_names
M2_incidence_data[[8]] <- NULL

M1_raw_dinos <- unlist(lapply(1:length(M1_Dino_List), function(i) {
  nrow(unique(M1_Dino_List[[i]][6]))
}))
  
## First make a list of quorum levels from 0.1-0.9
quorum_levels <- round(seq(from = 0.1, to = 0.9, by = 0.1), 1)

## Now, run the loop using the package pbmclapply
## If you have a large dataset, this part will put your computer under pressure...

  M2_estD_output3 <- estimateD(M2_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[3])
  M2_estD_output4 <- estimateD(M2_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[4])
  M2_estD_output5 <- estimateD(M2_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[5])
  M2_estD_output6 <- estimateD(M2_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[6])
  M2_estD_output7 <- estimateD(M2_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[7])
  
  M2_estD_output4 <- M2_estD_output4[M2_estD_output5$order == 0, ]
  M2_estD_output5 <- M2_estD_output5[M2_estD_output5$order == 0, ] #filter to just species richness
  M2_estD_output6 <- M2_estD_output6[M2_estD_output6$order == 0, ] #filter to just species richness
  M2_estD_output7 <- M2_estD_output7[M2_estD_output7$order == 0, ] #filter to just species richness
  
  estD_output4$reference_t <- sapply(incidence_data, sum)
  estD_output5$reference_t <- sapply(incidence_data, sum)
  estD_output6$reference_t <- sapply(incidence_data, sum)
  estD_output7$reference_t <- sapply(incidence_data, sum) #tally total occurrences in each bin

  estD_output5[which(estD_output4$t >= 2 * estD_output5$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3)
  estD_output5[which(estD_output5$t >= 2 * estD_output5$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) #no more than twice reference sample size
  estD_output6[which(estD_output6$t >= 2 * estD_output6$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) 
  estD_output7[which(estD_output7$t >= 2 * estD_output7$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) 
  
  estD_output4$quorum_level <- quorum_levels[4]  
  estD_output5$quorum_level <- quorum_levels[5] #create new column
  estD_output6$quorum_level <- quorum_levels[6] #create new column
  estD_output7$quorum_level <- quorum_levels[7] #create new column

## Take the output and make it more manageable:
subsampled_richness <- bind_rows(estD_output4, estD_output5, estD_output6, estD_output7) #binds lists of dataframes

## add intervals names
subsampled_richness <- subsampled_richness %>% rename(bin = site) %>% full_join(., binlist, by = "bin") #NB: if this throws up a warning message its not a problem!

## ensuring the variables below are treated as factors is important for easy of plotting!
subsampled_richness$stages <- as.factor(subsampled_richness$bin) 
subsampled_richness$quorum_level <- as.factor(subsampled_richness$quorum_level) 

## Pullinng the quorum levels of interest (typically 0.3-0.7) out into separate lists can make plotting in ggplot much easier
subsampled_richness_4 <- subset(subsampled_richness, quorum_level == 0.4)
subsampled_richness_5 <- subset(subsampled_richness, quorum_level == 0.5)
subsampled_richness_6 <- subset(subsampled_richness, quorum_level == 0.6)

subsampled_richnes_plot <- ggplot(filter(subsampled_richness, quorum_level %in% quorum_levels[4:6]), aes(x = bin_min, y = qD, ymin = qD.LCL, ymax = qD.UCL, xmin = bin_min, xmax = bin_max)) + 
  geom_ribbon(data=subsampled_richness_4, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "royalblue", alpha = 0.2) +
  geom_ribbon(data=subsampled_richness_5, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "steelblue1", alpha = 0.2) +
  geom_ribbon(data=subsampled_richness_6, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "skyblue", alpha = 0.2) +
  geom_line(data=subsampled_richness_4, aes(x = midpoint, y = qD), colour = 'royalblue', size = 1.2) +
  geom_line(data=subsampled_richness_5, aes(x = midpoint, y = qD), colour = 'steelblue1', size = 1.2) +
  geom_line(data=subsampled_richness_6, aes(x = midpoint, y = qD), colour = 'skyblue', size = 1.2) +
  #geom_point(aes(pch = method), size = 4.5) +
  scale_x_reverse(expand=c(0,0)) +
  scale_y_continuous(trans = "log10", expand=c(0,0)) +
  labs(x = "", y = "Coverage rarified richness (log scale)") +
  theme(panel.background = element_blank(),
        #legend.position="none",
        #plot.margin = margin(2, 2, 2, 2, "cm"),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=14, angle=0, hjust=0.5),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=12),
        aspect.ratio=1)
subsampled_richnes_plot

### METHOD 3: ASSIGN % OF DINOSAUR OCCS AT RANDOM TO A BIN BASED ON % OF FOPRMATION IN BIN - REPEAT AND AVERAGE ###

FormBin_M3<- function(formations, binlist, Form_list, times=10) {
  Bin_Dino_List3 <- list()
  M3_Dino_List <- list()
  # make an empty list of the dinos in each bin
  av_species <- as.list(rep(0, length(binlist)))
  
  for(n in 1:times){
    cat("cycle: ", n)

    for(b in 1:length(newbinlist)){ # for each new bin
      print(b)
      temp_dino_recs <- data.frame()
      for (f in 1:nrow(formations)) { # for each formation
        if (formations$max_age[f] < newbinlist[[b]][1] |
            newbinlist[[b]][2] < formations$min_age[f]) { # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
          print("not in bin")
          next # Just print, don't do anything else.
        }
        else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
          print(formations$Formation[f])  
          if (formations$max_age[f] < newbinlist[[b]][2] && 
              formations$min_age[f] > newbinlist[[b]][1]){ # If formation sits within boundaries
            temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]])
            print("Within bin") # Add dinos from that formation to dino list
            next
          }
          if (formations$max_age[f] > newbinlist[[b]][1] && 
              formations$min_age[f] < newbinlist[[b]][1] && 
              formations$max_age[f] < newbinlist[[b]][2]){#If formation crosses upper age limit
            x <- as.numeric(as.numeric(formations$max_age[f] - newbinlist[[b]][1]))
            y <- as.numeric(formations$max_age[f] - formations$min_age[f])    
            perc <- (x/y)
            ran_num <- perc*nrow(Form_list[[f]])
            temp_dframe <- Form_list[[f]]
            sampled <- temp_dframe[sample(nrow(temp_dframe), ran_num, replace = FALSE),]
            temp_dino_recs <- rbind(temp_dino_recs, sampled) # Add dinos from that formation to dino list
            print("Top bin boundary")
            next
          } # This bit (above) should be done
          if (formations$max_age[f] > newbinlist[[b]][2] && 
              formations$min_age[f] < newbinlist[[b]][2] &&
              formations$min_age[f] > newbinlist[[b]][1]){
            x <- as.numeric(as.numeric(newbinlist[[b]][2]) - formations$min_age[f])
            y <- as.numeric(formations$max_age[f] - formations$min_age[f])
            perc <- (x/y)
            ran_num <- perc*nrow(Form_list[[f]])
            temp_dframe <- Form_list[[f]]
            sampled <- temp_dframe[sample(nrow(temp_dframe), ran_num, replace = FALSE),]
            temp_dino_recs <- rbind(temp_dino_recs, sampled) # Add dinos from that formation to dino list
            print("Bottom bin boundary")
            next
          } # This bit (above) should be done
          if (formations$max_age[f] > newbinlist[[b]][2] &&
              formations$min_age[f] < newbinlist[[b]][1]){
            x <- as.numeric(formations$max_age[f] - formations$min_age[f])
            y <- as.numeric(as.numeric(newbinlist[[b]][2]) - as.numeric(newbinlist[[b]][1]))
            perc <- (y/x)
            ran_num <- perc*nrow(Form_list[[f]])
            temp_dframe <- Form_list[[f]]
            sampled <- temp_dframe[sample(nrow(temp_dframe), ceiling(ran_num), replace = FALSE),]
            temp_dino_recs <- rbind(temp_dino_recs, sampled) # Add dinos from that formation to dino list
            print("Top and Bottom boundary")
            next
          } # This bit (above) should be done
        }
      }
      M3_Dino_List[[b]] <- temp_dino_recs
    }
    
M3_raw_dinos <- unlist(lapply(1:length(M3_Dino_List), function(i) {
  nrow(unique(M3_Dino_List[[i]][6]))
}))
  

M3_incidence_data <- list()
  for (n in 1:length(M3_Dino_List)){
    incidence_raw <- M3_Dino_List[[n]] %>% .[,c("collection_no","occurrence.genus_name")] %>% distinct %>% table %>% t
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1]))
    M3_incidence_data[[n]] <- incidence_raw
  }
  

M3_incidence_data <- type.convert(M3_incidence_data)
names(M3_incidence_data) <- my_names

M3_incidence_data[[8]] <- NULL
M3_incidence_data[[11]] <- NULL

incidence_data2 <- M3_incidence_data
for (n in 1:length(incidence_data2)){
  if (length(M3_incidence_data[[n]]) == 0){
    M3_incidence_data[[n]] <- NULL
  }
}

  
  ## First make a list of quorum levels from 0.1-0.9
quorum_levels <- round(seq(from = 0.1, to = 0.9, by = 0.1), 1)

## Now, run the loop using the package pbmclapply
## If you have a large dataset, this part will put your computer under pressure...

  M3_estD_output3 <- estimateD(M3_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[3])
  M3_estD_output4 <- estimateD(M3_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[4])
  M3_estD_output5 <- estimateD(M3_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[5])
  M3_estD_output6 <- estimateD(M3_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[6])
  M3_estD_output7 <- estimateD(M3_incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[7])
  
  estD_output4 <- estD_output4[estD_output5$order == 0, ]
  estD_output5 <- estD_output5[estD_output5$order == 0, ] #filter to just species richness
  estD_output6 <- estD_output6[estD_output6$order == 0, ] #filter to just species richness
  estD_output7 <- estD_output7[estD_output7$order == 0, ] #filter to just species richness
  
  estD_output4$reference_t <- sapply(incidence_data, sum)
  estD_output5$reference_t <- sapply(incidence_data, sum)
  estD_output6$reference_t <- sapply(incidence_data, sum)
  estD_output7$reference_t <- sapply(incidence_data, sum) #tally total occurrences in each bin

  estD_output5[which(estD_output4$t >= 2 * estD_output5$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3)
  estD_output5[which(estD_output5$t >= 2 * estD_output5$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) #no more than twice reference sample size
  estD_output6[which(estD_output6$t >= 2 * estD_output6$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) 
  estD_output7[which(estD_output7$t >= 2 * estD_output7$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) 
  
  estD_output4$quorum_level <- quorum_levels[4]  
  estD_output5$quorum_level <- quorum_levels[5] #create new column
  estD_output6$quorum_level <- quorum_levels[6] #create new column
  estD_output7$quorum_level <- quorum_levels[7] #create new column

## Take the output and make it more manageable:
subsampled_richness <- bind_rows(estD_output4, estD_output5, estD_output6, estD_output7) #binds lists of dataframes

## add intervals names
subsampled_richness <- subsampled_richness %>% rename(bin = site) %>% full_join(., binlist, by = "bin") #NB: if this throws up a warning message its not a problem!

## ensuring the variables below are treated as factors is important for easy of plotting!
subsampled_richness$stages <- as.factor(subsampled_richness$bin) 
subsampled_richness$quorum_level <- as.factor(subsampled_richness$quorum_level) 

## Pullinng the quorum levels of interest (typically 0.3-0.7) out into separate lists can make plotting in ggplot much easier
subsampled_richness_4 <- subset(subsampled_richness, quorum_level == 0.4)
subsampled_richness_5 <- subset(subsampled_richness, quorum_level == 0.5)
subsampled_richness_6 <- subset(subsampled_richness, quorum_level == 0.6)
  
        if(n == 1){ # If the first loop round, just put data in this list
        print("YES")
        M3_Dino_List[[b]] <- temp_dino_recs
      }
      else { # Otherwise, add to what's previously there.
        print("NO")
        M3_Dino_List[[b]] <- c(M3_Dino_List[[b]], temp_dino_recs)
  
    
    for (q in 1:length(Bin_Dino_List3)){
      av_species[[q]] <- c(av_species[[q]],length(Bin_Dino_List3[[q]])) # Find number of species per bin for this cycle
    }   
  }

  av_species <- lapply(av_species, function(x) x[-1])
  Species_M3 <- c() # METHOD 3
  for (t in 1:length(av_species)) {
    Species_M3 <- c(Species_M3, mean(av_species[[t]]))
  }
  Species_M3 <<- Species_M3
  M3_Dino_List <<- M3_Dino_List
