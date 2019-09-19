# FORMATION BINNER 4
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


### Necessary Packages ###
install.packages("tidyverse")
install.packages("pbmcapply")
install.packages("divDyn")
install.packages("dplyr")

library(pbmcapply)
library(dplyr)
library(tidyverse)
library(divDyn)

select <- dplyr::select # ensure the select function is coming from tidyverse package

### Setting up Data ###

# Data input #
#setwd(dir = "/Users/ChristopherDean/MEGA/PhD/PROJECTS/FRM_BIN/Data/") # Set working Directory - Mac Version
setwd(dir = "F:/PhD Work/MEGA/MEGA2/PhD/PROJECTS/FRM_BIN/Data") # Set working Directory - PC version
formations <- read.csv (file = "Formations_test2.csv")  #Read in formations
occs <- read.csv(file = "NADINOS-occs-edit.csv") # Read in occurrences
#Berriasian <- c(145, 139.4)
#Valanginian <- c(139.4, 134.7)
#Hauterivian <- c(134.7, 130.8)
#Barremian <- c(130.8, 126.3)
#Aptian <- c(126.3, 113.1)
#Albian <- c(113.1, 100.5)
Cenomanian <- c(100.5, 93.9)
Turonian <- c(93.9, 89.8)
Coniacian <- c(89.8, 86.3)
Santonian <- c(86.3, 83.6)
Campanian <- c(83.6, 72.1)
Maastrichtian <- c(72.1, 66)
#stages <- c(145, 139.4, 134.7, 130.8, 126.3, 113.1, 100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66)
stages <- c(100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66)
stagelist <- list(Cenomanian, Turonian, Coniacian, Santonian, Campanian, Maastrichtian)

# Data Cleaning #
formations$max_age <- as.numeric(as.character(formations$max_age)) # Make Numeric
formations$min_age <- as.numeric(as.character(formations$min_age)) # Make Numeric 
formations <- formations[which(formations$Location=='WI'),] # Only formations from Western Interior
myformations <- sort(as.vector(formations$Formation)) # Organise
testoccs <- occs[occs$formation %in% myformations,] # Only include organisms from formation list
testoccs <- droplevels.data.frame(testoccs) # Remove old levels
Form_list <- split(testoccs, testoccs$formation) # Makes inputted occ data into lists from Formations
formations <- formations[order(formations$Formation),] # Reorganise formations
formations$forbinning <- 1:nrow(formations) # Number formations for easy plotting later

### FUNCTIONS ###

# GetInfo - basic information about Formations

GetInfo <- function(formations){
  mean <- rowMeans(subset(formations, select = c(3, 4)), na.rm = TRUE)
  range <- formations$max_age - formations$min_age
  par(mfrow=c(2,1))
  hist(mean, main = "Mean Age of Formations",xlab = "Mean Age (Ma)")
  hist(range, main = "Range of Formations (Ma)", xlab = "Range (Ma)")
  cat ('Mean range of formations: ', mean(range), fill = TRUE)
  cat ('Median range of formations: ', median(range),fill = TRUE)
}

# Scoring_Grid1 - Create a scoring grid using info from all formations

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

# Scoring Grid2 - Create a scorign grid ignoring formations with length longer than mean formation length. 

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


# newBins - Create new bins based on scoring grid, outputted in a variety of ways 

newBins <- function(score_grid, formations, bin_size, bins, stages){
  score_grid2<- as.data.frame(score_grid)
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- 66
  newbins_test <- c(min_age)
  testbin <- seq(min_age, max_age, binsize)
  for (i in 1:length((testbin)-1)){
    seqs <- seq(testbin[i],testbin[i]+binsize, 1)
    pasting <- c()
    for (n in 1:binsize){
      pasting <- c(pasting, paste("^",seqs[n],"|", sep = ""))
    }
    testmatch2 <- paste(pasting, collapse = "")
    testmatch2 <- substr(testmatch2, 1, nchar(testmatch2)-1) 
    a <- score_grid2[grep(testmatch2, names(score_grid2))] 
    z <- apply(a,1,which.max) 
    newbins_test[i+1] <- names(a)[z][nrow(a)]
  }
  newbins_test <- unique(newbins_test)
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
  prefix <- "FB."
  suffix <- seq(1:(length(newbinlist)))
  my_names <- paste(prefix, suffix, sep = "")
  binlist <- data.frame(bin = my_names, 
                      bottom = as.numeric(newbins_test[1:(length(newbins_test)-1)]), 
                      top = as.numeric(newbins_test[2:(length(newbins_test))]))
  binlist$mid <- (binlist$bottom + binlist$top) / 2
  binlist <<- binlist
}

# FormationGraph - Show what formations look like through time in comparison to Stages and new Bins.

FormationGraph <- function(formations, newbins_test, stages){
  fp <- data.frame("x1" = formations$min_age, "y1" = formations$forbinning, 'x2' = formations$max_age, 'y2' = formations$forbinning)
  plot(range(fp$x1, fp$x2), range(fp$y1, fp$y2),  xlab="Age (Ma)", ylab="Individual formation",
  type = 'n', xlim=rev(range.default(66:max(fp$x2))))
  segments(fp$x1, fp$y1, fp$x2, fp$y2, lwd = 2)
  segments
  abline(v=newbins_test, lty=2, col="blue")
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

### Running first results ###

binsize <- 3 # Set user defined bin size
Scoring_Grid1(formations, res = 0.5)
newBins(score_grid, formations, bin_size, bins, stages)
FormationGraph(formations, newbins_test, stages)

### METHOD 1: ASSIGN FORMATIONS AND SPECIES TO ALL BINS THEY APPEAR IN ###
FormBin_M1 <- function(formations, newbinlist, binlist, Form_list, Quorum=0.4) {
  M1_Dino_List <- list()# make an empty list of the dinos in each bin
  
  for(b in 1:length(newbinlist)){ # for each new bin
    temp_dino_recs <- data.frame()
    for (f in 1:nrow(formations)){ 
      if (formations$max_age[f] > newbinlist[[b]][1] && formations$min_age[f] < newbinlist[[b]][2]){ # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
        temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
      }
    }
    temp_dino_recs$bin_no <- b
    M1_Dino_List[[b]] <- temp_dino_recs
  }
  
  df <- do.call("rbind", M1_Dino_List)
  SQS <- subsample(df,iter=100, q=Quorum,tax="occurrence.genus_name", bin="bin_no", 
                   coll = 'collection_no', output="dist", type="sqs", useFailed = TRUE)
  tsplot(binlist,
  ylab="range-through diversity (genera)", ylim=c(0,25))
  shades(binlist$mid, SQS$divSIB, col="black")
}

### METHOD 2: ASSIGN FORMATIONS AND SPECIES TO 1 BIN, BASED ON % OF FORMATION IN BIN ###
FormBin_M2<- function(formations, newbinlist, binlist, Form_list, Quorum = 0.4) {
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
    if (nrow(temp_dino_recs) > 0){
      temp_dino_recs$bin_no <- b
    }
    M2_Dino_List[[b]] <- temp_dino_recs
  }
  df <- do.call("rbind", M2_Dino_List)
  SQS <- subsample(df,iter=100, q=Quorum,tax="occurrence.genus_name", bin="bin_no", 
                   coll = 'collection_no', output="dist", type="sqs", useFailed = TRUE)
  if (nrow(SQS$divSIB) < length(newbinlist) ){ # If last bin has failed
      SQS$divSIB[nrow(SQS$divSIB) + 1,] <- 0 # Add additional row of 0's for that bin
  }
  SQS <<- SQS
  tsplot(binlist,
  ylab="range-through diversity (genera)", ylim=c(0,25))
  shades(binlist$mid, SQS$divSIB, col="black")
}

### METHOD 3: ASSIGN % OF DINOSAUR OCCS AT RANDOM TO A BIN BASED ON % OF FOPRMATION IN BIN - REPEAT AND AVERAGE ###
times <- 10
quorum <- 0.4
FormBin_M3<- function(formations, newbinlist, binlist, Form_list, times=10, quorum=0.4) {
  M3_Dino_List <- list()
  for(n in 1:times){
    tryCatch({ # Skips run of a loop if it encounters the known SQS error (resulting from single occurrences in 1 bin)
    cat("cycle:", n)
    for(b in 1:length(newbinlist)){ # for each new bin
      print(b)
      temp_dino_recs <- data.frame()
      for (f in 1:nrow(formations)) { # for each formation
        if (formations$max_age[f] < newbinlist[[b]][1] |
        newbinlist[[b]][2] < formations$min_age[f]) { # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
          next # Just print, don't do anything else.
        }
        else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
          if (formations$max_age[f] < newbinlist[[b]][2] && 
          formations$min_age[f] > newbinlist[[b]][1]){ # If formation sits within boundaries
            temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
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
            next
          } # This bit (above) should be done
        }
      }
      if (nrow(temp_dino_recs) > 0){ # If there are occurrences in this bin
        temp_dino_recs$bin_no <- b # Label those occurrences with bin number
      }
      M3_Dino_List[[b]] <- temp_dino_recs # Add temp records to permanent list
    }

    df <- do.call("rbind", M3_Dino_List) # Create data.frame from list
    SQS <- subsample(df,iter=100, q=quorum,tax="occurrence.genus_name", bin="bin_no", 
                   coll = 'collection_no', output="dist", type="sqs", useFailed = TRUE) # Run SQS
  
    if (nrow(SQS$divSIB) < length(newbinlist) ){ # If last bin has failed
      SQS$divSIB[nrow(SQS$divSIB) + 1,] <- 0 # Add additional row of 0's for that bin
    }
    if(n == 1){ # If the first loop round, just put data in this list
      Final_SQS <- SQS$divSIB
    }
    else { # Otherwise, add to what's previously there.
      Final_SQS <- cbind(Final_SQS, SQS$divSIB) # Add to complete list of SQS results
    }
    }, error=function(e){print("ERROR: NOT ENOUGH DATA")})
  }
  tsplot(binlist,
  ylab="Sampled in Bin diversity (genera)", ylim=c(0,25))
  shades(binlist$mid, Final_SQS, col="black") # Create plot
  Final_SQS <<- Final_SQS
}

### Running Methods ###
FormBin_M1(formations, newbinlist, binlist, Form_list)
FormBin_M2(formations, newbinlist, binlist, Form_list)
FormBin_M3(formations, newbinlist, binlist, Form_list)

# Testing Resolution

Scoring_Grid1(formations, res=0.05)
newBins(score_grid, formations, bin_size, bins, stages)

Scoring_Grid1(formations, res=0.1)
newBins(score_grid, formations, bin_size, bins, stages)

Scoring_Grid1(formations, res=0.5)
newBins(score_grid, formations, bin_size, bins, stages)

Scoring_Grid1(formations, res=1)
newBins(score_grid, formations, bin_size, bins, stages)

# Comparing Scoring Grid Methods
Scoring_Grid1(formations, res=0.5)
newBins(score_grid, formations, bin_size, bins, stages)
FormationGraph(formations, newbins_test, stages)

Scoring_Grid2(formations, res=0.5)
newBins(score_grid, formations, bin_size, bins, stages)
FormationGraph(formations, newbins_test, stages)
