# FORMATION BINNER
# Written by Christopher Dean, 2016-2017
# 04/05/17
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
#                The function looks through time at intervals specified by the user (res) and gives a suitability
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


# Data input #
formations <- read.csv(file = "0.Data/Formations2.csv")  #Read in formations
occs <- read.csv(file = "0.Data/Occurrences.csv") # Read in occurrences

Form_list <- split(occs, occs$Formation) # Makes inputted occ data into lists from Formations

# FUNCTIONS #
Scoring_Grid <- function(formations, res) { 
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  bins <- seq(min_age, max_age+1, res) # Makes 1ma bins in sequence based on max/min ages.

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

Bins <- function(score_grid, thresh, bins) {
    threshold <- thresh
    threshold_list <- c()
    output <- vector()
   for(m in 2:(ncol(score_grid)-1)) {
     if (score_grid[26,m] > threshold) {
       if (score_grid[26,(m)] >= score_grid[26,(m-1)] && score_grid[26,(m)] >= score_grid[26,(m+1)]) {
         if (score_grid[26,(m)] != score_grid[26,(m-1)]) {
           threshold_list <- c(threshold_list,(colnames(score_grid)[m]))
         }
        }
      }
    }
    cat ('Time bins above selected threshold are: ')
    cat (threshold_list, sep=', ')
  
    plot(bins, colMeans(score_grid), type = 'l', xlim=rev(range.default(bins))) # plot, obviously. 
    abline(v=threshold_list, lty=2, col="red")

    threshold_list <- rev(threshold_list)         
    binlist <- list()
    for (i in 1:length(threshold_list)) {
      binlist[[i]] <- c(as.numeric(threshold_list[i]), as.numeric(threshold_list[i+1]))
    }
    binlist[[length(binlist)]][2] <- as.numeric(min(formations$min_age))
    binlist <<- binlist
}

### METHOD 1: ASSIGN FORMATIONS AND SPECIES TO ALL BINS THEY APPEAR IN ###
FormBin_M1 <- function(formations, binlist, Form_list) {
  Bin_Dino_List1 <- list() # make an empty list of the dinos in each bin
  Max_Outcrop_list1 <- list()
  Min_Outcrop_list1 <- list()
  Formation_Counter_List1 <- list() # make an empty list that's for counting formations

  for(b in 1:length(binlist)){ # for each new bin
    temp_dino_list <- c()
    min_outcrop_list <- c()
    max_outcrop_list <- c()
    form_counter <- 0 
    for (f in 1:nrow(formations)){ # for each formation
      if (formations$max_age[f] > binlist[[b]][1] && formations$min_age[f] > binlist[[b]][1] | 
      binlist[[b]][2] > formations$max_age[f]  && binlist[[b]][2] > formations$min_age[f]){ # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
        print(f) # Just print, don't do anything else.
      }
      else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
        temp_dino_list <- c(temp_dino_list, as.character(Form_list[[f]][[1]])) # Add dinos from that formation to dino list
        max_outcrop_list <- c(max_outcrop_list, formations$max_Outcrop[f])
        min_outcrop_list <- c(min_outcrop_list, formations$min_Outcrop[f])
        form_counter <- form_counter + 1
      }
    }
    Bin_Dino_List1[[b]] <- unique(temp_dino_list) # This bit kinda broken I think. Supposed to add all the unique dinos from a time bin to a final list.
    Max_Outcrop_list1[[b]] <- max_outcrop_list
    Min_Outcrop_list1[[b]] <- min_outcrop_list
    Formation_Counter_List1[[b]] <- form_counter
  }
  Area_M1_Max <<- unlist(lapply(Max_Outcrop_list1,sum))
  Area_M1_Min <<- unlist(lapply(Min_Outcrop_list1,sum))
  Species_M1 <- c() # METHOD 1
  for (f in 1:length(Bin_Dino_List1)) {
    Species_M1 <- c(Species_M1, length(Bin_Dino_List1[[f]]))
  }
  Species_M1 <<- Species_M1
  binpoints <- unlist(lapply(binlist, mean))

  par(mar=c(5, 4, 4, 6) + 0.1)
  plot(binpoints, Species_M1, type = 'l', yaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "black")
  axis(2, ylim=c(0,1),col="black",las=1)
  mtext("Raw Dinosaur Diversity",side=2,line=2.5)
  box()
  par(new=TRUE)
  plot(binpoints, Area_M1_Max, type = 'l', yaxt='n', xaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "red", lty=2)
  lines(binpoints, Area_M1_Min, col="blue", lty=2)
  mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
  axis(4, ylim=c(0,7000), col="black",col.axis="black",las=1)
  mtext("Age (Ma)",side=1,col="black",line=2.5) 
  legend("topleft",legend=c("Dinosaur Diversity","Max. Outcrop Area", "Min. Outcrop Area"),
  text.col=c("black","red","blue"),lty=c(1),col=c("black","red", "blue"),cex=0.7)

  DPSK_M1_Max <- Species_M1/Area_M1_Max
  DPSK_M1_Min <- Species_M1/Area_M1_Min
  
  par(mar=c(5, 4, 4, 6) + 0.1)
  plot(binpoints, DPSK_M1_Max, type = 'l', yaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "black")
  lines(binpoints, DPSK_M1_Min, col = "black")
  axis(2, ylim=c(0,1),col="black")
  mtext("DPSK",side=2,line=2.5)
  box()
  par(new=TRUE)
  plot(binpoints, Area_M1_Max, type = 'l', yaxt='n', xaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "red", lty=2)
  lines(binpoints, Area_M1_Min, col="blue", lty=2)
  mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
  axis(4, ylim=c(0,7000), col="black",col.axis="black",las=1)
  mtext("Age (Ma)",side=1,col="black",line=2.5) 
  
}


### METHOD 2: ASSIGN FORMATIONS AND SPECIES TO 1 BIN, BASED ON % OF FORMATION IN BIN ###


FormBin_M2<- function(formations, binlist, Form_list) {
  Bin_Form_List2 <- list()
  Bin_Dino_List2 <- list() # make an empty list of the dinos in each bin
  Max_Outcrop_list2 <- list()
  Min_Outcrop_list2 <- list()
  Formation_Counter_List2 <- list() 

  for(b in 1:length(binlist)){ # for each new bin
    temp_form_list <- c()
    temp_dino_list <- c()
    min_outcrop_list <- c()
    max_outcrop_list <- c()
    form_counter <- 0
    for (f in 1:nrow(formations)) { # for each formation
      if (formations$max_age[f] > binlist[[b]][1] && formations$min_age[f] > binlist[[b]][1] | 
          binlist[[b]][2] > formations$max_age[f]  && binlist[[b]][2] > formations$min_age[f]) { # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
        print(f) # Just print, don't do anything else.
      }
      else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
        if (formations$max_age[f] <= binlist[[b]][1] && formations$min_age[f] >= binlist[[b]][2]){ # If formation sits within boundaries
          temp_form_list <- c(temp_form_list, as.character(formations$Formation[f])) 
          temp_dino_list <- c(temp_dino_list, as.character(Form_list[[f]][[1]]))
          max_outcrop_list <- c(max_outcrop_list, formations$max_Outcrop[f])
          min_outcrop_list <- c(min_outcrop_list, formations$min_Outcrop[f]) # Add dinos from that formation to dino list
          form_counter <- form_counter + 1 
        }
        if (formations$max_age[f] >= binlist[[b]][1] && formations$min_age[f] >= binlist[[b]][2]){ #If formation crosses upper age limit
          x <- as.numeric(formations$max_age[f] - as.numeric(binlist[[b]][1]))
          y <- as.numeric(as.numeric(binlist[[b]][1]) - formations$min_age[f])
          if (y > x){
            temp_form_list <- c(temp_form_list, as.character(formations$Formation[f]))
            temp_dino_list <- c(temp_dino_list, as.character(Form_list[[f]][[1]])) # Add dinos from that formation to dino list
            max_outcrop_list <- c(max_outcrop_list, formations$max_Outcrop[f])
            min_outcrop_list <- c(min_outcrop_list, formations$min_Outcrop[f])
            form_counter <- form_counter + 1 
          }
        }
        if (formations$max_age[f] <= binlist[[b]][1] && formations$min_age[f] <= binlist[[b]][2]){
          x <- as.numeric(formations$max_age[f] - as.numeric(binlist[[b]][2]))
          y <- as.numeric(as.numeric(binlist[[b]][2]) - formations$min_age[f])
          if (x > y){
            temp_form_list <- c(temp_form_list, as.character(formations$Formation[f]))
            temp_dino_list <- c(temp_dino_list, as.character(Form_list[[f]][[1]])) # Add dinos from that formation to dino list
            max_outcrop_list <- c(max_outcrop_list, formations$max_Outcrop[f])
            min_outcrop_list <- c(min_outcrop_list, formations$min_Outcrop[f])
            form_counter <- form_counter + 1 
          }
        }
        if (formations$max_age[f] >= binlist[[b]][1] && formations$min_age[f] <= binlist[[b]][2]){
          x <- as.numeric(formations$max_age[f] - formations$min_age[f])
          y <- as.numeric(as.numeric(binlist[[b]][1]) - as.numeric(binlist[[b]][2]))
          if (y > (x/2)){
            temp_form_list <- c(temp_form_list, as.character(formations$Formation[f]))
            temp_dino_list <- c(temp_dino_list, as.character(Form_list[[f]][[1]])) # Add dinos from that formation to dino list
            max_outcrop_list <- c(max_outcrop_list, formations$max_Outcrop[f])
            min_outcrop_list <- c(min_outcrop_list, formations$min_Outcrop[f])
            form_counter <- form_counter + 1 
          }
        }
      }
    }
    Bin_Form_List2[[b]] <- unique(temp_form_list)
    Bin_Dino_List2[[b]] <- unique(temp_dino_list) # This bit kinda broken I think. Supposed to add all the unique dinos from a time bin to a final list.
    Max_Outcrop_list2[[b]] <- max_outcrop_list
    Min_Outcrop_list2[[b]] <- min_outcrop_list
    Formation_Counter_List2[[b]] <- form_counter
  }
  Area_M2_Max <- unlist(lapply(Max_Outcrop_list2,sum))
  Area_M2_Min <- unlist(lapply(Min_Outcrop_list2,sum))
  Species_M2 <- c() # METHOD 2
  for (f in 1:length(Bin_Dino_List2)) {
    Species_M2 <- c(Species_M2, length(Bin_Dino_List2[[f]]))
  }
  Species_M2 <<- Species_M2
  binpoints <- unlist(lapply(binlist, mean))

  par(mar=c(5, 4, 4, 6) + 0.1)
  plot(binpoints, Species_M2, type = 'l', yaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "black")
  axis(2, ylim=c(0,1),col="black",las=1)
  mtext("Raw Dinosaur Diversity",side=2,line=2.5)
  box()
  par(new=TRUE)
  plot(binpoints, Area_M2_Max, type = 'l', yaxt='n', xaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "red", lty=2)
  lines(binpoints, Area_M2_Min, col="blue", lty=2)
  mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
  axis(4, ylim=c(0,7000), col="black",col.axis="black",las=1)
  mtext("Age (Ma)",side=1,col="black",line=2.5) 
  legend("topleft",legend=c("Dinosaur Diversity","Max. Outcrop Area", "Min. Outcrop Area"),
  text.col=c("black","red","blue"),lty=c(1),col=c("black","red", "blue"),cex=0.7)

  DPSK_M2_Max <- Species_M2/Area_M2_Max
  DPSK_M2_Min <- Species_M2/Area_M2_Min
  
  par(mar=c(5, 4, 4, 6) + 0.1)
  plot(binpoints, DPSK_M2_Max, type = 'l', yaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "black")
  lines(binpoints, DPSK_M2_Min, col = "black")
  axis(2, ylim=c(0,1),col="black")
  mtext("DPSK",side=2,line=2.5)
  box()
  par(new=TRUE)
  plot(binpoints, Area_M2_Max, type = 'l', yaxt='n', xaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "red", lty=2)
  lines(binpoints, Area_M2_Min, col="blue", lty=2)
  mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
  axis(4, ylim=c(0,7000), col="black",col.axis="black",las=1)
  mtext("Age (Ma)",side=1,col="black",line=2.5) 
}


### METHOD 3: ASSIGN % OF DINOSAUR OCCS AT RANDOM TO A BIN BASED ON % OF FOPRMATION IN BIN - REPEAT AND AVERAGE ###

FormBin_M3<- function(formations, binlist, Form_list, times) {
  Bin_Form_List3 <- list()
  Bin_Dino_List3 <- list()
  Max_Outcrop_list3 <- list()
  Min_Outcrop_list3 <- list() 
  # make an empty list of the dinos in each bin
  Formation_Counter_List3 <- list() 

  times <- times
  av_species <- as.list(rep(0, length(binlist)))
  
  for(n in 1:times){
    for(b in 1:length(binlist)){ # for each new bin
      print(b)
      temp_form_list <- c()
      temp_dino_list <- c()
      min_outcrop_list <- c()
      max_outcrop_list <- c()
      form_counter <- 0
      for (f in 1:nrow(formations)) { # for each formation
        if (formations$max_age[f] > binlist[[b]][1] && formations$min_age[f] > binlist[[b]][1] | 
            binlist[[b]][2] > formations$max_age[f]  && binlist[[b]][2] > formations$min_age[f]) { # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
              cat("bin:", b, "\n") # Just print, don't do anything else.
        } 
        else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
          print(1)
          if (formations$max_age[f] <= binlist[[b]][1] && formations$min_age[f] >= binlist[[b]][2]){ # If formation sits within boundaries
            temp_form_list <- c(temp_form_list, as.character(formations$Formation[f])) 
            temp_dino_list <- c(temp_dino_list, as.character(Form_list[[f]][[1]])) # Add dinos from that formation to dino list
            max_outcrop_list <- c(max_outcrop_list, formations$max_Outcrop[f])
            min_outcrop_list <- c(min_outcrop_list, formations$min_Outcrop[f])
            form_counter <- form_counter + 1 
          }
          if (formations$max_age[f] >= binlist[[b]][1] && formations$min_age[f] >= binlist[[b]][2]){ #If formation crosses upper age limit
            print(2)
            x <- as.numeric(as.numeric(binlist[[b]][1]) - formations$min_age[f])
            y <- as.numeric(formations$max_age[f] - formations$min_age[f])         
            perc <- (x/y)
            ran_num <- perc*length(Form_list[[f]][[1]])
            sampled <- sample(Form_list[[f]][[1]], ran_num, replace = FALSE)
            temp_form_list <- c(temp_form_list, as.character(formations$Formation[f]))
            temp_dino_list <- c(temp_dino_list, as.character(sampled)) # Add dinos from that formation to dino list
            max_outcrop_list <- c(max_outcrop_list, perc*as.numeric(formations$max_Outcrop[f]))
            min_outcrop_list <- c(min_outcrop_list, perc*as.numeric(formations$min_Outcrop[f])) 
            form_counter <- form_counter + 1 
          } # This bit (above) should be done
          if (formations$max_age[f] <= binlist[[b]][1] && formations$min_age[f] <= binlist[[b]][2]){
            print(3)
            x <- as.numeric(formations$max_age[f] - as.numeric(binlist[[b]][2]))
            y <- as.numeric(formations$max_age[f] - formations$min_age[f])
            perc <- (x/y)
            ran_num <- perc*length(Form_list[[f]][[1]])
            sampled <- sample(Form_list[[f]][[1]], ran_num, replace = FALSE)
            temp_form_list <- c(temp_form_list, as.character(formations$Formation[f]))
            temp_dino_list <- c(temp_dino_list, as.character(sampled))# Add dinos from that formation to dino list
            max_outcrop_list <- c(max_outcrop_list, perc*as.numeric(formations$max_Outcrop[f]))
            min_outcrop_list <- c(min_outcrop_list, perc*as.numeric(formations$min_Outcrop[f]))
            form_counter <- form_counter + 1 
          } # This bit (above) should be done
          if (formations$max_age[f] >= binlist[[b]][1] && formations$min_age[f] <= binlist[[b]][2]){
            print(4)
            x <- as.numeric(formations$max_age[f] - formations$min_age[f])
            y <- as.numeric(as.numeric(binlist[[b]][1]) - as.numeric(binlist[[b]][2]))
            perc <- (y/x)
            ran_num <- perc*length(Form_list[[f]][[1]])
            sampled <- sample(Form_list[[f]][[1]], ran_num, replace = FALSE)
            temp_form_list <- c(temp_form_list, as.character(formations$Formation[f]))
            temp_dino_list <- c(temp_dino_list, as.character(sampled)) # Add dinos from that formation to dino list
            max_outcrop_list <- c(max_outcrop_list, perc*as.numeric(formations$max_Outcrop[f]))
            min_outcrop_list <- c(min_outcrop_list, perc*as.numeric(formations$min_Outcrop[f]))
            form_counter <- form_counter + 1    
          } # This bit (above) should be done
        }
      }
      Bin_Form_List3[[b]] <- unique(temp_form_list)
      Bin_Dino_List3[[b]] <- unique(temp_dino_list) 
      Max_Outcrop_list3[[b]] <- max_outcrop_list
      Min_Outcrop_list3[[b]] <- min_outcrop_list
      Formation_Counter_List3[[b]] <- form_counter
    }
    for (q in 1:length(Bin_Dino_List3)){
      av_species[[q]] <- c(av_species[[q]],length(Bin_Dino_List3[[q]])) # Find number of species per bin for this cycle
    }   
  }
  Area_M3_Max <- unlist(lapply(Max_Outcrop_list3,sum))
  Area_M3_Min <- unlist(lapply(Min_Outcrop_list3,sum))
  av_species <- lapply(av_species, function(x) x[-1])
  Species_M3 <- c() # METHOD 3
  for (t in 1:length(av_species)) {
    Species_M3 <- c(Species_M3, mean(av_species[[t]]))
  }
  Species_M3 <<- Species_M3
  binpoints <- unlist(lapply(binlist, mean))

  par(mar=c(5, 4, 4, 6) + 0.1)
  plot(binpoints, Species_M3, type = 'l', yaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "black")
  axis(2, ylim=c(0,1),col="black",las=1)
  mtext("Raw Dinosaur Diversity",side=2,line=2.5)
  box()
  par(new=TRUE)
  plot(binpoints, Area_M3_Max, type = 'l', yaxt='n', xaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "red", lty=2)
  lines(binpoints, Area_M3_Min, col="blue", lty=2)
  mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
  axis(4, ylim=c(0,7000), col="black",col.axis="black",las=1)
  mtext("Age (Ma)",side=1,col="black",line=2.5) 
  legend("topleft",legend=c("Dinosaur Diversity","Max. Outcrop Area", "Min. Outcrop Area"),
  text.col=c("black","red","blue"),lty=c(1),col=c("black","red", "blue"),cex=0.7)

  DPSK_M3_Max <- Species_M3/Area_M3_Max
  DPSK_M3_Min <- Species_M3/Area_M3_Min
  
  par(mar=c(5, 4, 4, 6) + 0.1)
  plot(binpoints, DPSK_M3_Max, type = 'l', yaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "black")
  lines(binpoints, DPSK_M3_Min, col = "black")
  axis(2, ylim=c(0,1),col="black")
  mtext("DPSK",side=2,line=2.5)
  box()
  par(new=TRUE)
  plot(binpoints, Area_M3_Max, type = 'l', yaxt='n', xaxt='n', xlim=rev(range.default(binpoints)), xlab="", ylab="", col = "red", lty=2)
  lines(binpoints, Area_M3_Min, col="blue", lty=2)
  mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
  axis(4, ylim=c(0,7000), col="black",col.axis="black",las=1)
  mtext("Age (Ma)",side=1,col="black",line=2.5) 
}