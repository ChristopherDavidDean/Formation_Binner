#=============================================== FUNCTIONS FOR FORMATION BINNER =================================================
#
#      Selection of functions that work with PBDB and Formation data to set up Formation Binning. Functions include those
#      that visualise or retrieve data, those that setup and run Formation Binning, and those which run diversity analyses
#      on inputted data with varying bin types. 
#
#=============================================== iPAK AND REQUIRED PACKAGES =====================================================

# function that automatically installs necessary packages that the user is lacking.
#===== iPAK =====
ipak <- function(pkg){ # Function to install packages. Read in character vector of any packages required. 
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

#=============================================== GET_INFO ============================================================

# Retrieves basic information about inputted formations, including mean age and range.
GetInfo <- function(formations){ # Require formation data
  mean <- rowMeans(subset(formations, select = c(3, 4)), na.rm = TRUE)
  range <- formations$max_age - formations$min_age
  par(mfrow=c(2,1))
  hist(mean, main = "Mean Age of Formations",xlab = "Mean Age (Ma)")
  hist(range, main = "Range of Formations (Ma)", xlab = "Range (Ma)")
  cat ('Mean range of formations: ', mean(range), fill = TRUE)
  cat ('Median range of formations: ', median(range),fill = TRUE)
}

#=============================================== SCORING_GRID ========================================================

# Functions to create a score grid to identify best position to draw boundaries between formations. The function 
# looks through time at intervals and gives a suitability score for each formation to draw a boundary at that point.
# If a formation does not cross that boundary, the formation is given a score of 100. If it does cross that boundary,
# it works out how much of the formation sits each side of the boundary, and downweights accordingly - e.g. if a 
# formation was found to have 10% of it's total range one side of a boundary and 90% the other, it would receive a
# higher score than a formation where 50% of it's range sat either side. Mean scores are generated for each time bin.
      
#==== Scoring_Grid_1 ====

# Creates a scoring grid using info from all formations.
Scoring_Grid_1 <- function(formations, res=0.01) { # Requires formation information. Resolution of time lines is set automatically at 0.01, but can be adjusted.
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  allbins <- seq(min_age-1.0045, max_age+1.0045, res) # Makes 1ma bins in sequence based on max/min ages. 0.0045 added to ensure formation is never exactly equivalent to a bin.
  
  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(allbins)) # makes a matrix for the scoring 
  # of each time line in terms of how good it is to be a bin boundary
  colnames(score_grid) <- allbins # All time lines
  rownames(score_grid) <- formations$Formation # All formations
  
  counter <- 0
  for(i in allbins) { # Go through each time line
    counter <- sum(counter,1) 
    for (f in 1:nrow(formations)){ # go through each formation 
      if (i <= formations$max_age[f] && i >= formations$min_age[f]){ # if timeline is between max/min age of a formation (i.e. formation crosses that line)
        a <- formations$max_age[f] - i # Work out how much of formation is older than timeline
        b <- i - formations$min_age[f] # Work out how much of formation is younger than timeline
        range <- formations$max_age[f] - formations$min_age[f] # Calculate range of formation
        if (a > b){
          score_grid[,counter][f] <- (a/range)*100 # Work out percentage that sits each side of line, reduce score by that amount.
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
  score_grid <<- score_grid # Outputs score_grid for later use
  allbins <<- allbins # Outputs bins for later use
}

#==== Scoring_Grid_2 ====

# Create a scoring grid ignoring formations with length longer than mean formation length. In his way, long ranging formations don't 
# bias the creation of bins, especially when they appear during the same time interval.
Scoring_Grid_2 <- function(formations, res=0.01) { # Requires formation information. Resolution of time lines is set automatically at 0.01, but can be adjusted.
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  allbins <- seq(min_age-1.0045, max_age+1.0045, res) # Makes 1ma bins in sequence based on max/min ages. 0.0001 added to ensure formation is never exactly equivalent to a bin.
  
  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(allbins)) # makes a matrix for the scoring of each time line in terms of how good it is to be a bin boundary
  colnames(score_grid) <- allbins # All time lines
  rownames(score_grid) <- formations$Formation # All formations
  
  counter <- 0
  for(i in allbins) { # Go through each time line
    counter <- sum(counter,1)
    for (f in 1:nrow(formations)){ # go through each formation 
      if (formations$max_age[f] - formations$min_age[f] < # If formation range is less than the mean formation range
          mean(formations$max_age - formations$min_age)) {
        if (i <= formations$max_age[f] && i >= formations$min_age[f]){ # if timeline is between max/min age of a formation (i.e. formation crosses that line)
          a <- formations$max_age[f] - i # Work out how much of formation is older than timeline
          b <- i - formations$min_age[f] # Work out how much of formation is younger than timeline
          range <- formations$max_age[f] - formations$min_age[f] # Calculate range of formation
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
      else { # If formation range is longer than mean formation range, skip (bin drawing isn't affected)
        next
      }
    }
  }  
  score_grid <- na.omit(score_grid) # Remove effect of formations longer than mean formation range
  means <- colMeans(score_grid) # Work out mean score for each time bin
  score_grid <- rbind(score_grid, means) # add to grid
  score_grid <<- score_grid
  allbins <<- allbins
}

#=============================================== PLOTMAKER =======================================================

# Generates plots through time with user inputted data and Formation_Bins, whilst providing traditional stage 
# data for comparison. NOTE: xlim is specified to fit the chosen time window of this study - as such, this would
# have to be adjusted if other data were to be used.
plotMaker <- function(rel_data, binlist, ulabel){ # Takes relevant data for plotting, binlist(dataframe of bins) and a user generated level for y axis. 
  useful_bins <- c(binlist$bottom, binlist$top[nrow(binlist)]) # Converts binlist into format useful for plotting
  tsplot(stages, boxes=c("short","system"), ylab = ulabel, # Creates plot using data from DivDyn package. 
         xlim=75:81,  ylim=c(0,(max(rel_data, na.rm = TRUE)+(max(rel_data, na.rm = TRUE)*0.1))), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  
  for(n in 1:length(useful_bins)){ # For each bin
    if(((n %% 2) == 0) == TRUE) next # If the bin is even, skip it (allows for alternating colours of bins)
    else {
      if(n == length(useful_bins)){
        if(nrow(binlist) %% 2 == 0){ # Skips colouring last bin if there are an even number of bins
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          (max(rel_data, na.rm = TRUE)+(max(rel_data, na.rm = TRUE)*0.1)), 
          col = "#32323232", border = NA)
        }
      }
      else{
        rect(useful_bins[n], 0, useful_bins[n+1], 
        (max(rel_data, na.rm = TRUE)+(max(rel_data, na.rm = TRUE)*0.1)),
        col = "#32323232", border = NA)
      }
    }
  }
  lines(binlist$mid, rel_data, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1) # Adds lines from relevant data provided.
}

#=============================================== OVERLAP_COUNTER ==============================================================

# Generates a quick high resolution plot of the number of formations present through geologic time. 

overlap_counter <- function(score_grid){ # Takes score_grid as input, generated from either Scoring_Grid_1 or Scoring_Grid_2 functions. 
  score_grid[score_grid < 100] <- 1 # Turn all score less than 100 (i.e. formation exists at that point) into 1's
  score_grid[score_grid == 100] <- 0 # Turn all scores of 100 (no formation present) into 0's
  total_forms <- colSums(score_grid) # Sum the number of formations occurring in each bin
  par(mar = c(4.1, 4.1, 1, 2.1))
   tsplot(stages, boxes=c("short","system"), # Generates plot using DivDyn package
         xlim=75:81,  ylim=c(min(colMeans(score_grid), na.rm = TRUE), max(total_forms)+1), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Number of Formations") 
  lines(allbins, total_forms)
}

#=============================================== NEWBINS ==============================================================

# Looks at the previously generated score_grid and generates appropriate new bins based on those scores. Boundaries are 
# outputted as a list (binlist). If bins are shorter than 0.5 Ma, they are amalgamated into the bins above and below, and 
# a warning is produced.

newBins <- function(score_grid, formations, bin_limits, allbins, stages, smallamalg = TRUE){ # Takes previously generated score grid, formations, 
  # allbins and stages from DivDyn package. Also require bin_limits, a user made vector of the following: 
  # 1) user chosen time window in which to look to draw bins. Advised to be set at 3 Ma. 
  # 2) Hard maximum age of bins
  # 3) Hard minimum age of bins

  score_grid<- as.data.frame(score_grid)
  bin_size <- bin_limits[1]
  max_age <- bin_limits[2]
  min_age <- bin_limits[3]
  form_bins <- c(min_age)
  testbin <- seq(min_age, max_age, bin_size) # Creates broader bins for testing best point to draw a bin
  
  # Drawing bins and giving form_bins (vector of bin boundaries)
  for (i in 1:length((testbin)-1)){
    seqs <- seq(testbin[i],testbin[i]+bin_size, 1) # Creates a sequence of ages to draw bins within
    pasting <- c()
    for (n in 1:bin_size){
      pasting <- c(pasting, paste("^",seqs[n],"|", sep = "")) # Sets up expression to match to score_grid
    }
    testmatch2 <- paste(pasting, collapse = "")
    testmatch2 <- substr(testmatch2, 1, nchar(testmatch2)-1) 
    a <- score_grid[grep(testmatch2, names(score_grid))] # Finds all scores within this time window
    z <- apply(a,1,which.max) # Finds maximum bin score within this time window
    form_bins[i+1] <- names(a)[z][nrow(a)] # Adds maximum bin score to a vector of bins
  }
  form_bins <- as.numeric(unique(form_bins)) # Finds all unique bins
  if (smallamalg == TRUE){ # If small bin amalgamation is turned on (is on automatically):
    range <- (diff(form_bins) < 0.5) # Finds all bins which are under 0.5 Ma in length
    range_checker <- c()
    for (r in 1:length(range)){
      if (range[r] == TRUE){ # If a bin is under 0.5 Ma in length
        difference <- diff(c(form_bins[r], form_bins[r+1])) # Find the length of that bin
        warning("Original bin ",  r, " removed due to small range: ~", signif(difference, digits = 3), " Ma. The difference in time has been added to the bins above and below.") # Generate warning about bin amalgamation
        form_bins[r] <- form_bins[r]+(difference/2) # Adds half length of old bin to bin below
        form_bins[r+1] <- form_bins[r+1]-(difference/2) # Adds half length of old bin to bin above
        range_checker <- c(range_checker, r) # Records which bin was too small
      }
    }
    if(length(range_checker) > 0){ # If there have been amalgamated bins
      form_bins <- form_bins[-range_checker] # Remove old amalgamated bins
    }
  }
  if (smallamalg == FALSE){
    warning("Small bin amalgamation is turned off. Bins may be too short to record occurrences. You are advised to check bins before running further analyses.") # Generate warning if small bin amalgamation is turned off. 
  }
  form_bins <<- form_bins
  
  # Creating binlist (data.frame of bins and appropriate age info)
  prefix <- "FB."
  suffix <- seq(1:(length(form_bins)-1))
  my_names <- paste(prefix, suffix, sep = "")
  binlist <- data.frame(bin = my_names, # Combines bin data to make dataframe of minimum, maximum and mid point of each new bin
                        bottom = as.numeric(form_bins[1:(length(form_bins)-1)]), 
                        top = as.numeric(form_bins[2:(length(form_bins))]))
  binlist$mid <- (binlist$bottom + binlist$top) / 2
  binlist <<- binlist

  #Plotting new bins
  par(mar = c(4.1, 4.1, 1, 2.1))
  tsplot(stages, boxes=c("short","system"), # Generates plot using DivDyn package
         xlim=75:81,  ylim=c(min(colMeans(score_grid), na.rm = TRUE), 100), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Bin Splitting Score") 
  lines(allbins, colMeans(score_grid)) # draws bin splitting score on plot
  
  
  for(n in 1:length(form_bins)){ # draws new bins as coloured boxes for comparison to traditional bins
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(form_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 100, col = "#32323232", border = NA)
        }
      }
      else{
        rect(form_bins[n], 0, form_bins[n+1], 100-0.01, col = "#32323232", border = NA)
      }
    }
  }
}

#=============================================== FORMATIONGRAPH =======================================================

# Shows what formations look like through time in comparison to Stages and new Bins.
FormationGraph <- function(formations, form_bins, stages){ # Requires formations, form_bins from newBins function and stages from DivDyn package
  fp <- data.frame("x1" = formations$min_age, "y1" = formations$forbinning, 'x2' = formations$max_age, 'y2' = formations$forbinning)
  par(mar = c(4.1, 4.1, 1, 2.1))
  tsplot(stages, boxes=c("short","system"), # Generates plot using Divdyn package
         xlim=75:81,  ylim=range(fp$y1, fp$y2), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Formations by Number")
  segments(fp$x1, fp$y1, fp$x2, fp$y2, lwd = 2) # Plots formations as lines showing their duration
  
  for(n in 1:length(form_bins)){ # draws new bins as coloured boxes for comparison to traditional bins
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(form_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          max(fp$y2, na.rm = TRUE), col = "#32323232", border = NA)
        }
      }
      else{
        rect(form_bins[n], 0, form_bins[n+1], 
        max(fp$y2, na.rm = TRUE), col = "#32323232", border = NA)
      }
    }
  }
}

#=============================================== DIVERSITY METHODS ====================================================

# A variety of methods which generate diversity curves based on the generated Formation Bins. 

#==== FormBin_M1 ====
# Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to all bins
# that they occur in. Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen
# Quorum levels. 

FormBin_M1 <- function(formations, binlist, Form_list, Quorum) {
  sqsmst <- list()
  for (q in 1:length(Quorum)){
    M1_List <- list()# make an empty list of the occurrences in each bin

    # Code for assigning formations and associated occurrences to bins
    for(b in 1:nrow(binlist)){ # for each new bin
      temp_recs <- data.frame()
      for (f in 1:nrow(formations)){ 
        if (formations$max_age[f] >= binlist[b,2] && formations$min_age[f] <= binlist[b,3]){ # If Formation max. age is greater than Bin min. age AND if Formation min. age is less then Bin max. age (i.e. falls within bin at some point)
          temp_recs <- rbind(temp_recs, Form_list[[f]]) # Add occurrences from that formation to occurrence list
        }
      }
      if (nrow(temp_recs) > 0){
        temp_recs$bin_no <- b
      }
      M1_List[[b]] <- temp_recs
    }
    df <- do.call("rbind", M1_List)

    # Code for subsampling
    niter <- 100 # number of SQS iterations
    SQS <- subsample(df,iter=100, q=Quorum[q],tax="occurrence.genus_name", bin="bin_no", 
                    coll = 'collection_no', output="dist", type="sqs", 
                    duplicates = TRUE, useFailed = TRUE)
    
    #calculate Mean, SD and 95% confidence intervals for SQS data
    sqsaverage <- rowMeans2(SQS$divSIB, na.rm = TRUE)
    sqssd <- rowSds(SQS$divSIB, na.rm = TRUE)
    sqsminsd <- sqsaverage - sqssd
    sqsmaxsd <- sqsaverage + sqssd
  
    #bind data
    combined.sqs <- cbind(sqsaverage, sqssd, sqsminsd, sqsmaxsd)
    rownames(combined.sqs) <- binlist$bin
    combined.sqs <- as.data.frame(combined.sqs)
    
    temp_name <- paste("q.",deparse(Quorum[q]),"_", "SQS_Results", sep = "") #Name files based on data entered to function
    assign(temp_name, combined.sqs, envir = .GlobalEnv)
    sqsmst[[q]] <- combined.sqs
    names(sqsmst)[[q]] <- deparse(Quorum[q])
  }
  sqsmst <<- sqsmst
  
  #Binstats
  bin_info <- binstat(df, tax="occurrence.genus_name", bin="bin_no", 
                      coll = 'collection_no')
  bin_info[bin_info==0] <- NA
  bin_info <<- bin_info
  
  # Recording Results
  dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
  write.csv(bin_info, file.path(paste("Results/M1_Bin_info.csv", sep="")))
  for (q in 1:length(Quorum)){
    temp_name <- paste("M1_SQS_", Quorum[q], sep = "")
    write.csv(sqsmst[q], file.path(paste("Results/", temp_name, ".csv", sep="")))
  }

  # Plotting Raw Div and Sampling Proxies
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
  par(mar = c(0, 4.1, 4.1, 2.1))
  with(binlist, tsplot(stages, ylab = "Raw Diversity",
              xlim=75:81,  ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
              shading=NULL, plot.args=list(xaxt = 'n')))
  useful_bins <- c(binlist$bottom, binlist$top[nrow(binlist)])
  for(n in 1:length(useful_bins)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(useful_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          (max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1)), 
          col = "#32323232", border = NA)
        }
      }
      else{
        rect(useful_bins[n], 0, useful_bins[n+1], 
        (max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1)), 
        col = "#32323232", border = NA)
      }
    }
  }
  lines(binlist$mid, bin_info$SIBs, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)
  par(mar = c(4.1, 4.1, 0, 2.1))
  with(plotMaker(bin_info$colls, binlist, "Number of Collections"),0)
  
  par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
  plotMaker(bin_info$u, binlist, "Good's u")
  
  # Plotting SQS
  tsplot(stages, boxes=c("short","system"), 
         xlim=75:81,  ylim=c(0,(max(sqsmst[length(sqsmst)][[1]][1])+(max(sqsmst[length(sqsmst)][[1]][1])*0.1))), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Subsampled Diversity")  
  for (q in 1:length(Quorum)){
    g.col <- gray.colors(length(Quorum), start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL)
    lines(binlist$mid, sqsmst[[q]]$sqsaverage, type = 'o', col = g.col[q], 
          pch = 21, bg = "grey")
    polygon(x = c(binlist$mid, rev(binlist$mid)),
            y = c(sqsmst[[q]]$sqsmaxsd, rev(sqsmst[[q]]$sqsminsd)), 
            col =  adjustcolor(g.col[q], alpha.f = 0.40), border = NA)
  }
  
  for(n in 1:length(form_bins)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(form_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1)), 
          col = "#32323232", border = NA)
        }
      }
      else{
        rect(form_bins[n], 0, form_bins[n+1], 
        (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1)), 
        col = "#32323232", border = NA)
      }
    }
  }
}

#==== FormBin_M2 ====
# Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to bins that
# the majority of the formation occurs in. Formations which have lengths greater than 2x the length of a bin are ignored.
# Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen Quorum levels. 

FormBin_M2<- function(formations, binlist, Form_list, Quorum) {
  sqsmst <- list()
  for (q in 1:length(Quorum)){
    M2_List <- list()# make an empty list of the occurrencess in each bin
    for(b in 1:nrow(binlist)){ # for each new bin
      temp_recs <- data.frame()
      for (f in 1:nrow(formations)) { # for each formation
        if (formations$max_age[f] < binlist[b,2] |
            binlist[b,3] < formations$min_age[f]) { # If the formation DOES NOT sit in this bin 
          next # Just print, don't do anything else.
        }
        else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
          if (formations$max_age[f] <= binlist[b,3] && 
              formations$min_age[f] >= binlist[b,2]){# If formation sits within boundaries
            temp_recs <- rbind(temp_recs, Form_list[[f]]) # Add occurrences from that formation to occurrence list
            next
          }
          if (formations$max_age[f] >= binlist[b,2] && 
              formations$min_age[f] <= binlist[b,2] && 
              formations$max_age[f] <= binlist[b,3]){ #If formation crosses upper age limit (minimum age)
            x <- as.numeric(formations$max_age[f] - as.numeric(binlist[b,2]))
            y <- as.numeric(as.numeric(binlist[b,2]) - formations$min_age[f])
            if (x > y){
              temp_recs <- rbind(temp_recs, Form_list[[f]]) # Add occurrences from that formation to occurrence list
              next
            }
          }
          if (formations$max_age[f] >= binlist[b,3] && 
              formations$min_age[f] <= binlist[b,3] &&
              formations$min_age[f] >= binlist[b,2]){ #If formation crosses lower age limit (maximum age)
            x <- as.numeric(formations$max_age[f] - as.numeric(binlist[b,3]))
            y <- as.numeric(as.numeric(binlist[b,3]) - formations$min_age[f])
            if (y > x){
              temp_recs <- rbind(temp_recs, Form_list[[f]]) # Add occurrences from that formation to occurrence list
              next
            }
          }
          if (formations$max_age[f] > binlist[b,3] &&
              formations$min_age[f] < binlist[b,2]){
            temp_recs <- rbind(temp_recs, Form_list[[f]]) # Add occurrences from that formation to occurrence list
            next
          }
        }
      }
      if (nrow(temp_recs) > 0){
        temp_recs$bin_no <- b
      }
      M2_List[[b]] <- temp_recs
    }
    df <- do.call("rbind", M2_List)
    
    # Code for subsampling
    niter <- 100 # number of SQS iterations
    SQS <- subsample(df,iter=100, q=Quorum[q],tax="occurrence.genus_name", bin="bin_no", 
                       coll = 'collection_no', output="dist", type="sqs", 
                       duplicates = TRUE, useFailed = TRUE)
    SQS$divSIB[SQS$divSIB == 0] <- NA
    while (nrow(SQS$divSIB) < nrow(binlist) ){ # If last bin has failed
      SQS$divSIB <- rbind(SQS$divSIB, rep(NA, 100)) # Add additional row of NA's for that bin. Temporary solution - need to match to ensure that 
    }
    
    #calculate Mean, SD and 95% confidence intervals for SQS data
    sqsaverage <- rowMeans2(SQS$divSIB, na.rm = TRUE)
    sqssd <- rowSds(SQS$divSIB, na.rm = TRUE)
    sqsminsd <- sqsaverage - sqssd
    sqsmaxsd <- sqsaverage + sqssd
  
    #bind data
    combined.sqs <- cbind(sqsaverage, sqssd, sqsminsd, sqsmaxsd)
    combined.sqs[is.nan(combined.sqs)] <- NA
    rownames(combined.sqs) <- binlist$bin
    combined.sqs <- as.data.frame(combined.sqs)
    
    temp_name <- paste("q.",deparse(Quorum[q]),"_", "SQS_Results", sep = "") #Name files based on data entered to function
    assign(temp_name, combined.sqs, envir = .GlobalEnv)
    sqsmst[[q]] <- combined.sqs
    names(sqsmst)[[q]] <- deparse(Quorum[q])
  }
  sqsmst <<- sqsmst
  
  #Binstats
  bin_info <- binstat(df, tax="occurrence.genus_name", bin="bin_no", 
                      coll = 'collection_no')
  while (nrow(bin_info) < nrow(binlist) ){ # If last bin has failed
      bin_info[nrow(bin_info)+1, ] <- NA # Add additional row of NA's for that bin. Temporary solution - need to match to ensure that 
      bin_info[nrow(bin_info), 1] <- nrow(bin_info)
  }
  bin_info[bin_info==0] <- NA
  bin_info <<- bin_info
  
  # Recording Results
  dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
  write.csv(bin_info, file.path(paste("Results/M2_Bin_info.csv", sep="")))
  for (q in 1:length(Quorum)){
    temp_name <- paste("M2_SQS_", Quorum[q], sep = "")
    write.csv(sqsmst[q], file.path(paste("Results/", temp_name, ".csv", sep="")))
  }
  
  # Plotting Raw Div and Sampling Proxies
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
  par(mar = c(0, 4.1, 4.1, 2.1))
  with(binlist, tsplot(stages, ylab = "Raw Diversity",
              xlim=75:81,  ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
              shading=NULL, plot.args=list(xaxt = 'n')))
  useful_bins <- c(binlist$bottom, binlist$top[nrow(binlist)])
  for(n in 1:length(useful_bins)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(useful_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          (max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1)), 
          col = "#32323232", border = NA)
        }
      }
      else{
        rect(useful_bins[n], 0, useful_bins[n+1], 
        (max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1)), 
        col = "#32323232", border = NA)
      }
    }
  }
  lines(binlist$mid, bin_info$SIBs, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)
  par(mar = c(4.1, 4.1, 0, 2.1))
  with(plotMaker(bin_info$colls, binlist, "Number of Collections"),0)
  
  par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
  plotMaker(bin_info$u, binlist, "Good's u")
  
  # Plotting SQS
  tsplot(stages, boxes=c("short","system"), 
         xlim=75:81,  ylim=c(0,(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+
                                  (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1))), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Subsampled Diversity")  
  for (q in 1:length(Quorum)){
    g.col <- gray.colors(length(Quorum), start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL)
    enc <- rle(!is.na(sqsmst[[q]]$sqsaverage))
    endIdxs <- cumsum(enc$lengths)
    for(i in 1:length(enc$lengths)){
      if(enc$values[i]){
        endIdx <- endIdxs[i]
        startIdx <- endIdx - enc$lengths[i] + 1
    
        subdat <- binlist$mid[startIdx:endIdx]
        submin <- sqsmst[[q]]$sqsminsd[startIdx:endIdx]
        submax <- sqsmst[[q]]$sqsmaxsd[startIdx:endIdx]
        subdepth <- sqsmst[[q]]$sqsaverage[startIdx:endIdx]
    
        x <- c(subdat, rev(subdat))
        y <- c(submax, rev(submin))
    
        polygon(x = x, y = y, col = adjustcolor(g.col[q], alpha.f = 0.40), border = NA)
        lines(binlist$mid, sqsmst[[q]]$sqsaverage, type = 'o', col = g.col[q], 
          pch = 21, bg = "grey")
      }
    }
  }
  for(n in 1:length(form_bins)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(form_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1)), 
          col = "#32323232", border = NA)
        }
      }
      else{
        rect(form_bins[n], 0, form_bins[n+1], 
        (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1)), 
        col = "#32323232", border = NA)
      }
    }
  }
}


#==== FormBin_M3 ====
# Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to all bins
# that they occur in, based on the percentage of the formation that sits within that bin. Occurrences are selected at
# random from the formation list and not replaced. The test is repeated according to user specified number of runs.
# Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen Quorum levels. 

FormBin_M3<- function(formations, binlist, Form_list, times=10, Quorum) {
  ptm <- proc.time()
  for (q in 1:length(Quorum)){
    allSQS <- data.frame(binlist$bin)
    allbininfo <- data.frame(binlist$bin)
    for(n in 1:times){
      M3_List <- list()
      for(b in 1:nrow(binlist)){ # for each new bin
        temp_recs <- data.frame() 
        for (f in 1:nrow(formations)) { # for each formation
          if (formations$max_age[f] < binlist[b,2] |
              binlist[b,3] < formations$min_age[f]) { # If the formation DOES NOT sit in this bin 
            next # Just print, don't do anything else.
          }
          else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
            if (formations$max_age[f] <= binlist[b,3] && 
              formations$min_age[f] >= binlist[b,2]){ # If formation sits within boundaries
              temp_recs <- rbind(temp_recs, Form_list[[f]]) # Add occurrences from that formation to occurrence list
              next
            }
            if (formations$max_age[f] >= binlist[b,2] && 
                formations$min_age[f] <= binlist[b,2] && 
                formations$max_age[f] <= binlist[b,3]){ #If formation crosses upper age limit (minimum age)
              x <- as.numeric(formations$max_age[f] - as.numeric(binlist[b,2])) # Inside bin
              y <- as.numeric(formations$max_age[f] - formations$min_age[f]) # Outside bin
              perc <- (x/y) # Percentage of formation inside bin
              ran_num <- perc*nrow(Form_list[[f]])
              temp_dframe <- Form_list[[f]]
              sampled <- temp_dframe[sample(nrow(temp_dframe), ran_num, replace = FALSE),]
              temp_recs <- rbind(temp_recs, sampled) # Add occurrences from that formation to occurrence list
              next
            } 
            if (formations$max_age[f] >= binlist[b,3] && 
                formations$min_age[f] <= binlist[b,3] &&
                formations$min_age[f] >= binlist[b,2]){ #If formation crosses lower age limit (maximum age)
              x <- as.numeric(as.numeric(binlist[b,3]) - formations$min_age[f]) # Inside bin
              y <- as.numeric(formations$max_age[f] - formations$min_age[f]) # Length of formation
              perc <- (x/y) # Percentage of formation inside bin
              ran_num <- perc*nrow(Form_list[[f]])
              temp_dframe <- Form_list[[f]]
              sampled <- temp_dframe[sample(nrow(temp_dframe), ran_num, replace = FALSE),]
              temp_recs <- rbind(temp_recs, sampled) # Add occurrences from that formation to occurrence list
              next
            }
            if (formations$max_age[f] >= binlist[b,3] &&
                formations$min_age[f] <= binlist[b,2]){
              x <- as.numeric(as.numeric(binlist[b,3]) - as.numeric(binlist[b,2])) # Length of bin
              y <- as.numeric(formations$max_age[f] - formations$min_age[f]) # Length of formation
              perc <- (x/y)
              ran_num <- perc*nrow(Form_list[[f]])
              temp_dframe <- Form_list[[f]]
              sampled <- temp_dframe[sample(nrow(temp_dframe), ceiling(ran_num), replace = FALSE),]
              temp_recs <- rbind(temp_recs, sampled) # Add occurrences from that formation to occurrence list
              next
            } 
          }
        }
        if (nrow(temp_recs) > 0){ # If there are occurrences in this bin
          temp_recs$bin_no <- b # Label those occurrences with bin number
        }
        M3_List[[b]] <- temp_recs # Add temp records to permanent list
      }
      df <- do.call("rbind", M3_List) # Create data.frame from list
        
      # Code for subsampling, and skipping errors with little data - Not very useful here, but useful for last method. Kept in just in case.
      niter <- 100 # number of SQS iterations
      SQS_test <- try(subsample(df,iter=niter, q=Quorum[q],tax="occurrence.genus_name", bin="bin_no", 
                         coll = 'collection_no', output="dist", type="sqs", 
                         duplicates = TRUE, useFailed = TRUE), silent = TRUE)
      if(class(SQS_test)=='try-error'){
        sqsaverage <- rep(NA, nrow(binlist))
        allSQS <- cbind.data.frame(allSQS, sqsaverage)
        next
      }
      else{
        SQS <- subsample(df,iter=100, q=Quorum[q],tax="occurrence.genus_name", bin="bin_no", 
                           coll = 'collection_no', output="dist", type="sqs", 
                           duplicates = TRUE, useFailed = TRUE)
      }
  
      SQS$divSIB[SQS$divSIB == 0] <- NA
      while (nrow(SQS$divSIB) < nrow(binlist) ){ # If last bin has failed
        SQS$divSIB <- rbind(SQS$divSIB, rep(NA, 100)) # Add additional row of NA's for that bin. Temporary solution - need to match to ensure that 
      }
      sqsaverage <- rowMeans2(SQS$divSIB, na.rm = TRUE)
      
      allSQS <- cbind(allSQS, sqsaverage)
      
      #Binstats
      bin_info <- binstat(df, tax="occurrence.genus_name", bin="bin_no", 
                            coll = 'collection_no')
      while (nrow(bin_info) < nrow(binlist) ){ # If last bin has failed
        bin_info[nrow(bin_info)+1, ] <- NA # Add additional row of NA's for that bin. Temporary solution - need to match to ensure that 
        bin_info[nrow(bin_info), 1] <- nrow(bin_info)
      }
      bin_info[bin_info==0] <- NA
      
      if (n == 1){
        M3_Bin_Info <- list(bin_info$occs, bin_info$colls, bin_info$SIBs, bin_info$u)
      }
      else{
        M3_Bin_Info[[1]] <- cbind(M3_Bin_Info[[1]], bin_info$occs)
        M3_Bin_Info[[2]] <- cbind(M3_Bin_Info[[2]], bin_info$colls)
        M3_Bin_Info[[3]] <- cbind(M3_Bin_Info[[3]], bin_info$SIBs)
        M3_Bin_Info[[4]] <- cbind(M3_Bin_Info[[4]], bin_info$u)
      }
    }
    #calculate Mean, SD and 95% confidence intervals for SQS data
    allSQS <- allSQS %>% remove_rownames %>% column_to_rownames(var="binlist.bin")
    allsqsaverage <- rowMeans2(as.matrix(allSQS), na.rm = TRUE)
    allsqssd <- rowSds(as.matrix(allSQS), na.rm = TRUE)
    nr <- min(rowSums(is.na(allSQS))) #number of unsuccessful runs
    nr <- times - nr # number of successful runs
    sqserror <- qnorm(0.975) * allsqssd/sqrt(nr) #2.5% each side of tail
    sqsminerror <- allsqsaverage - sqserror
    sqsmaxerror <- allsqsaverage + sqserror
  
    #bind data
    combined.sqs <- cbind(allsqsaverage, allsqssd, nr, sqserror, sqsminerror, sqsmaxerror)
    rownames(combined.sqs) <- binlist$bin
    combined.sqs <- as.data.frame(combined.sqs)
    
    temp_name <- paste("q.",deparse(Quorum[q]),"_", "SQS_Results", sep = "") #Name files based on data entered to function
    assign(temp_name, combined.sqs, envir = .GlobalEnv)
    sqsmst[[q]] <- combined.sqs
    names(sqsmst)[[q]] <- deparse(Quorum[q])
    
    # Calculate Averaged Bin_Info
    occsaveraged <- rowMeans2(M3_Bin_Info[[1]])
    collsaveraged <- rowMeans2(M3_Bin_Info[[2]])
    divaveraged <- rowMeans2(M3_Bin_Info[[3]])
    goodsuaveraged <- rowMeans2(M3_Bin_Info[[4]])
    combined.bin_info <- cbind(occsaveraged, collsaveraged, divaveraged, goodsuaveraged)
    rownames(combined.bin_info) <- binlist$bin
    combined.bin_info <- as.data.frame(combined.bin_info)
    combined.bin_info <<- combined.bin_info
  }

  # Recording Results
  dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
  write.csv(combined.bin_info, file.path(paste("Results/M3_Bin_info.csv", sep="")))
  for (q in 1:length(Quorum)){
    temp_name <- paste("M3_SQS_", Quorum[q], sep = "")
    write.csv(sqsmst[q], file.path(paste("Results/", temp_name, ".csv", sep="")))
  }
  
  # Plotting Raw Div and Sampling Proxies - IN PROGRESS
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
  par(mar = c(0, 4.1, 4.1, 2.1))
  with(binlist, tsplot(stages, ylab = "Raw Diversity",
              xlim=75:81,  ylim=c(0,(max(combined.bin_info[["divaveraged"]], na.rm = TRUE)
              +(max(combined.bin_info[["divaveraged"]], na.rm = TRUE)*0.1))), 
              shading=NULL, plot.args=list(xaxt = 'n')))
  useful_bins <- c(binlist$bottom, binlist$top[nrow(binlist)])
  for(n in 1:length(useful_bins)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(useful_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          (max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1)), 
          col = "#32323232", border = NA)
        }
      }
      else{
        rect(useful_bins[n], 0, useful_bins[n+1], 
        (max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1)), 
        col = "#32323232", border = NA)
      }
    }
  }
  lines(binlist$mid, combined.bin_info$divaveraged, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)
  par(mar = c(4.1, 4.1, 0, 2.1))
  with(plotMaker(combined.bin_info$collsaveraged, binlist, "Number of Collections"),0)
  
  par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
  plotMaker(combined.bin_info$goodsuaveraged, binlist, "Good's u")
  
  # Plotting SQS
  tsplot(stages, boxes=c("short","system"), 
         xlim=75:81,  ylim=c(0,(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+
                                  (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1))), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Subsampled Diversity")  
  for (q in 1:length(Quorum)){
    g.col <- gray.colors(length(Quorum), start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL)
    enc <- rle(!is.na(sqsmst[[q]]$allsqsaverage))
    endIdxs <- cumsum(enc$lengths)
    for(i in 1:length(enc$lengths)){
      if(enc$values[i]){
        endIdx <- endIdxs[i]
        startIdx <- endIdx - enc$lengths[i] + 1
    
        subdat <- binlist$mid[startIdx:endIdx]
        submin <- sqsmst[[q]]$sqsminerror[startIdx:endIdx]
        submax <- sqsmst[[q]]$sqsmaxerror[startIdx:endIdx]
        subdepth <- sqsmst[[q]]$allsqsaverage[startIdx:endIdx]
    
        x <- c(subdat, rev(subdat))
        y <- c(submax, rev(submin))
    
        polygon(x = x, y = y, col = adjustcolor(g.col[q], alpha.f = 0.40), border = NA)
        lines(binlist$mid, sqsmst[[q]]$allsqsaverage, type = 'o', col = g.col[q], 
          pch = 21, bg = "grey")
      }
    }
  }
  for(n in 1:length(form_bins)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(form_bins)){
        if(nrow(binlist) %% 2 == 0){
          next
        }
        else{
          rect(useful_bins[n], 0, useful_bins[n-1], 
          (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1)), 
          col = "#32323232", border = NA)
        }
      }
      else{
        rect(form_bins[n], 0, form_bins[n+1], 
        (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1)), 
        col = "#32323232", border = NA)
      }
    }
  }
  beepr::beep(sound = 3)
  proc.time() - ptm
}

    