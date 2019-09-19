#=============================================== FUNCTIONS FOR FORMATION BINNER =================================================
#
#      Selection of functions that work with PBDB and Formation data to set up Formation Binning. Functions include those
#      that visualise or retrieve data, those that setup and run Formation Binning, and those which run diversity analyses
#      on inputted data with varying bin types. 
#
#=============================================== iPAK AND REQUIRED PACKAGES =====================================================

# fucntion that automatically installs necessary packages that the user is lacking.

#===== iPAK =====
ipak <- function(pkg){ # Function to install packages. Read in character vector of any packages required. 
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} 

#=============================================== GET_INFO ============================================================

# Retrieves basic information about inputted formations, including mean age and range.
GetInfo <- function(formations){
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
# If a formation does not cross that boundary, the formation is given a score of 100.If it does cross that boundary,
# it works out how much of the formation sits each side of the boundary, and downweights accordingly - e.g. if a 
# formation was found to have 10% of it's total range one side of a boundary and 90% the other, it would receive a
# higher score than a formation where 50% of it's range sat either side. Mean scores are generated for each time bin.
      
#==== Scoring_Grid_1 ====

# Creates a scoring grid using info from all formations.
Scoring_Grid_1 <- function(formations, res=0.01) { #includes all formations
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
  score_grid <<- score_grid # Outputs score_grid for later use
  allbins <<- allbins # Outputs bins for later use
}

#==== Scoring_Grid_2 ====

# Create a scoring grid ignoring formations with length longer than mean formation length. 
Scoring_Grid_2 <- function(formations, res=0.01) { # ignores formations longer than mean range
  max_age <- max(formations$max_age) #finds max age of all formations
  min_age <- min(formations$min_age) #finds min age of all formations
  allbins <- seq(min_age-1.0045, max_age+1.0045, res) # Makes 1ma bins in sequence based on max/min ages. 0.0001 added to ensure formation is never exactly equivalent to a bin.
  
  score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(allbins)) # makes a matrix for the scoring 
  # of each time line in terms of how good it is to be a bin boundary
  colnames(score_grid) <- allbins # All time lines
  rownames(score_grid) <- formations$Formation # All formations
  
  counter <- 0
  for(i in allbins) { # Go through each time line
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
  allbins <<- allbins
}

#=============================================== PLOTMAKER =======================================================

# Generates plots through time with user inputted data and Formation_Bins, whilst providing traditional stage 
# data for comparison. NOTE: xlim is specified to fit the chosen time window of this study - as such, this would
# have to be adjusted if other data were to be used.
plotMaker <- function(rel_data, binlist, ulabel){
  useful_bins <- c(binlist$bottom, binlist$top[nrow(binlist)])
  tsplot(stages, boxes=c("short","system"), ylab = ulabel,
         xlim=75:81,  ylim=c(0,(max(rel_data, na.rm = TRUE)+(max(rel_data, na.rm = TRUE)*0.1))), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  
  for(n in 1:length(useful_bins)){
    if(((n %% 2) == 0) == TRUE) next
    else {
      if(n == length(useful_bins)){
        if(nrow(binlist) %% 2 == 0){
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
  lines(binlist$mid, rel_data, type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)
}

#=============================================== NEWBINS ==============================================================

# Looks at the previously generated score_grid and graphically highlights all places where it would be suitable to draw
# a boundary, above a user specified threshold (thresh). Boundaries are outputted as a list (binlist).

newBins <- function(score_grid, formations, bin_limits, allbins, stages){
  # Set up
  score_grid<- as.data.frame(score_grid)
  bin_size <- bin_limits[1]
  max_age <- bin_limits[2]
  min_age <- bin_limits[3]
  form_bins <- c(min_age)
  testbin <- seq(min_age, max_age, bin_size) # Creates broader bins for testing best point to draw a bin
  
  # Drawing bins and giving form_bins (vector of bin boundaries)
  for (i in 1:length((testbin)-1)){
    seqs <- seq(testbin[i],testbin[i]+bin_size, 1)
    pasting <- c()
    for (n in 1:bin_size){
      pasting <- c(pasting, paste("^",seqs[n],"|", sep = ""))
    }
    testmatch2 <- paste(pasting, collapse = "")
    testmatch2 <- substr(testmatch2, 1, nchar(testmatch2)-1) 
    a <- score_grid[grep(testmatch2, names(score_grid))] 
    z <- apply(a,1,which.max) 
    form_bins[i+1] <- names(a)[z][nrow(a)]
  }
  form_bins <- as.numeric(unique(form_bins))
  range <- (diff(form_bins) < 0.5)
  range_checker <- c()
  for (r in 1:length(range)){
    if (range[r] == TRUE){
      difference <- diff(c(form_bins[r], form_bins[r+1]))
      warning("Original bin ",  r, " removed due to small range: ~", signif(difference, digits = 3), " Ma. The difference in time has been added to the bins above and below.")
      form_bins[r] <- form_bins[r]+(difference/2)
      form_bins[r+1] <- form_bins[r+1]-(difference/2)
      range_checker <- c(range_checker, r)
    }
  }
  if(length(range_checker) > 0){
    form_bins <- form_bins[-range_checker]
  }
  form_bins <<- form_bins
  
  # Creating binlist (data.frame of bins and appropriate age info)
  prefix <- "FB."
  suffix <- seq(1:(length(form_bins)-1))
  my_names <- paste(prefix, suffix, sep = "")
  binlist <- data.frame(bin = my_names, 
                        bottom = as.numeric(form_bins[1:(length(form_bins)-1)]), 
                        top = as.numeric(form_bins[2:(length(form_bins))]))
  binlist$mid <- (binlist$bottom + binlist$top) / 2
  binlist <<- binlist

  #Plotting new bins
  par(mar = c(4.1, 4.1, 1, 2.1))
  tsplot(stages, boxes=c("short","system"),
         xlim=75:81,  ylim=c(min(colMeans(score_grid), na.rm = TRUE), 100), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Bin Splitting Score") 
  lines(allbins, colMeans(score_grid))
  
  
  for(n in 1:length(form_bins)){
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
FormationGraph <- function(formations, form_bins, stages){
  fp <- data.frame("x1" = formations$min_age, "y1" = formations$forbinning, 'x2' = formations$max_age, 'y2' = formations$forbinning)
  par(mar = c(4.1, 4.1, 1, 2.1))
  tsplot(stages, boxes=c("short","system"),
         xlim=75:81,  ylim=range(fp$y1, fp$y2), 
         shading=NULL, boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
         ylab = "Formations by Number")
  segments(fp$x1, fp$y1, fp$x2, fp$y2, lwd = 2)
  
  for(n in 1:length(form_bins)){
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
# that they occur in. Produces graphs showing raw dinosaur occurrences, rock outcrop area, and DPSK through time.

FormBin_M1 <- function(formations, binlist, Form_list, Quorum) {
  sqsmst <- list()
  for (q in 1:length(Quorum)){
    M1_Dino_List <- list()# make an empty list of the dinos in each bin

    # Code for assigning formations and associated occurrences to bins
    for(b in 1:nrow(binlist)){ # for each new bin
      temp_dino_recs <- data.frame()
      for (f in 1:nrow(formations)){ 
        if (formations$max_age[f] >= binlist[b,2] && formations$min_age[f] <= binlist[b,3]){ # If Formation max. age is greater than Bin min. age AND if Formation min. age is less then Bin max. age (i.e. falls within bin at some point)
          temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
        }
      }
      if (nrow(temp_dino_recs) > 0){
        temp_dino_recs$bin_no <- b
      }
      M1_Dino_List[[b]] <- temp_dino_recs
    }
    df <- do.call("rbind", M1_Dino_List)

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
# Produces graphs showing raw dinosaur occurrences, rock outcrop area, and DPSK through time.
FormBin_M2<- function(formations, binlist, Form_list, Quorum) {
  sqsmst <- list()
  for (q in 1:length(Quorum)){
    M2_Dino_List <- list()# make an empty list of the dinos in each bin
    for(b in 1:nrow(binlist)){ # for each new bin
      temp_dino_recs <- data.frame()
      for (f in 1:nrow(formations)) { # for each formation
        if (formations$max_age[f] < binlist[b,2] |
            binlist[b,3] < formations$min_age[f]) { # If the formation DOES NOT sit in this bin 
          next # Just print, don't do anything else.
        }
        else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
          if (formations$max_age[f] <= binlist[b,3] && 
              formations$min_age[f] >= binlist[b,2]){# If formation sits within boundaries
            temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
            next
          }
          if (formations$max_age[f] >= binlist[b,2] && 
              formations$min_age[f] <= binlist[b,2] && 
              formations$max_age[f] <= binlist[b,3]){ #If formation crosses upper age limit (minimum age)
            x <- as.numeric(formations$max_age[f] - as.numeric(binlist[b,2]))
            y <- as.numeric(as.numeric(binlist[b,2]) - formations$min_age[f])
            if (x > y){
              temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
              next
            }
          }
          if (formations$max_age[f] >= binlist[b,3] && 
              formations$min_age[f] <= binlist[b,3] &&
              formations$min_age[f] >= binlist[b,2]){ #If formation crosses lower age limit (maximum age)
            x <- as.numeric(formations$max_age[f] - as.numeric(binlist[b,3]))
            y <- as.numeric(as.numeric(binlist[b,3]) - formations$min_age[f])
            if (y > x){
              temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
              next
            }
          }
          if (formations$max_age[f] > binlist[b,3] &&
              formations$min_age[f] < binlist[b,2]){
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
# Produces graphs showing raw dinosaur occurrences, rock outcrop area, and DPSK through time. 
FormBin_M3<- function(formations, binlist, Form_list, times=10, Quorum) {
  ptm <- proc.time()
  for (q in 1:length(Quorum)){
    allSQS <- data.frame(binlist$bin)
    allbininfo <- data.frame(binlist$bin)
    for(n in 1:times){
      M3_Dino_List <- list()
      for(b in 1:nrow(binlist)){ # for each new bin
        temp_dino_recs <- data.frame() 
        for (f in 1:nrow(formations)) { # for each formation
          if (formations$max_age[f] < binlist[b,2] |
              binlist[b,3] < formations$min_age[f]) { # If the formation DOES NOT sit in this bin 
            next # Just print, don't do anything else.
          }
          else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
            if (formations$max_age[f] <= binlist[b,3] && 
              formations$min_age[f] >= binlist[b,2]){ # If formation sits within boundaries
              temp_dino_recs <- rbind(temp_dino_recs, Form_list[[f]]) # Add dinos from that formation to dino list
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
              temp_dino_recs <- rbind(temp_dino_recs, sampled) # Add dinos from that formation to dino list
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
              temp_dino_recs <- rbind(temp_dino_recs, sampled) # Add dinos from that formation to dino list
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
              temp_dino_recs <- rbind(temp_dino_recs, sampled) # Add dinos from that formation to dino list
              next
            } 
          }
        }
        if (nrow(temp_dino_recs) > 0){ # If there are occurrences in this bin
          temp_dino_recs$bin_no <- b # Label those occurrences with bin number
        }
        M3_Dino_List[[b]] <- temp_dino_recs # Add temp records to permanent list
      }
      df <- do.call("rbind", M3_Dino_List) # Create data.frame from list
        
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

    