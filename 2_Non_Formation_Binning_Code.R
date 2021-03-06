#=============================================================================================================================
#============================================ 2. TRADITIONAL BINNING METHODS =================================================
#=============================================================================================================================

# Christopher D. Dean, A. Alessandro Chiarenza & Susannah Maidment 
# 2019

#================================================ DATA SETUP =================================================================

#===== 0. Packages to load =====
library(matrixStats)
library(divDyn)
library(colorspace)

#=== Data entry and setup ===
occs <- read.csv(file = "Data/Occurrences_Final.csv") # Read in occurrences
formations <- read.csv (file = "Data/Formations_Final.csv")  #Read in formations
data(stages)

# Make parts Numeric
formations$max_age <- as.numeric(as.character(formations$max_age)) # Make Numeric
formations$min_age <- as.numeric(as.character(formations$min_age)) # Make Numeric 

# Select appropriate formations and order
formations <- formations[which(formations$Location=='WI'),] # Only formations from Western Interior
myformations <- sort(as.vector(formations$Formation)) # Organise

# Select appropriate occurrences
occs <- occs[occs$formation %in% myformations,] # Only include occurrences from formation list

# Create Formation/occurrences list
Form_list <- split(occs, occs$formation)

colnames(occs)[38] <- "new_bin"
occs$new_bin <- as.numeric(as.character(occs$new_bin))

plotrect <- function(binDframe, ydata){
  for(n in 1:nrow(binDframe)){ # For each bin
    if(((n %% 2) == 0) == TRUE) next # If the bin is even, skip it (allows for alternating colours of bins)
    else {
      if(n == nrow(binDframe)){
        if(nrow(binDframe) %% 2 == 0){ # Skips colouring last bin if there are an even number of bins
          next
        }
        else{
          rect(binDframe$LAD[n], 0, binDframe$FAD[n], # Dimensions of the coloured bin
               (max(ydata, na.rm = TRUE)+(max(ydata, na.rm = TRUE)*0.1)), 
               col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
        }
      }
      else{
        rect(binDframe$LAD[n], 0, binDframe$FAD[n], # Dimensions of the coloured bin
             (max(ydata, na.rm = TRUE)+(max(ydata, na.rm = TRUE)*0.1)),
             col = rgb(0.89,0.89,0.89,alpha=0.5), border = NA)
      }
    }
  }
}

#================================================ BIN SETUP ==================================================================

#=== Binning at Stage level ===

bins <- c(113.0, 100.5, 93.9, 89.8, 86.3, 83.6, 72.1, 66)

binDframe <- data.frame(bin = c("Albian", 
                                "Cenomanian",
                                "Turonian",
                                "Coniacian",
                                "Santonian",
                                "Campanian",
                                "Maastrichtian"), # Combines bin data to make dataframe of minimum, maximum and mid point of each new bin
                        FAD = as.numeric(bins[1:(length(bins)-1)]), 
                        LAD = as.numeric(bins[2:(length(bins))]))
binDframe$mid <- (binDframe$FAD+binDframe$LAD)/2

#=== Binning at sub-Stage level ===

bins <- c(113.0, 110.2, 107.6, 100.5, 98.5, 96.3, 93.9, 92.9, 91.5, 89.8, 88.5, 88.1, 86.3, 86, 83.6, 80.6, 76.3, 72.1, 69.9, 66)

binDframe <- data.frame(bin = c("Lower Albian", "Middle Albian", "Upper Albian", 
                                "Lower Cenomanian", "Middle Cenomanian", "Upper Cenomanian",
                                "Lower Turonian", "Middle Turonian", "Upper Turonian",
                                "Lower Coniacian", "Middle Coniacian", "Upper Coniacian",
                                "Lower Santonian", "Upper Santonian", 
                                "Lower Campanian", "Middle Campanian", "Upper Campanian",
                                "Lower Maastrichtian", "Upper Maastrichtian"), # Combines bin data to make dataframe of minimum, maximum and mid point of each new bin
                        FAD = as.numeric(bins[1:(length(bins)-1)]), 
                        LAD = as.numeric(bins[2:(length(bins))]))
binDframe$mid <- (binDframe$FAD+binDframe$LAD)/2

#=============================================================================================================================
#================================================ PBDB AGES ==================================================================
#=============================================================================================================================

#================================================ BINNING METHOD =============================================================

#=== Bin by midpoint ===
for(s in 1:nrow(binDframe)){
  for(o in 1:nrow(occs)){
    if(occs$ma_mid[o] <= binDframe$FAD[s] && occs$ma_mid[o] >= binDframe$LAD[s]){
      occs$new_bin[o] <- stages$stg[s]
    }
  }
}
df <- occs

#=== Bin in any bin it falls within ===
All_Bin_List <- list()
for(s in 1:nrow(binDframe)){ # for each new bin
  temp_recs <- data.frame()
  for (o in 1:nrow(occs)){ 
    if (occs$ma_max[o] >= binDframe$LAD[s] && occs$ma_min[o] <= binDframe$FAD[s]){ # If occurrence max. age is greater than Bin min. age AND if occurrence min. age is less then Bin max. age (i.e. falls within bin at some point)
      temp_recs <- rbind(temp_recs, occs[o,]) # Add that occurrence to binlist
    }
  }
  if (nrow(temp_recs) > 0){
    temp_recs$new_bin <- s
  }
  All_Bin_List[[s]] <- temp_recs
}
df <- do.call("rbind", All_Bin_List)

#================================================ RAW DIVERSITY AND SQS ======================================================

bin_info <- binstat(df, tax="occurrence.genus_name", bin="new_bin", 
                    coll = 'collection_no')
bin_info[bin_info==0] <- NA

sqsmst <- list()
Quorum <- c(0.4, 0.6, 0.8)
for (q in 1:length(Quorum)){
  SQS <- subsample(df,iter=100, q=Quorum[q], tax="occurrence.genus_name", bin="new_bin", 
                   coll = 'collection_no', output="dist", type="sqs", 
                   duplicates = TRUE, useFailed = TRUE)
  SIBSQS <- SQS$divSIB
  MeanSQS <- rowMeans2(SIBSQS)
  SDSQS <- rowSds(SIBSQS)
  MaxSDSQS <- MeanSQS + SDSQS
  MinSDSQS <- MeanSQS - SDSQS
  
  combined.sqs <- cbind(MeanSQS, SDSQS, MaxSDSQS, MinSDSQS)
  rownames(combined.sqs) <- binDframe$bin
  combined.sqs <- as.data.frame(combined.sqs)
  
  temp_name <- paste("q.",deparse(Quorum[q]),"_", "SQS_Results", sep = "") #Name files based on data entered to function
  assign(temp_name, combined.sqs, envir = .GlobalEnv)
  sqsmst[[q]] <- combined.sqs
  names(sqsmst)[[q]] <- deparse(Quorum[q])
}

#================================================ PLOTTING AND RECORDING RESULTS =============================================

# Diversity and Collections
par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Raw Diversity",
       xlim=75:81,  ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
        boxes.col=c("col","systemCol"), labels.args=list(cex=1.8))
plotrect(binDframe, bin_info$SIBs)
lines(binDframe$mid, bin_info$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1.5)
box(lwd=2)

tsplot(stages, boxes=c("short","system"), ylab = "Number of Collections",
       xlim=75:81,  ylim=c(0,(max(bin_info$colls, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
      boxes.col=c("col","systemCol"), labels.args=list(cex=1.8))
plotrect(binDframe, bin_info$colls)
lines(binDframe$mid, bin_info$colls, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1.5)
box(lwd=2)

# Goods u
par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Good's U", # Creates plot using data from DivDyn package. 
       xlim=75:81,  ylim=c(0,(max(bin_info$u, na.rm = TRUE)+(max(bin_info$u, na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
       boxes.col=c("col","systemCol"), labels.args=list(cex=1.8))  
plotrect(binDframe, bin_info$u)
lines(binDframe$mid, bin_info$u, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1.5)
box(lwd=2)

#=== Plotting SQS ===
tsplot(stages, boxes=c("short","system"), 
       xlim=75:81,  ylim=c(0,(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+
                                (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
       boxes.col=c("col","systemCol"), labels.args=list(cex=1.8),
       ylab = "Subsampled Diversity")  
plotrect(binDframe, sqsmst[[length(sqsmst)]]$MeanSQS)
for (q in 1:length(Quorum)){
  g.col <- (sequential_hcl(5, "YlGn"))
  enc <- rle(!is.na(sqsmst[[q]]$MeanSQS))
  endIdxs <- cumsum(enc$lengths)
  for(i in 1:length(enc$lengths)){
    if(enc$values[i]){
      endIdx <- endIdxs[i]
      startIdx <- endIdx - enc$lengths[i] + 1
      
      subdat <- binDframe$mid[startIdx:endIdx]
      submin <- sqsmst[[q]]$MinSDSQS[startIdx:endIdx]
      submax <- sqsmst[[q]]$MaxSDSQS[startIdx:endIdx]
      subdepth <- sqsmst[[q]]$MeanSQS[startIdx:endIdx]
      
      x <- c(subdat, rev(subdat))
      y <- c(submax, rev(submin))
      
      polygon(x = x, y = y, col = adjustcolor(g.col[q], alpha.f = 0.40), border = NA)
      lines(binDframe$mid, sqsmst[[q]]$MeanSQS, type = 'o', col = g.col[q], 
            pch = 21, bg = g.col[q], lwd = 2, cex = 1.5)
    }
  }
}
box(lwd=2)

#=== Recording Results ===
dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
write.csv(bin_info, file.path(paste("Results/StandardBinned_Bin_info.csv", sep="")))
for (q in 1:length(Quorum)){
  temp_name <- paste("StandardBinned_SQS_", Quorum[q], sep = "")
  write.csv(sqsmst[q], file.path(paste("Results/", temp_name, ".csv", sep="")))
}

#=============================================================================================================================
#================================================ FORMATION AGES =============================================================
#=============================================================================================================================

#================================================ BINNING METHOD =============================================================

#=== Setup for midpoint ===
occs$new_mid_age <- NA # Make new column for new mid ages
for(o in 1:nrow(occs)){
  for(f in 1:nrow(formations)){
    if(occs$formation[o] == formations$Formation[f]){ # If occurrence formation matches current formation in loop
      occs$new_mid_age[o] <- mean(c(formations$max_age[f], formations$min_age[f])) # Add formation mid point to new mid point column
    }
  }
}

#=== Binning by midpoint === 
for(s in 1:nrow(binDframe)){
  for(o in 1:nrow(occs)){
    if(occs$new_mid_age[o] <= binDframe$FAD[s] && occs$new_mid_age[o] >= binDframe$LAD[s]){
      occs$new_bin[o] <- s
    }
  }
}

for(o in 1:nrow(occs)){
  if(occs$new_mid_age[o] < 66){
    occs$new_bin[o] <- nrow(binDframe)
  }
}

df <- occs
df$new_bin <- as.numeric(df$new_bin)

#=== Setup for any bin ===
occs <- droplevels.data.frame(occs)
formations <- droplevels.data.frame(formations)

for(o in 1:nrow(occs)){
  for(f in 1:nrow(formations)){
    if(occs$formation[o] == formations$Formation[f]){
      occs$form_max_ma[o] <- formations$max_age[f]
      occs$form_min_ma[o] <- formations$min_age[f]
    }
  }
}

#=== Bin in any bin it falls within ===
All_Bin_List <- list()
for(s in 1:nrow(binDframe)){ # for each new bin
  temp_recs <- data.frame()
  for (o in 1:nrow(occs)){ 
    if (occs$form_max_ma[o] >= binDframe$LAD[s] && occs$form_min_ma[o] <= binDframe$FAD[s]){ # If occurrence max. age is greater than Bin min. age AND if occurrence min. age is less then Bin max. age (i.e. falls within bin at some point)
      temp_recs <- rbind(temp_recs, occs[o,]) # Add that occurrence to binlist
    }
  }
  if (nrow(temp_recs) > 0){
    temp_recs$new_bin <- s
  }
  All_Bin_List[[s]] <- temp_recs
}
df <- do.call("rbind", All_Bin_List)

#================================================ RAW DIVERSITY AND SQS ======================================================

bin_info <- binstat(df, tax="occurrence.genus_name", bin="new_bin", 
                    coll = 'collection_no')
bin_info[bin_info==0] <- NA

sqsmst <- list()
Quorum <- c(0.4, 0.6, 0.8)
for (q in 1:length(Quorum)){
  SQS <- subsample(df,iter=100, q=Quorum[q], tax="occurrence.genus_name", bin="new_bin", 
                   coll = 'collection_no', output="dist", type="sqs", 
                   duplicates = TRUE, useFailed = TRUE)
  SIBSQS <- SQS$divSIB
  MeanSQS <- rowMeans2(SIBSQS)
  SDSQS <- rowSds(SIBSQS)
  MaxSDSQS <- MeanSQS + SDSQS
  MinSDSQS <- MeanSQS - SDSQS
  
  combined.sqs <- cbind(MeanSQS, SDSQS, MaxSDSQS, MinSDSQS)
  rownames(combined.sqs) <- binDframe$bin
  combined.sqs <- as.data.frame(combined.sqs)
  
  temp_name <- paste("q.",deparse(Quorum[q]),"_", "SQS_Results", sep = "") #Name files based on data entered to function
  assign(temp_name, combined.sqs, envir = .GlobalEnv)
  sqsmst[[q]] <- combined.sqs
  names(sqsmst)[[q]] <- deparse(Quorum[q])
}
sqsmst[[1]][sqsmst[[1]] == 0] <- NA
sqsmst[[2]][sqsmst[[2]] == 0] <- NA
sqsmst[[3]][sqsmst[[3]] == 0] <- NA

#================================================ PLOTTING AND RECORDING RESULTS =============================================

# Diversity and Collections
par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Raw Diversity",
       xlim=75:81,  ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
       boxes.col=c("col","systemCol"), labels.args=list(cex=1.8))
plotrect(binDframe, bin_info$SIBs)
lines(binDframe$mid, bin_info$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1.5)
box(lwd=2)


tsplot(stages, boxes=c("short","system"), ylab = "Number of Collections",
       xlim=75:81,  ylim=c(0,(max(bin_info$colls, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
       boxes.col=c("col","systemCol"), labels.args=list(cex=1.8))
plotrect(binDframe, bin_info$colls)
lines(binDframe$mid, bin_info$colls, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1.5)
box(lwd=2)

par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Good's U", # Creates plot using data from DivDyn package. 
       xlim=75:81,  ylim=c(0,(max(bin_info$u, na.rm = TRUE)+(max(bin_info$u, na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
       boxes.col=c("col","systemCol"), labels.args=list(cex=1.8)) 
plotrect(binDframe, bin_info$u)
lines(binDframe$mid, bin_info$u, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1.5)
box(lwd=2)

#=== Plotting SQS ===
tsplot(stages, boxes=c("short","system"), 
       xlim=75:81,  ylim=c(0,(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+
                                (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1))), 
       prop = 0.08, plot.args = list(cex.lab = 2, cex.axis = 2),
       boxes.col=c("col","systemCol"), labels.args=list(cex=1.8),
       ylab = "Subsampled Diversity") 
plotrect(binDframe, sqsmst[[length(sqsmst)]]$MeanSQS)
for (q in 1:length(Quorum)){
  g.col <- (sequential_hcl(5, "YlGn"))
  enc <- rle(!is.na(sqsmst[[q]]$MeanSQS))
  endIdxs <- cumsum(enc$lengths)
  for(i in 1:length(enc$lengths)){
    if(enc$values[i]){
      endIdx <- endIdxs[i]
      startIdx <- endIdx - enc$lengths[i] + 1
      
      subdat <- binDframe$mid[startIdx:endIdx]
      submin <- sqsmst[[q]]$MinSDSQS[startIdx:endIdx]
      submax <- sqsmst[[q]]$MaxSDSQS[startIdx:endIdx]
      subdepth <- sqsmst[[q]]$MeanSQS[startIdx:endIdx]
      
      x <- c(subdat, rev(subdat))
      y <- c(submax, rev(submin))
      
      polygon(x = x, y = y, col = adjustcolor(g.col[q], alpha.f = 0.40), border = NA)
      lines(binDframe$mid, sqsmst[[q]]$MeanSQS, type = 'o', col = g.col[q], 
            pch = 21, bg = g.col[q], lwd = 2, cex = 1.5)
    }
  }
}
box(lwd=2)

#=== Recording Results ===
dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
write.csv(bin_info, file.path(paste("Results/SB_NewFormAge_Bin_info.csv", sep="")))
for (q in 1:length(Quorum)){
  temp_name <- paste("SB_NewFormAge_SQS_", Quorum[q], sep = "")
  write.csv(sqsmst[q], file.path(paste("Results/", temp_name, ".csv", sep="")))
}

