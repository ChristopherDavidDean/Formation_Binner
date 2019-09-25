#===== Running Traditional Binning Methods in SQS =====
# christopher D. Dean, Alfio A. Chiarenza, Susannah Maidment 

# Packages to load
library(matrixStats)
library(divDyn)

#===== 1. Diversity based on PBDB ages =====

#=== Set Up ===
occs <- read.csv(file = "Data/NADINOS-occs-edit.csv") # Read in occurrences
formations <- read.csv (file = "Data/Formations_test2.csv")  #Read in formations
data(stages)

# Make parts Numeric
formations$max_age <- as.numeric(as.character(formations$max_age)) # Make Numeric
formations$min_age <- as.numeric(as.character(formations$min_age)) # Make Numeric 

# Select appropriate formations and order
formations <- formations[which(formations$Location=='WI'),] # Only formations from Western Interior
myformations <- sort(as.vector(formations$Formation)) # Organise

# Select appropriate occurrences
occs <- occs[occs$formation %in% myformations,] # Only include occurrences from formation list
occs <- droplevels.data.frame(testoccs) # Remove old levels

colnames(occs)[38] <- "new_bin"
occs$new_bin <- as.numeric(as.character(occs$new_bin))

#=== Binning ===
for(s in 1:nrow(stages)){
  for(o in 1:nrow(occs)){
    if(occs$ma_mid[o] < stages$bottom[s] && occs$ma_mid[o] > stages$top[s]){
      occs$new_bin[o] <- stages$stg[s]
    }
  }
}

#=== Raw Diversity and SQS ===
bin_info <- binstat(occs, tax="occurrence.genus_name", bin="new_bin", 
                    coll = 'collection_no')
bin_info[bin_info==0] <- NA

sqsmst <- list()
Quorum <- c(0.4, 0.6, 0.8)
for (q in 1:length(Quorum)){
  SQS <- subsample(occs,iter=100, q=Quorum[q], tax="occurrence.genus_name", bin="new_bin", 
                   coll = 'collection_no', output="dist", type="sqs", 
                   duplicates = TRUE, useFailed = TRUE)
  SIBSQS <- SQS$divSIB[75:81,]
  MeanSQS <- rowMeans2(SIBSQS)
  SDSQS <- rowSds(SIBSQS)
  MaxSDSQS <- MeanSQS + SDSQS
  MinSDSQS <- MeanSQS - SDSQS
  
  combined.sqs <- cbind(MeanSQS, SDSQS, MaxSDSQS, MinSDSQS)
  rownames(combined.sqs) <- stages$stage[75:81]
  combined.sqs <- as.data.frame(combined.sqs)
  
  temp_name <- paste("q.",deparse(Quorum[q]),"_", "SQS_Results", sep = "") #Name files based on data entered to function
  assign(temp_name, combined.sqs, envir = .GlobalEnv)
  sqsmst[[q]] <- combined.sqs
  names(sqsmst)[[q]] <- deparse(Quorum[q])
}

#=== Plotting raw diversity, collections and Good's U ===
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
par(mar = c(0, 4.1, 4.1, 2.1))
with(stages[75:81,], tsplot(stages, ylab = "Raw Diversity",
                     xlim=75:81,  ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
                     shading="stage", plot.args=list(xaxt = 'n')))

lines(stages$mid[75:81], bin_info$SIBs[75:81], type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)
par(mar = c(4.1, 4.1, 0, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Number of Collections", # Creates plot using data from DivDyn package. 
       xlim=75:81,  ylim=c(0,(max(bin_info$colls[75:81], na.rm = TRUE)+(max(bin_info$colls[75:81], na.rm = TRUE)*0.1))), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  
lines(stages$mid[75:81], bin_info$colls[75:81], type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)

par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Good's U", # Creates plot using data from DivDyn package. 
       xlim=75:81,  ylim=c(0,(max(bin_info$u[75:81], na.rm = TRUE)+(max(bin_info$u[75:81], na.rm = TRUE)*0.1))), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  
lines(stages$mid[75:81], bin_info$u[75:81], type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)

#=== Plotting SQS ===
tsplot(stages, boxes=c("short","system"), 
       xlim=75:81,  ylim=c(0,(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+
                                (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1))), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
       ylab = "Subsampled Diversity")  
for (q in 1:length(Quorum)){
  g.col <- gray.colors(length(Quorum), start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL)
  enc <- rle(!is.na(sqsmst[[q]]$MeanSQS))
  endIdxs <- cumsum(enc$lengths)
  for(i in 1:length(enc$lengths)){
    if(enc$values[i]){
      endIdx <- endIdxs[i]
      startIdx <- endIdx - enc$lengths[i] + 1
      
      subdat <- stages$mid[75:81][startIdx:endIdx]
      submin <- sqsmst[[q]]$MinSDSQS[startIdx:endIdx]
      submax <- sqsmst[[q]]$MaxSDSQS[startIdx:endIdx]
      subdepth <- sqsmst[[q]]$MeanSQS[startIdx:endIdx]
      
      x <- c(subdat, rev(subdat))
      y <- c(submax, rev(submin))
      
      polygon(x = x, y = y, col = adjustcolor(g.col[q], alpha.f = 0.40), border = NA)
      lines(stages$mid[75:81], sqsmst[[q]]$MeanSQS, type = 'o', col = g.col[q], 
            pch = 21, bg = "grey")
    }
  }
}

#=== Recording Results ===
bin_info <- bin_info[75:81,]
dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
write.csv(bin_info, file.path(paste("Results/StandardBinned_Bin_info.csv", sep="")))
for (q in 1:length(Quorum)){
  temp_name <- paste("StandardBinned_SQS_", Quorum[q], sep = "")
  write.csv(sqsmst[q], file.path(paste("Results/", temp_name, ".csv", sep="")))
}

#===== 1. Diversity based on updated formation ages =====

#=== Setup ===
occs$new_mid_age <- NA # Make new column for new mid ages
for(o in 1:nrow(occs)){
  for(f in 1:nrow(formations)){
    if(occs$formation[o] == formations$Formation[f]){ # If occurrence formation matches current formation in loop
      occs$new_mid_age[o] <- mean(c(formations$max_age[f], formations$min_age[f])) # Add formation mid point to new mid point column
    }
  }
}

#=== Binning === 
for(s in 1:nrow(stages)){
  for(o in 1:nrow(occs)){
    if(occs$new_mid_age[o] < stages$bottom[s] && occs$new_mid_age[o] > stages$top[s]){
      occs$new_bin[o] <- stages$stg[s]
    }
  }
}

#=== Raw Diversity and SQS ===
bin_info <- binstat(occs, tax="occurrence.genus_name", bin="new_bin", 
                    coll = 'collection_no')
bin_info[bin_info==0] <- NA

sqsmst <- list()
Quorum <- c(0.4, 0.6, 0.8)
for (q in 1:length(Quorum)){
  SQS <- subsample(occs,iter=100, q=Quorum[q], tax="occurrence.genus_name", bin="new_bin", 
                   coll = 'collection_no', output="dist", type="sqs", 
                   duplicates = TRUE, useFailed = TRUE)
  SIBSQS <- SQS$divSIB[75:81,]
  MeanSQS <- rowMeans2(SIBSQS)
  SDSQS <- rowSds(SIBSQS)
  MaxSDSQS <- MeanSQS + SDSQS
  MinSDSQS <- MeanSQS - SDSQS
  
  combined.sqs <- cbind(MeanSQS, SDSQS, MaxSDSQS, MinSDSQS)
  rownames(combined.sqs) <- stages$stage[75:81]
  combined.sqs <- as.data.frame(combined.sqs)
  
  temp_name <- paste("q.",deparse(Quorum[q]),"_", "SQS_Results", sep = "") #Name files based on data entered to function
  assign(temp_name, combined.sqs, envir = .GlobalEnv)
  sqsmst[[q]] <- combined.sqs
  names(sqsmst)[[q]] <- deparse(Quorum[q])
}
sqsmst[[1]][sqsmst[[1]] == 0] <- NA
sqsmst[[2]][sqsmst[[2]] == 0] <- NA
sqsmst[[3]][sqsmst[[3]] == 0] <- NA

#=== Plotting raw diversity, collections and Good's U ===
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
par(mar = c(0, 4.1, 4.1, 2.1))
with(stages[75:81,], tsplot(stages, ylab = "Raw Diversity",
                            xlim=75:81,  ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))), 
                            shading="stage", plot.args=list(xaxt = 'n')))

lines(stages$mid[75:81], bin_info$SIBs[75:81], type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)
par(mar = c(4.1, 4.1, 0, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Number of Collections", # Creates plot using data from DivDyn package. 
       xlim=75:81,  ylim=c(0,(max(bin_info$colls[75:81], na.rm = TRUE)+(max(bin_info$colls[75:81], na.rm = TRUE)*0.1))), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  
lines(stages$mid[75:81], bin_info$colls[75:81], type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)

par(mfrow=c(1,1), mar = c(4.1, 4.1, 1, 2.1))
tsplot(stages, boxes=c("short","system"), ylab = "Good's U", # Creates plot using data from DivDyn package. 
       xlim=75:81,  ylim=c(0,(max(bin_info$u[75:81], na.rm = TRUE)+(max(bin_info$u[75:81], na.rm = TRUE)*0.1))), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75))  
lines(stages$mid[75:81], bin_info$u[75:81], type = "o", pch = 21, col = "black", bg = "grey", lwd = 1)

#=== Plotting SQS ===
tsplot(stages, boxes=c("short","system"), 
       xlim=75:81,  ylim=c(0,(max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)+
                                (max(sqsmst[length(sqsmst)][[1]][1], na.rm = TRUE)*0.1))), 
       shading="stage", boxes.col=c("col","systemCol"), labels.args=list(cex=0.75),
       ylab = "Subsampled Diversity")  
for (q in 1:length(Quorum)){
  g.col <- gray.colors(length(Quorum), start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL)
  enc <- rle(!is.na(sqsmst[[q]]$MeanSQS))
  endIdxs <- cumsum(enc$lengths)
  for(i in 1:length(enc$lengths)){
    if(enc$values[i]){
      endIdx <- endIdxs[i]
      startIdx <- endIdx - enc$lengths[i] + 1
      
      subdat <- stages$mid[75:81][startIdx:endIdx]
      submin <- sqsmst[[q]]$MinSDSQS[startIdx:endIdx]
      submax <- sqsmst[[q]]$MaxSDSQS[startIdx:endIdx]
      subdepth <- sqsmst[[q]]$MeanSQS[startIdx:endIdx]
      
      x <- c(subdat, rev(subdat))
      y <- c(submax, rev(submin))
      
      polygon(x = x, y = y, col = adjustcolor(g.col[q], alpha.f = 0.40), border = NA)
      lines(stages$mid[75:81], sqsmst[[q]]$MeanSQS, type = 'o', col = g.col[q], 
            pch = 21, bg = "grey")
    }
  }
}

#=== Recording Results ===
bin_info <- bin_info[75:81,]
dir.create(paste0("Results"), showWarnings = FALSE) #stops warnings if folder already exists
write.csv(bin_info, file.path(paste("Results/SB_NewFormAge_Bin_info.csv", sep="")))
for (q in 1:length(Quorum)){
  temp_name <- paste("SB_NewFormAge_SQS_", Quorum[q], sep = "")
  write.csv(sqsmst[q], file.path(paste("Results/", temp_name, ".csv", sep="")))
}
