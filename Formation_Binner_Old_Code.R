# FORMATION BINNER TEST CODE
# Christopher Dean & Susannah Maidment
# 03/06/16

####################         ####################
####################   OLD   ####################
####################         ####################

max_age <- max(formations$max_age) #finds max age of all formations
min_age <- min(formations$min_age) #finds min age of all formations
bins <- seq(min_age, max_age+1, 1) # Makes 1ma bins in sequence based on max/min ages.

score_grid <- matrix(data = NA, nrow = nrow(formations), ncol = length(bins)) # makes a matrix for the scoring 
                               # of each time line in terms of how good it is to be a bin boundary
colnames(score_grid) <- bins # All time lines
rownames(score_grid) <- formations$Formation # All formations


# WITHOUT DPSK weighting
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

## Working out/plotting optimal bin positions
threshold <- 95
threshold_list <- c()

for(m in 2:(ncol(score_grid)-1)) {
  if (score_grid[26,m] > threshold) {
    if (score_grid[26,(m)] >= score_grid[26,(m-1)] && score_grid[26,(m)] >= score_grid[26,(m+1)]) {
       threshold_list <- c(threshold_list,(colnames(score_grid)[m]))
      }
  }
}

threshold_list # Reports all instances above a desired threshold that a time bin could be placed at.

plot(bins, means, type = 'l', xlim=rev(range.default(bins))) # plot, obviously. 
abline(v=threshold_list, lty=2, col="red")

#### EXTRA/TESTS ####

# Test for assigning formations to all relevant bins
Bin_Dino_List <- list(0,0,0,0,0,0,0) # make an empty list of the dinos in each bin
Formation_Counter_List <- list() 
Bin_Form_List <- list(0,0,0,0,0,0,0)

for(b in 1:length(binlist)){ # for each new bin
  temp_form_list <- c()
  for (f in 1:nrow(formations)) { # for each formation
    if (formations$max_age[f] > binlist[[b]][1] && formations$min_age[f] > binlist[[b]][1] | 
        binlist[[b]][2] > formations$max_age[f]  && binlist[[b]][2] > formations$min_age[f]) { # If the formation DOES NOT sit in this bin (you can work out max/min things, isn't too hard)
      print(f) # Just print, don't do anything else.
    }
    else { # Otherwise (i.e. if a formation DOES sit in/cross this bin in any way)
      temp_form_list <- c(temp_form_list, as.character(formations$Formation[f]))
      temp_dino_list <- c(temp_dino_list, as.character(Form_list[[f]][[1]])) # Add dinos from that formation to dino list
      form_counter <- form_counter + 1 # Add dinos from that formation to dino list
    }
  }
  Bin_Form_List[[b]] <- unique(temp_form_list) # Adds all formations found in one bin to that bin's list.
  Bin_Dino_List[[b]] <- unique(temp_dino_list) # This bit kinda broken I think. Supposed to add all the unique dinos from a time bin to a final list.
  Formation_Counter_List[[b]] <- form_counter
}

max_outcrop_lis <- lapply(Max_Outcrop_list3,mean)
average_outcrop_max <- mapply(c, average_outcrop_max, max_outcrop_lis, SIMPLIFY = FALSE)

###### EXTRA - this just shows how you could edit to downweight certain formations. Skip to next section! ########

# WITH DPSK weighting
counter <- 0
for(i in bins) { # Go through each time line
  counter <- sum(counter,1)
  for (f in 1:nrow(formations)){ # go through each formation 
    if (i <= formations$max_age[f] && i >= formations$min_age[f]){ # if timeline is between max/min age of a formation (i.e. formation crosses that line)
      # Need to translate i to position in grid - think this should work
      a <- formations$max_age[f] - i
      b <- i - formations$min_age[f]
      range <- formations$max_age[f] - formations$min_age[f] # Work out percentage that sits each side of line, reduce score by that amount.
      if (a > b){
        score_grid[,counter][f] <- ((a/range)*100)*formations[,9][f]
      }
      else{ 
        score_grid[,counter][f] <- ((b/range)*100)*formations[,9][f]
      }
    }
    else {
      score_grid[,counter][f] = 100*formations[,9][f] # else gets a normal score.
    }
  }
}



  average_outcrop_max <- lapply(Max_Outcrop_list3,mean)
  average_outcrop_min <- lapply(Min_Outcrop_list3,mean)



Frontier <- read.csv("0.Data/Frontier.csv")
Milk_River <- read.csv("0.Data/Milk_River.csv")
Two_Medicine <- read.csv("0.Data/Two_Medicine.csv")
Mesaverde <- read.csv("0.Data/Mesaverde.csv")
Claggett <- read.csv("0.Data/Claggett.csv")
Foremost <- read.csv("0.Data/Foremost.csv")
Judith_River <- read.csv("0.Data/Judith_River.csv")
Pierre <- read.csv("0.Data/Pierre.csv")
Wapiti <- read.csv("0.Data/Wapiti.csv")
Oldman <- read.csv("0.Data/Oldman.csv")
Dino_Park <- read.csv("0.Data/Dinosaur_Park.csv")
Williams_Fork <- read.csv("0.Data/Williams_Fork.csv")
Bearpaw <- read.csv("0.Data/Bearpaw.csv")
Horseshoe_Canyon <- read.csv("0.Data/Horseshoe_Canyon.csv")
Almond <- read.csv("0.Data/Almond.csv")
St_Mary_River <- read.csv("0.Data/St_Mary_River.csv")
Laramie <- read.csv("0.Data/Laramie.csv")
Pinyon_Conglomerate <- read.csv("0.Data/Pinyon_Conglomerate.csv")
Hell_Creek <- read.csv("0.Data/Hell_Creek.csv")
Lance <- read.csv("0.Data/Lance.csv")
Scollard <- read.csv("0.Data/Scollard.csv")
Denver <- read.csv("0.Data/Denver.csv")
Willow_Creek <- read.csv("0.Data/Willow_Creek.csv")
Frenchman <- read.csv("0.Data/Frenchman.csv")
Ferris <- read.csv("0.Data/Ferris.csv")

Form_list <- list(Frontier, Milk_River, Two_Medicine,Mesaverde,Claggett,Foremost,Judith_River,Pierre,Wapiti,Oldman,
                  Dino_Park,Williams_Fork,Bearpaw,Horseshoe_Canyon,Almond,St_Mary_River,Laramie,Pinyon_Conglomerate,
                  Hell_Creek,Lance,Scollard,Denver,Willow_Creek,Frenchman,Ferris)

# Counting in Bins OLD TEST
#a <- c(97.043,86.043)
#b <- c(86.043,83.043)
#c <- c(83.043,75.043)
#d <- c(75.043,73.043)
#e <- c(73.043,71.043)
#f <- c(71.043,68.043)
#g <- c(68.043,66.043)

#binlist <- list(a,b,c,d,e,f,g) # make a list of the bins

# Make list using values from threshold
threshold_list <- rev(threshold_list)         
binlist <- list()
for (i in 1:length(threshold_list)) {
    binlist[[i]] <- c(as.numeric(threshold_list[i]), as.numeric(threshold_list[i+1]))
 }
binlist[[length(binlist)]][2] <- as.numeric(min(formations$min_age))


### GENERATING PLOTS ###

Species_M1 <- c() # METHOD 1
for (f in 1:length(Bin_Dino_List1)) {
    Species_M1 <- c(Species_M1, length(Bin_Dino_List1[[f]]))
}
Species_M2 <- c() # METHOD 2
for (k in 1:length(Bin_Dino_List2)) {
    Species_M2 <- c(Species_M2, length(Bin_Dino_List2[[k]]))
}
Species_M3 <- c() # METHOD 3
for (t in 1:length(av_species)) {
    Species_M3 <- c(Species_M3, mean(av_species[[t]]))
}

plot_bins <- unlist(lapply(binlist, mean))

DPSK_M1_Max <- Species_M1/Area_M1_Max
DPSK_M2_Max <- Species_M2/Area_M2_Max
DPSK_M3_Max <- Species_M3/Area_M3_Max

plot(plot_bins,DPSK_M2_Max, xlim=rev(range.default(plot_bins)), type = 'l') 

plot(Species_M1, type = 'l') # plot, obviously. 
lines(Species_M2, col="red")
lines(Species_M3, col="blue")

plot(Area_M1_Max, type = 'l', col="red")
lines(Area_M1_Min, type = 'l', col='blue')
lines(Area_M2_Max, col="red")
lines(Area_M3_Max, pch=22, col="blue")
  
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)
plot(Species_M1, type = 'l', axes=FALSE, xlab="", ylab="", col = "black")
axis(2, ylim=c(0,1),col="black",las=1)
mtext("Raw Dinosaur Diversity",side=2,line=2.5)
box()
par(new=TRUE)
plot(Area_M1_Max, type = 'l', axes=FALSE, xlab="", ylab="", col = "black", lty=2)
lines(Area_M1_Min, col="red", lty=2)
mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
axis(4, ylim=c(0,7000), col="black",col.axis="red",las=1)


lines(Area_M3_Max, col="blue", lty=2)



## Plot the second plot and put axis scale on right
plot(time, cell.density, pch=15,  xlab="", ylab="", ylim=c(0,7000), 
    axes=FALSE, type="b", col="red")
## a little farther out (line=4) to make room for labels
mtext("Outcrop Area (square km)",side=4,col="black",line=4) 
axis(4, ylim=c(0,7000), col="black",col.axis="red",las=1)

## Draw the time axis
axis(1,at = bins,7)
mtext("Time (Hours)",side=1,col="black",line=2.5)  

## Add Legend
legend("topleft",legend=c("Beta Gal","Cell Density"),
  text.col=c("black","red"),pch=c(16,15),col=c("black","red"))
# NOTE: need to add something in to count up DPSK as well, not just standard diversity.