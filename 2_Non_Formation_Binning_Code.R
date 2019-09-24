# Running Traditional Binning Methods in SQS

library(divDyn)

# Running based on PBDB ages

#Set Up
occs <- read.csv(file = "Data/NADINOS-occs-edit.csv") # Read in occurrences
data(stages)
colnames(occs)[38] <- "new_bin"
occs$new_bin <- as.numeric(as.character(occs$new_bin))

# Binning
for(s in 1:nrow(stages)){
  for(o in 1:nrow(occs)){
    if(occs$ma_mid[o] < stages$bottom[s] && occs$ma_mid[o] > stages$top[s]){
      occs$new_bin[o] <- stages$stg[s]
    }
  }
}

# Div and SQS
bin_info <- binstat(occs, tax="occurrence.genus_name", bin="new_bin", 
                    coll = 'collection_no')
SQS <- subsample(occs,iter=100, q=0.6,tax="occurrence.genus_name", bin="new_bin", 
                 coll = 'collection_no', output="dist", type="sqs", 
                 duplicates = TRUE, useFailed = TRUE)



# Running based on updated Formation ages