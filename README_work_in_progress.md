# Formation_Binner

## Authors

* **Christopher D. Dean** - [ChristopherDavidDean](https://github.com/ChristopherDavidDean)
* Alfio A. Chiarenza
* Susannah Maidment

---
## Introduction
### What is it?
Formation Binner consists of script to generate time bins based on regional geology and subsequently produce raw and subsampled diversity estimates based on both original and newly generated bins.
### Why should I use it?
Formation Binner allows for 

---
## Requirements

### Datasets

Formation Binner requires two initial datasets:
* A cleaned PBDB download of selected occurrences.
* A dataset of chronostratigraphic Formation information, including start and ends dates of all formations that appear within the occurrence dataset. 

### R Packages

The following R packages are required to run the script for Formation Binner:

```
library(pbmcapply)
library(dplyr)
library(tidyverse)
library(divDyn)
library(rowr)
library(matrixStats)
library(beepr)
```
---
## How it works

Formation binner...

---

## Files
### 0_Functions_Form_Binner
This file contains all the necessary functions to generate formation bins and assess diversity within those bins through time. It contains the following functions:

`ipak(pkg)` - Function to install packages. Read in character vector of any packages required.

`GetInfo(formations)` - Retrieves basic information about inputted formations, including mean age and range.

`Scoring_Grid_1(formations, res=0.01)` - Creates a scoring grid using info from all formations.

`Scoring_Grid_2(formations, res=0.01)` - Create a scoring grid ignoring formations with length longer than mean formation length. In his way, long ranging formations don't 
bias the creation of bins, especially when they appear during the same time interval.

`plotMaker(rel_data, binlist, ulabel)` - Generates plots through time with user inputted data and Formation_Bins, whilst providing traditional stage 
data for comparison. NOTE: xlim is specified to fit the chosen time window of this study - as such, this would
have to be adjusted if other data were to be used.

`overlap_counter(score_grid)` - Generates a quick high resolution plot of the number of formations present through geologic time. 

`newBins(score_grid, formations, bin_limits, allbins, stages, smallamalg = TRUE)` - Looks at the previously generated score_grid and generates appropriate new bins based on those scores. Boundaries are 
outputted as a list (binlist). If bins are shorter than 0.5 Ma, they are amalgamated into the bins above and below, and 
a warning is produced.

`FormationGraph(formations, form_bins, stages)` - Shows what formations look like through time in comparison to Stages and new Bins.

`FormBin_M1(formations, binlist, Form_list, Quorum)` - Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to all bins
that they occur in. Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen
Quorum levels. 

`FormBin_M2(formations, binlist, Form_list, Quorum)` - Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to bins that
the majority of the formation occurs in. Formations which have lengths greater than 2x the length of a bin are ignored.
Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen Quorum levels. 

`FormBin_M3(formations, binlist, Form_list, times=10, Quorum)` - Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to all bins
that they occur in, based on the percentage of the formation that sits within that bin. Occurrences are selected at
random from the formation list and not replaced. The test is repeated according to user specified number of runs.
Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen Quorum levels. 


### 1_Setup_Formation_Binning

### 2_Non_Formation_Binning

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```
---
## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
