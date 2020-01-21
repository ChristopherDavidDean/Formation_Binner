# Formation_Binner

## Authors

* **Christopher D. Dean** - [ChristopherDavidDean](https://github.com/ChristopherDavidDean)
* Alfio A. Chiarenza - [AlfioAlessandroChiarenza](https://github.com/AlfioAlessandroChiarenza)
* Susannah Maidment

---
## Introduction
### What is it?
Formation Binner consists of script to generate time bins based on regional geology and subsequently produce raw and subsampled diversity estimates based on both original and newly generated bins.
### Why should I use it?
Formation Binner can:
* Increase the resolution of diversity studies within regional contexts.
* Provide the opportunity to easily compare multiple temporal frameworks and observe how they affect the shape of the diversity curve.
* Show the distributions of formations and taxic diversity through time and space, giving the user an idea of:
  * The spatial architecture behind the distribution of formations, and the fossil record they subsequently provide;
  * Whether their data is of suitable quality to answer their chosen research question.

---
## Requirements

### Datasets

Formation Binner requires two initial datasets:
* A cleaned PBDB download of selected occurrences.
* A dataset of chronostratigraphic Formation information, including start and ends dates of all formations that appear within the occurrence dataset. 

### R Packages

The following R packages are required to run the script for Formation Binner:

```
library(pbapply)
library(dplyr)
library(tidyverse)
library(divDyn)
library(rowr)
library(matrixStats)
library(beepr)
library(iNEXT)
library(reshape2)
library(RColorBrewer)
```
---
## How it works

![Schematic diagram of how Formation_Binner works.](https://github.com/ChristopherDavidDean/Formation_Binner/blob/master/Images/Fig.%201%20(v2.1).png)

Formation binner is made up of a series of functions, which are stored within the 0_Functions_Form_Binner.R file and loaded through the 1_Setup_Formation_Binning.R file. The functions take data in the form of a .csv file of occurrences downloaded from the [Paleobiology Database](www.paleobiodb.org) and a .csv file of formations, including the neccessary fields Formation, max_age and min_age. 

Before running the main functions, a few things need to be setup. First, run `data(stages)` to load stage data from the package `divDyn`, and then select only the rows which contain the stages of interest for your study. It is also neccessary to run the line `Form_list <- split(testoccs, testoccs$formation)` for later use, as well as setting up your chosen quorums for running the SQS analysis and your binlimits which define the resolution that you are intending to work at (these settings can be founnd in the 1_Setup_Formation_Binning.R file). 

Once both occurrence and formation data have been cleaned, the first step is to run either the Scoring_Grid_1 or Scoring_Grid_2 function **(A)**, inputting your formations (the difference between the two is discussed below). Scoring grid repeatedly checks the suitability of drawing a time bin boundary in 0.01 Ma intervals, with the aim of minimizing the number of formations that are split by bins. This is accomplished by assessing the chronologic position of each formation relative to the proposed time bin boundary. If the formation does not cross the boundary, it gets an automatic maximum binning score of 100. If the formation crosses the boundary, the script assesses the smallest percentage of the duration of the formation that sits either side of that boundary, and reduces the binning score by that appropriate percentage. For example, if a formation spanned from 80 to 90 Ma and the proposed boundary was 81 Ma, 10% of the formation would fall over this boundary, and so the binning score for that formation for that boundary would be reduced by 10%, equaling a final score of 90. Scores for each formation are recorded in a score grid **(B)**, which is outputted as `score_grid` in the global environment. The vector of all possible bins is also outputted as `allbins` in the global environment. 

Next, the function `newBins` is run with the variables `score_grid`, `formations`, `bin_limits`, `allbins`, and `stages`. `newBins` generates the mean score across all formations for each increment of time that intends to reflect the suitability to draw a boundary at that point. Mean scores are then compared throughout the time series in user defined windows of time as specified in the `bin_limits` argument **(C)**, and the maximum binning values within each are outputted as bins **(D)**. if the argument `smallamalg` is set to `TRUE`, then in the event that a bin is under 0.5 Ma in duration, the script will amalgamate the bin equally into the bins above and below it, whilst generating a warning for the user. `newBins` also draws a plot of the mean score through time for the user, and outputs `form_bins`, a vector of bin ages, and `binlist`, a dataframe of bin maximum, minimum and midpoint ages. 

Once bins have been generated, you can use either `FormBin_M1`, `FormBin_M2` or `FormBin_M3` to bin occurrences and generate plots of raw diversity, number of collections, Good's u and SQS diversity through time. These functions require `formations`, `binlist`, `Form_list` and `Quorum`. 

---

## Files
### 0_Functions_Form_Binner
This file contains all the necessary functions to generate formation bins and assess diversity within those bins through time. It contains the following functions:

`ipak(pkg)` - Function to install packages. Read in character vector of any packages required.

`GetInfo(formations)` - Retrieves basic information about inputted formations, including mean age and range.

`Scoring_Grid_1(formations, res=0.01)` - Creates a scoring grid using info from all formations.

`Scoring_Grid_2(formations, res=0.01)` - Create a scoring grid ignoring formations with length longer than mean formation length. In his way, long ranging formations don't bias the creation of bins, especially when they appear during the same time interval.

`plotMaker(rel_data, binlist, ulabel)` - Generates plots through time with user inputted data and Formation_Bins, whilst providing traditional stage data for comparison. NOTE: xlim is specified to fit the chosen time window of this study - as such, this would
have to be adjusted if other data were to be used.

`overlap_counter(score_grid)` - Generates a quick high resolution plot of the number of formations present through geologic time. Will get different results depending on whether `Scoring_Grid_1` or `Scoring_Grid_2` was used to generate `score_grid`.

`newBins(score_grid, formations, bin_limits, allbins, stages, smallamalg = TRUE)` - Looks at the previously generated score_grid and generates appropriate new bins based on those scores. Boundaries are outputted as a list (binlist). When `smallamalg` is set to `TRUE`, if bins are shorter than 0.5 Ma they are amalgamated into the bins above and below, and a warning is produced for the user. 

`FormationGraph(formations, form_bins, stages, score_grid_2 = FALSE, draw_by = "Lat", Col = "None", legend = TRUE, STAGE = FALSE)` - Shows what formations look like through time in comparison to Stages and new Bins, with a number of choices to personalise output. When `score_grid_2` is set to `TRUE`, long ranging formations are greyed out. `draw_by` can be set to `Lat`, `Max_age` or `Number`, which determines the y axis for the graph (either arranging formations by latitude, maximum age or artibrary number respectively). `Col` can be set to `Diversity`, `Occurrences` or `None`, which determines whether the colour of the formations is arranged by their taxic diversity, number of occurrences, or not coloured. `Legend` allows for a colour legend to be added to the graph. `STAGE` allows for grey block to conform to traditional stages rather than new formation bins. 

`FormBin_M1(formations, binlist, Form_list, Quorum)` - Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to all bins that they occur in. Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen Quorum levels. 

`FormBin_M2(formations, binlist, Form_list, Quorum)` - Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to bins that the majority of the formation occurs in. Formations which have lengths greater than 2x the length of a bin are ignored. Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen Quorum levels. 

`FormBin_M3(formations, binlist, Form_list, times=10, Quorum, run_SQS = TRUE)` - Uses the generated boundaries from Bins to assign user specified occurrences (Form_list) and formations to all bins that they occur in, based on the percentage of the formation that sits within that bin. Occurrences are selected at random from the formation list and not replaced. The test is repeated according to user specified number of runs. Produces graphs showing raw diversity, number of collections, Good's U and SQS results at chosen Quorum levels. `run_SQS` can be set to `FALSE` to skip generating SQS results, for faster speed.

### 1_Setup_Formation_Binning
This file is set up for running binning by formation. It contains the following sections:

`DATA SETUP` - This section contains script for loading functions, packages and data, as well as cleaning and setting up data ready for running formation binning tests.

`RUNNING TESTS` - This section contains script for binning by formation, visualising formation data, running the different diversity methods and testing how different resolutions of formation_binning affect diversity results. 

### 2_Non_Formation_Binning
This file contains script for producing diversity curves using traditional bins, either at stage or substage level. Binning occurrences is done using either traditional PBDB ages or updated formation ages which were used for formation binner to compare their relative accuracy and how this influences diversity curves. 

---
