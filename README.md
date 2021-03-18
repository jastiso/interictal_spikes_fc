# Changes in functional connectivity associated with interictal epileptiform discharges (IEDs)

*Note*: This project is still in progress, and the structure and content of analyses is likely to change

Code for project investigating how functional connectivity changes in the context of IEDs. As a neuroscientist who uses a lot of iEEG data, it was never quite clear to me how to handle the presence of IED artifcats in data. This project was started in part to make those preprocessing decisions more clear

Special thanks to collaborators Lorenzo Cacialgi, Peter Hadar, Kathryn Davis, Timothy Lucas, and Danielle Bassett.

# Packages
## MATLAB (R2018A)
- fieldtrip https://www.fieldtriptoolbox.org/
- ARfit https://www.mathworks.com/matlabcentral/fileexchange/174-arfit

## Python (v3.6, all availible with pip)
- sys
- pandas
- seaborn >= v0.9.0
- matplotlib
- itertools
- statsmodels
- sklearn
- pallettable

## R (3.4)
- ggplot2
- dplyr
- lmerTest
- lm.beta
- RColorBrewer
- wesandersen
- nationalparkcolors (https://github.com/katiejolly/nationalparkcolors)
- rjson
- reticulate
- MASS
- car
- coin
- gridExtra
- stringr
- lmPerm
- glmnet
- gglasso

# Data
Most data is availible in the public RAM release http://memory.psych.upenn.edu/RAM. The only exception is additional clinical information that is not shared publically due to privacy concerns
1. Get permission to access the data from the link above, and download all three releases from the private box link.
2. Put all releases in a folder called RAM.
3. Use this RAM folder as your top_dir in all scripts

# Code
All files references here are in the "script" folder
## Preprocessing
1. preproc_wrapper_0.m
  - Preprocesses all data. This includes concatenation of individual contact files into matrices, filtering, channel rejection, resampling, rereferencing, epoching, demeaning, and detrending. After this script, all data will be in the form of a fieldtrip compatible MATLAB struct
  - This script is set up to run in parrallel for different subjects
  - All error messages are saved from this analysis
2. spike_wrapper_1.m
  - Automatically detects IEDs, by one of two methods, either Janca et. al ("") or Delphos Detector ("_delphos")
  - This will save another stuct with spike times and locations, as well as all errors
3. temporal_artifact_2.m
  - Automatically detects sharp transients and periods of flatlining
  - This saves a vector that indicates which time points have each artifcat repsectively
  - This script generates plots that will show you the data being flagged
4. reject_subj_3.m
  - Use the number of artifacts and PSD shape to reject entire datasets from further analysis
  - this script generates plots of the data being rejected
5. functional_connectivity.m
  - calculates 5 functional connectivity netrics in 5 bands and broadband data.
  - saves large CSV files giving the strength of each measure in each window for each contact
  - this script is the bottleneck of the whole pipeline. It can be parellezed on a dual or higher core machine as is
6. fit_model_5.Rmd
  - takes data from the previous steps, and calculates regression coefficients that describe how much connectivity strength changes in the context of IEDs

## Analysis
7. analysis_6.ipynb
  - includes EDA, and data frame formatting for most anayses in the paper
  - most figures are made in this script

## Statistics
8. stats_7.Rmd
  - collects effect sizes, degrees of freedom, and p-values for all statistical tests
  - also tests the normality of data distributions to justify the use of non-parametric tests

All custom functions are included in the "functions" folder. Additional scripts are for supplemental analyses
