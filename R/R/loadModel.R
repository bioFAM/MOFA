
#############################################
## Functions to load a trained scGFA model ##
#############################################

# scGFA is trained using Python, but post-analysis and plotting are simpler in R.
# This script provides functions to load the trained model into R

# The trained model is saved in an HDF5 file with the following groups:
#   expectations
#   parameters
#   training_options
#   training_stats
#   data

library(rhdf5)
library(stringr)

setwd("/Users/ricard/git/scGFA/R/R")
source("GFAobject.R")

loadModel <- function(file) {
  
  # Load expectations and parameters
  expectations <- h5read(file,"expectations")
  parameters <- h5read(file,"parameters")
  
  # Load training statistics
  training_stats <- h5read(file,"training_stats")
  colnames(training_stats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  
  # Load training options
  training_opts <- as.list(h5read(file,"training_opts", read.attributes=T))
  
  # Load training data
  # training_data <- h5read(file,"training_data", read.attributes=T)
  # for (view in names(training_data)) {
  #   tmp <- h5read(file,str_c("training_data/",view), read.attributes=T)
  #   rownames(training_data[[view]]) <- attr(tmp,"rownames")
  #   colnames(training_data[[view]]) <- attr(tmp,"colnames")
  # }
  
  # Load training data
  data <- h5read(file,"data")
  featuredata <- h5read(file,"features")
  sampledata <- h5read(file,"samples")
  for (m in names(data)) {
    rownames(data[[m]]) <- sampledata
    colnames(data[[m]]) <- featuredata[[m]]
  }
  
  # Define GFATrainedModel
  gfa <- new("GFATrainedModel", Data=data, TrainStats=training_stats, TrainOpts=training_opts, Expectations=expectations, Parameters=parameters)
  
  return(gfa)
}

