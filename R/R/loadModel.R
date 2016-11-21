
#############################################
## functions to load a trained scGFA model ##
#############################################

# scGFA is trained using Python, but post-analysis and plotting are simpler in R.
# This file provides functions to load the trained model into R

# Ideally, we would save the entire model in one file using pickle or something similar. However,
# We found no easy way to exchange the data with R.
# Therefore, we followed a different approach: we save the parameters and/or expectations in separate NumPy files (.npy),
# which can be read into R using the RcppCNPy package

# Functions:
# - loadModelfromNpy: load the parameters and/or expectations of a model from .npy files
# - loadTrainingStats: load the training statistics from a .txt file

library(RcppCNPy)
library(stringr)

# setwd("/Users/ricard/git/scGFA/R/R")

# source("GFAobject.R")


# Function to load training statistics
# Input: 
# Output: 
loadTrainingStats <- function(infolder) {
  # files = list.files(path=infolder,pattern="*.txt")
  # stats <- list()
  # for (file in files) {
  #   stat_name <- strsplit(file, "\\.")[[1]][1]
  #   stats[[stat_name]] <- read.table(str_c(infolder,"/",file))
  # }
  stats <- list()
  stats[["activeK"]] <- read.table(str_c(infolder,"/activeK.txt"))[,1]
  stats[["elbo"]] <- read.table(str_c(infolder,"/elbo.txt"))[,1]
  stats[["elbo_terms"]] <- read.table(str_c(infolder,"/elbo_terms.txt"), header=T)
  return(stats)
}



# Function to load training options
# Input: 
# Output: 
loadTrainingOpts <- function(infolder) {
  # files = list.files(path=infolder,pattern="*.txt")
  # stats <- list()
  # for (file in files) {
  #   stat_name <- strsplit(file, "\\.")[[1]][1]
  #   stats[[stat_name]] <- read.table(str_c(infolder,"/",file))
  # }
  tmp <- read.table(str_c(infolder,"/opts.txt"), sep=":", stringsAsFactors=F)
  opts <- as.list(tmp[,2])
  names(opts) <- tmp[,1]
  return(opts)
}

# Function to load data used for training
# Input: 
# Output: 
loadTrainingData <- function(infolder) {
  files = list.files(path=infolder,pattern="*.npy")
  data <- list()
  for (file in files) {
    view <- strsplit(file, "\\.")[[1]][1]
    data[[view]] <- npyLoad(str_c(infolder,"/",file))
  }
  return(data)
}


# Function to load model expectations
# Input: the name of the directory where the parameters or expectations are stored as follows:
#   .....
# Output: a list where the name of the items are the name of the node expectation/parameter and the value is the corresponding matrix
loadModelfromNpy <- function(infolder, only_first_moments=F) {
  data = list()
  # List all NumPy files in the directory
  files = str_replace( list.files(path=infolder,pattern="*.npy"), ".npy", "" )
  # Extract the name of the nodes (W,Z,...)
  vars = unique( sapply(files, function(f) strsplit(f, "\\_")[[1]][1]) )
  # Loop over nodes
  for (v in vars) {
    v_files = files[str_detect(files,str_c(v,"_"))]
    # Extract the name of the expectation/parameter (E1,E2,cov,mean,...)
    params = unique( sapply(v_files, function(f) strsplit(f, "\\_")[[1]][2]) )
    data[[v]] <- sapply(params, function(x) NULL)
    # Loop over expectations/parameters
    for (p in params) {
      p_files = v_files[str_detect(v_files,p)]
      views = as.numeric( unique( sapply(p_files, function(f) strsplit(f, "\\_")[[1]][3]) ) )
      if (all(is.na(views))) {
        # Single-view nodes
        data[[v]][[p]] = npyLoad(sprintf("%s/%s_%s.npy",infolder,v,p))
      } else {
        # Multi-view nodes
        data[[v]][[p]] = vector("list",length=length(views))
        # Loop over views
        for (m in views)
          data[[v]][[p]][[m]] = npyLoad(sprintf("%s/%s_%s_%s.npy",infolder,v,p,m))
      }
    }
  }
  return(data)
}

# infolder="/Users/ricard/git/scGFA/scMT/tmp/model"

# Create GFA object

# infolder="/Users/ricard/git/scGFA/scMT/tmp/model/first_moments"
# model = loadModelfromNpy(infolder)
# infolder <- "/Users/ricard/git/scGFA/scMT/tmp/opts"
# opts <- loadTrainingOpts(infolder)
# infolder <- "/Users/ricard/git/scGFA/scMT/tmp/stats"
# stats <- loadTrainingStats(infolder)
# 
# gfa <- new("GFATrainedModel", TrainStats=stats, TrainOpts=opts, Expectations=model)
# getTrainOpts(gfa)
# getTrainStats(gfa)
