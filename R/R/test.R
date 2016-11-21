

setwd("/Users/ricard/git/scGFA/R/R")

source("loadModel.R")
source("GFAobject.R")

infolder="/tmp/test/model"
model = loadModelfromNpy(infolder)
infolder <- "/Users/ricard/git/scGFA/scMT/tmp/opts"
opts <- loadTrainingOpts(infolder)
infolder <- "/Users/ricard/git/scGFA/scMT/tmp/stats"
stats <- loadTrainingStats(infolder)

gfa <- new("GFATrainedModel", TrainStats=stats, TrainOpts=opts, Expectations=model)
getTrainOpts(gfa)
getTrainStats(gfa)
