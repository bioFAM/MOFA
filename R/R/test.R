

setwd("/Users/ricard/git/scGFA/R/R")

source("loadModel.R")
model <- loadModel("/tmp/test/asd.hd5")

getTrainOpts(model)
getTrainStats(model)

source("proportion_variance.R")
CalculateResidualVariance(model,active_genes_threshold=0.5, plot=T)
