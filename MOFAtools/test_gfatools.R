
library(MOFAtools)

# file = "/Users/ricard/data/CLL/out/robustness/trained_model_batch2_1.hdf5"
# file = "/Users/ricard/data/CLL/out/imputation/model_0.05.hdf5"
file = "/Users/ricard/test.hdf5"
model <- loadModel(file)

# FactorsCorPlot(model)
p <- CalculateVariance_Views(model, views="all", factors="all")
# ViewFactorPlot(model)

FactorsCorPlot(model)
CorrplotLFvsallPC(model, method="svd", noPCs=5)

## training curves ##
trainCurve(model, statistic="activeK", xlabel=" ", ylabel="")
trainCurve(model, statistic="elbo", xlabel=" ", ylabel="")


## test mutation ##

# asd < /Users/ricard/data/CLL/Pre-Analysis/views/assembled/model_EDMGI/Covariates.txt
load("/Users/ricard/data/CLL/Pre-Analysis/views/assembled/CLLOmicsList.RData")

samples <- sampleNames(model)
metadata <- OmicsList$metadata_view[samples,]
cor(metadata$IGHV, getExpectations(model,"Z","E"), use="complete.obs")

scatterPlot(model, 2, 3, colour_by=metadata$IGHV)
            
mut <- OmicsList$mutation_view[samples,]

cor(mut, getExpectations(model,"Z","E"), use="complete.obs") %>% View

## Test GSEA ##

# # Simulate a test data set
# N = 100
# D = 50
# K = 10
# M = 20
# 
# data <- matrix(rnorm(N*D),N,D)
# W <- matrix(rnorm(D*K),D,K)
# Z <- matrix(rnorm(N*K),N,K)
# sets <- matrix(sample(c(0,1),size=M*D,replace=T,prob=c(0.8,0.2)),M,D)
# local.statistic="loading"
# transformation="abs.value"
# global.statistic="mean.diff"
# statistical.test="parametric"
# nperm=NA
# factor.indexes=c(1,2,3)
# 

# Read binary matrix of annotations


b <- readRDS("/Users/ricard/data/reactome/v59/homo_sapiens/out/human_reactome.rds")
view <- "mRNA"
factor.indexes="all"
gene.sets=b
local.statistic="z"
transformation="abs.value"
global.statistic="mean.diff"
statistical.test="permutation"
nperm=3

p <- GSEA(model, view, factor.indexes, gene.sets=b, local.statistic=local.statistic, transformation=transformation, 
          global.statistic=global.statistic, statistical.test=statistical.test, nperm=nperm, min_size=15)


##################

sigmoid <- function(x) 1/(1+exp(-x))
Ypred <- sigmoid(t(model@Expectations$SW$mut$E %*% t(model@Expectations$Z$E)))
# Ypred <- t(model@Expectations$SW$mut$E %*% t(model@Expectations$Z$E))
Yobs <- model@TrainData$mut

# impute and pca
impute <- function(d, margin) {
  if (margin == 1)
    means <- rowMeans(d, na.rm=T)
  else if (margin == 2)
    means <- colMeans(d, na.rm=T)
  else
    stop("Margin has to be either 1 (rows) or 2 (cols)")
  
  if (any(is.na(means))) {
    stop('Insufficient data for mean imputation!')
  }
  
  for (i in 1:length(means)) {
    if (margin == 1)
      d[i,is.na(d[i,])] <- means[i]
    else if (margin == 2)
      d[is.na(d[,i]), i] <- means[i]
  }
  return (d)
}

Yobs_imput <- impute(Yobs,margin=2)
p <- prcomp(Yobs_imput)

asd <- cor(model@Expectations$Z$E,p$x)

