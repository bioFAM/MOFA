
library(MOFAtools)

# file = "/Users/ricard/git/britta/scGFA/surv_expr/model0.hdf5"
file = "/Users/ricard/data/CLL/out/expr_nocovariates.hdf5"
model <- loadModel(file)
model

# require(gridExtra)
# p1 <- FeaturesCorPlot(model, view=1, regress_factors=NULL, legend=TRUE)
# p2 <- FeaturesCorPlot(model, view=1, regress_factors="all", legend=TRUE)
# grid.arrange(p1, p2, ncol=2)
# arrangeGrob(qp,hm, layout_matrix=lm)

FactorsCorPlot(model)
CalculateProportionVariance(model,plot=T)
ViewFactorPlot(model)


N <- model@Dimensions[["N"]]

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
view <- "expr_mRNA"
factor.indexes="all"
gene.sets=b_filt
local.statistic="loading"
transformation="abs.value"
global.statistic="mean.diff"
statistical.test="cor.adj.parametric"
nperm=10

p <- GSEA(model, view, factor.indexes, gene.sets=b, local.statistic=local.statistic, transformation=transformation, 
          global.statistic=global.statistic, statistical.test=statistical.test, nperm=nperm, min_size=15)
