
setwd("/Users/ricard/git/scGFA/R/R")


source("loadModel.R")
file = "/tmp/test/asd.hd5"
model <- loadModel(file)

source("proportion_variance.R")
# prvar_mk <- CalculateProportionResidualVariance(model, plot=T)

# View vs Factor plot
source("view_vs_factor.R")
PlotViewVsFactor(model, color=colorRampPalette(c("grey100", "grey0"))(100), title="", cluster_cols=F, cluster_rows=F, show_rownames=T, show_colnames=T,
                 legend=T, treeheight_row=20, treeheight_col=20, fontsize_row=20, fontsize_col=20, cellheight=NA, cellwidth=NA, outfile=NA)
  

# Training statistics
source("training_statistics.R")
elbo_trainCurve(model)
activeK_trainCurve(model)

# Scatterplot
source("scatterplot.R")
scatterPlot(model, idx=1, idy=2)

# Corplot
source("corPlot.R")
FactorsCorrPlot(model)

  
# Y <- model@TrainData
e <- model@Expectations
p <- model@Parameters
data <- model@Data


e$Alpha$bernoulli$E
colMeans(abs(e$SW$bernoulli$EW))
e$Theta$bernoulli$E
abs(e$SW$bernoulli$EW)[,1]
e$SW$bernoulli$EW[e$SW$bernoulli$ES[,1] > 0.75,1]

i <- which(e$SW$bernoulli$ES[,1] > 0.75)
apply(e$Y$bernoulli$E,2,var)[i]
mean(apply(e$Y$bernoulli$E,2,var))


getTrainOpts(model)
getTrainStats(model)

m = "foo"
e$Alpha[[m]]$E
e$Theta[[m]]$E
k=1
e$Alpha[[m]]$E[k]
e$Theta[[m]]$E[k]

p$Alpha[[m]]$a[k]
p$Alpha[[m]]$b[k]

p$Theta[[m]]$a[k]
p$Theta[[m]]$b[k]

mean(e$SW[[m]]$ES[,k])
mean(e$SW[[m]]$ES[,k] > 0.75)

e$SW[[m]]$ES[,k]
plot(e$Z$E[,k])

