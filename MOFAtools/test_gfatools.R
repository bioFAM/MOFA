
library(MOFAtools)

# file = "/Users/ricard/git/britta/scGFA/surv_expr/model0.hdf5"
file = "/Users/ricard/git/gastrulation/join/MOFA/out/model0.hdf5"
model <- loadModel(file)
model

FactorsCorPlot(model)
CalculateProportionVariance(model,plot=T)
ViewFactorPlot(model)


N <- model@Dimensions[["N"]]
colour_by <- c(rep(T,N/2),rep(F,N/2),FALSE)
shape_by <- c(rep(T,N/2),rep(F,N/2),FALSE)
scatterPlot(model, 1, 2, title="", titlesize=16, xlabel="Latent variable 1", ylabel="Latent variable 2",
dotsize=2.5, colour_by=colour_by, shape_by=shape_by, xlim_down=NULL, xlim_up=NULL, ylim_down=NULL, ylim_up=NULL)

trainCurve(model, statistic="activeK", xlabel=" ", ylabel="")

