
setwd("/Users/ricard/git/scGFA/R/R")

source("loadModel.R")
file = "/Users/ricard/git/britta/scGFA/fullmodel300.hdf5"
model <- loadModel(file)
meta <- read.table("/Users/ricard/git/britta/processed_data/joined/metadata.txt", sep="\t", header=T)

Y <- model@Data
e <- model@Expectations
p <- model@Parameters

# prvar_mk <- CalculateProportionResidualVariance(model, plot=T)

# View vs Factor plot
source("proportion_variance.R")
source("view_vs_factor.R")
PlotViewVsFactor(model, color=colorRampPalette(c("grey100", "grey0"))(100), title="", cluster_cols=T, cluster_rows=F, show_rownames=T, show_colnames=T,
                 legend=T, treeheight_row=20, treeheight_col=20, fontsize_row=20, fontsize_col=15, cellheight=NA, cellwidth=NA, outfile=NA)

# Training statistics
source("training_statistics.R")
elbo_trainCurve(model)
# activeK_trainCurve(model)

# Scatterplot
source("scatterplot.R")

bool <- c(FALSE,TRUE)[meta$ighv]

scatterPlot(model, idx=10, idy=4, colour_by=bool, shape_by=bool )

# Look which genes are responsible for thsi variation
# Option 1: Differential expression 
# Option 2: weights
asd <- readRDS("/Users/ricard/git/britta/processed_data/expr/filt2/expr_matrix.rds")
genes <- colnames(asd)
w <- e$SW$expr$ESW[,10]
names(w) <- genes
colnames(Y$expr) <- genes

sort(abs(w))

ENSG00000075651 -> HAS BEEN ASSOCIATED WITH CANCER
mean(Y$expr[meta$ighv=="U","ENSG00000075651"], na.rm=T)
mean(Y$expr[meta$ighv=="M","ENSG00000075651"], na.rm=T)

ENSG00000205683
mean(Y$expr[meta$ighv=="U","ENSG00000205683"], na.rm=T)
mean(Y$expr[meta$ighv=="M","ENSG00000205683"], na.rm=T)

ENSG00000196263

ENSG00000205730
# Corplot
# source("corPlot.R")
# FactorsCorrPlot(model)


M = length(Y)
K = model@TrainStats$activeK[1]
alpha <- matrix(nr=K, nc=M); colnames(alpha) <- names(Y)
for (m in 1:M)
  alpha[,m] = e$Alpha[[m]]$E
View(alpha)

View(e$Z$E)
e$Theta$expr$E

head(e$Tau$expr$E)
head(1/apply(Y$expr,2,var))

d=4
hist(Y$expr[,d], xlim=c(-2.5,2.5), ylim=c(0,1), freq=F)
par(new=T)
x <- seq(-2.5,2.5,length=1000)
plot(x=x, y=dnorm(x,mean=0, sd=sqrt(1/4.072509)), type="l", lwd=1, xlim=c(-2.5,2.5), ylim=c(0,1))

1/var(Y$expr[,1])


## Correlate PC with components ##

pca <- prcomp(Y$expr, scale=F, center=T)
# pca <- prcomp(data, scale=F, center=T)
pca$pvar <- pca$sdev**2/sum(pca$sdev**2)

# scatter plot
library(ggplot2)
x <- 1
y <- 2
tmp <- data.frame(sample=rownames(pca$x), pcx=pca$x[,x], pcy=pca$x[,y])
tmp$ighv <- meta[tmp$sample,"ighv"]

ggplot(tmp, aes(x=pcx, y=pcy)) + 
  geom_point(aes(color=ighv), size=2) +
  xlab(sprintf('pc%d',x)) + ylab(sprintf('pc%d',y)) +
  theme(panel.background=element_blank())


