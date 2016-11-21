# Function to generate the view vs factor plot

library(pheatmap)
library(RColorBrewer)

setwd("/Users/ricard/git/scGFA/R/R")
# source("GFAobject.R")
# source("loadModel.R")

# Create GFA object
basedir="/tmp/test"
model = loadModelfromNpy(paste(basedir,"/model",sep=""))
opts <- loadTrainingOpts(paste(basedir,"/opts",sep=""))
stats <- loadTrainingStats(paste(basedir,"/stats",sep=""))
data <- loadTrainingData(paste(basedir,"/data",sep=""))

gfa <- new("GFATrainedModel", TrainStats=stats, TrainOpts=opts, TrainData=data)

# Spike and Slab expectations
gfa@Expectations$tau <- model$tau$E
gfa@Expectations$S <- model$SW$ES
gfa@Expectations$W <- model$SW$EW
gfa@Expectations$Z <- model$Z$E
gfa@Expectations$alpha <- do.call(rbind,model$alpha$E)
gfa@Expectations$Zeta <- model$Zeta$E
gfa@Expectations$Y <- model$Y$E



GFA_plotViewVsFactor <- function(object, color=colorRampPalette(c("grey100", "grey0"))(100), title="", cluster_cols=F, cluster_rows=F, show_rownames=T, show_colnames=T,
                           legend=T, treeheight_row=20, treeheight_col=20, fontsize_row=20, fontsize_col=20, cellheight=NA, cellwidth=NA, outfile=NA) {
  
  # Calculate proportion of residual variation explained by each factor in each view
  # HOW DOES THIS WORK WITH PSEUDODATA NODES?
  
  Y <- gfa@Expectations$Y
  # Y <- gfa@TrainData
  residual_var <- sapply(1:length(Y), function(m) sum(apply(Y[[m]],2,var,na.rm=T) - 1/model$tau$E[[m]]))
  total_var <- sapply(1:model_opts$M, function(m) sum(apply(Y[[m]],2,var)))
  
  mat = object@Expectations$alpha
  p <- pheatmap(t(mat), border_color="black", main=title, color=color,
           cluster_cols=cluster_cols, cluster_rows=cluster_rows, show_rownames=show_rownames, show_colnames=show_colnames,
           legend=legend, treeheight_row=treeheight_row, treeheight_col=treeheight_col,
           fontsize_row=fontsize_row, fontsize_col=fontsize_col, cellheight=cellheight, filename=outfile)
  return(p)
}

# p <- GFA_plotViewVsFactor(gfa)
# print(p)

setGeneric(name="plotViewVsFactor", def=function(object) { standardGeneric("plotViewVsFactor") })
setMethod("plotViewVsFactor", signature("GFATrainedModel"),
          function(object, color=colorRampPalette(c("grey100", "grey0"))(100), title="", cluster_cols=T, cluster_rows=F, show_rownames=T, show_colnames=T,
                   legend=T, treeheight_row=20, treeheight_col=20, fontsize_row=20, fontsize_col=20, cellheight=NA, cellwidth=NA, outfile=NA) {
            GFA_plotViewVsFactor(object, color, title, cluster_cols, cluster_rows, show_rownames, show_colnames,
            legend, treeheight_row, treeheight_col, fontsize_row, fontsize_col, cellheight, cellwidth, outfile)
          })
            