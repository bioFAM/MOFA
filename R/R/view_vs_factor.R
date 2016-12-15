
##################################################
## Function to generate the view vs factor plot ##
##################################################

library(pheatmap)
library(RColorBrewer)

PlotViewVsFactor <- function(model, color=colorRampPalette(c("grey100", "grey0"))(100), title="", cluster_cols=F, cluster_rows=F, show_rownames=T, show_colnames=T,
                           legend=T, treeheight_row=20, treeheight_col=20, fontsize_row=20, fontsize_col=20, cellheight=NA, cellwidth=NA, outfile=NA) {
  
  # Calculate proportion of residual variation explained by each factor in each view
  prvar_mk <- CalculateProportionResidualVariance(model, plot=F)

  p <- pheatmap(prvar_mk, 
           border_color="black", main=title, color=color,
           cluster_cols=cluster_cols, cluster_rows=cluster_rows, 
           show_rownames=show_rownames, show_colnames=show_colnames,
           legend=legend, treeheight_row=treeheight_row, treeheight_col=treeheight_col,
           fontsize_row=fontsize_row, fontsize_col=fontsize_col, cellheight=cellheight, 
           filename=outfile)
  return(p)
}

# setGeneric(name="plotViewVsFactor", def=function(object) { standardGeneric("plotViewVsFactor") })
# setMethod("plotViewVsFactor", signature("GFATrainedModel"),
#           function(object, color=colorRampPalette(c("grey100", "grey0"))(100), title="", cluster_cols=T, cluster_rows=F, show_rownames=T, show_colnames=T,
#                    legend=T, treeheight_row=20, treeheight_col=20, fontsize_row=20, fontsize_col=20, cellheight=NA, cellwidth=NA, outfile=NA) {
#             GFA_plotViewVsFactor(object, color, title, cluster_cols, cluster_rows, show_rownames, show_colnames,
#             legend, treeheight_row, treeheight_col, fontsize_row, fontsize_col, cellheight, cellwidth, outfile)
#           })
            