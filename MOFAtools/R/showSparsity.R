
#####################################
## Functions to visualise sparsity ##
#####################################


#' @title showSparsity: show the sparsity of each factor in each view
#' @name showSparsity
#' @description Function to visualize the sparsity of a given factor(s) for a given view(s)
#' @param model a MOFA model
#' @param views character vector of view names (default: "all")
#' @param factors character vector of factor names (default: "all")
#' @details fill this
#' @return fill this
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal 
#' @importFrom grDevices colorRampPalette
#' @export
showSparsity <- function(model, views = "all", factors = "all", main = NULL) {
  
  # Sanity checks
  if (class(model) != "MOFAmodel") stop("'model' has to be an instance of MOFAmodel")

  # Define views
  if (paste0(views,collapse="")=="all") { 
    views <- viewNames(model) 
  } else {
    stopifnot(all(views %in% viewNames(model)))  
  }
  
  # Define factors
  factors <- as.character(factors)
  if (paste0(factors,collapse="")=="all") { 
    factors <- factorNames(model) 
  } else {
    stopifnot(all(factors %in% factorNames(model)))  
  }
  
  # Collect relevant data
  theta <- do.call("rbind",lapply(getExpectations(model,"Theta","E"), function(x) x[1,]))
  if (model@ModelOpts$learnMean) { theta <- theta[,-1]; factors <- factors[factors != "intercept"] }
  theta <- theta[views,factors]
  
  # Heatmap
  col <- colorRampPalette(c("gray97","darkblue"))(n=100)
  if(is.null(main)) main <- "Overview of sparsity levels per factor and view \n (0=sparse, 1=nonsparse)"
  # pheatmap::pheatmap(t(theta), main=main, cluster_rows = F, cluster_cols = F, col=col, ...)
  pheatmap::pheatmap(t(theta), main=main, cluster_rows = F, cluster_cols = F, col=col)
}
  