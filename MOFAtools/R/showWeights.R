
#' @title showWeights: get the loadings in a specific view
#' @name showWeights
#' @description Function to extract the weights that each feautre has on a factor in the view specified and visualize it as a heatmap.
#' @param model a fitted MOFA model
#' @param viewnm name of view from which to get the corresponding weights
#' @param showHeatmap boolean, wether to plot a heatmap of the weights (default: True)
#' @param ... further arguments that can be passed to pheatmap
#' @details fill this
#' @return a weight matrix of dimension d (feautres) x k (number of latent factors)
#' @import pheatmap
#' @export

showWeights <- function(model, view, features="all", factors="all", main=NULL, ...) {
  
  # Sanity checks
  if (class(model) != "MOFAmodel")
    stop("'model' has to be an instance of MOFAmodel")
  stopifnot(all(view %in% viewNames(model)))  
    
  # Define factors
  if (factors=="all") { 
    factors <- factorNames(model) 
  } else {
    stopifnot(all(factors %in% factorNames(model)))  
  }

  # Define features
  if (features=="all") { 
    features <- featureNames(model)[[view]]
  } else {
    stopifnot(all(features %in% featureNames(model)[[view]]))  
  }
  
  # Collect expectations
  # Z <- getExpectations(model,"Z","E")[,factors]
  # if(!is.null(colnames(Z))) namesLF <- colnames(Z) else namesLF <-  1:ncol(Z)
  
  # COMMENT Britta: use featureNames and factorNames, they automatically rename all dimensions of all variables
  
  W <- getExpectations(model,"SW","E")[[view]][features,factors]
  # colnames(W) <- namesLF
  # rownames(W) <- colnames(model@TrainData[[viewnm]])
  
  # Plot heatmap
  if (is.null(main)) main <- paste("W of Latent Factors on", view)
  pheatmap::pheatmap(W, main=main,...)
}

