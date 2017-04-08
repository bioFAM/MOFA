
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

showWeights <- function(model, viewnm, showHeatmap =T, main=NULL,...) {
stopifnot(viewnm %in% viewNames(model)) 
  
Z<-model@Expectations$Z$E

#get weights and appropriate names
weights<-model@Expectations$SW[[viewnm]]$E

if(!is.null(colnames(Z))) namesLF<-colnames(Z) else namesLF<- 1:ncol(Z)
colnames(weights)<-namesLF
rownames(weights)<-colnames(model@TrainData[[viewnm]])

#plot heatmap
if(showHeatmap){
if(is.null(main)) main <- paste("Weights of Latent Factors on", viewnm)
pheatmap(weights, main=main,...)
}

#return weight matrix
return(weights)
}

