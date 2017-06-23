
########################################
## Functions to visualise the weights ##
########################################


#' @title showWeightHeatmap: show heatmap of the factor loadings in a specific view
#' @name showWeightHeatmap
#' @description Function to visualize the weights for each feature in each factor for a given view
#' @param model a MOFA model
#' @param view view name
#' @param features name of the features to include (default: "all")
#' @param factors name of the factors to include (default: "all")
#' @param RemoveIntercept boolean wether to include the intercept factor if present in the model (by default not included)
#' @param main asd 
#' @param color asd 
#' @param breaks asd 
#' @param ... further arguments that can be passed to pheatmap
#' @details fill this
#' @return fill this
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal 
#' @importFrom grDevices colorRampPalette
#' @export
showWeightHeatmap <- function(model, view, features = "all", factors = "all", RemoveIntercept = T, main=NULL, color=NULL, breaks=NULL, ...) {
  
  # Sanity checks
  if (class(model) != "MOFAmodel") stop("'model' has to be an instance of MOFAmodel")
  stopifnot(all(view %in% viewNames(model)))  
  
  # Define factors
  factors <- as.character(factors)
  if (paste0(factors,collapse="")=="all") { 
    factors <- factorNames(model) 
    # if(is.null(factors)) factors <- 1:ncol(getExpectations(model,"Z","E"))
  } else {
    stopifnot(all(factors %in% factorNames(model)))  
  }
  
  # Remove intercept factor
  if(model@ModelOpts$learnMean==T & RemoveIntercept) { factors <- factors[factors != factorNames(model)[1]] }
  
  # Define features
  if (paste(features,collapse="")=="all") { 
    features <- featureNames(model)[[view]]
  } else {
    stopifnot(all(features %in% featureNames(model)[[view]]))  
  }

  # Get relevant data
  W <- getExpectations(model,"SW","E")[[view]][features,factors]

  # Set title
  if (is.null(main)) { main <- paste("Loadings of Latent Factors on", view) }
  
  # set colors and breaks if not specified
  if(is.null(color) & is.null(breaks)){
    palLength <- 100
    minV <- min(W)
    maxV <- max(W)
    color <- colorRampPalette(colors=c("black", "blue", "white", "orange","red"))(palLength)
    breaks <- c(seq(minV, 0, length.out=ceiling(palLength/2) + 1), 
                seq(maxV/palLength,maxV, length.out=floor(palLength/2)))
  }
  
  # Plot heatmap
  pheatmap::pheatmap(W, fontsize=15, ...)
}



#' @title showWeights: visualise the loadings for a certain factor in a given view
#' @name showWeights
#' @description Function to visualize weights in a dotplot highlighting features with highest loadings
#' @param model a fitted MOFA model
#' @param view name of view from which to get the corresponding weights
#' @param factor name of the factor
#' @param ntop number of highest positive weights to label
#' @param ntail number of lowest negative weights to label
#' @param manual feature names which are labeled in the plot
#' @param th threshold on the absolute value beyond which weigths are labeled
#' @param ntotal number of features with highest absolute value to label
#' @param FeatureGroups  a named vector of same lengh as number of feautres annotating feature to groups that 
#'  will be colored in the dotplot, names have to be in accordance with feature names
#' @param FeatureGroupsCols  a named vector speciying for each elemnt in FeautreGroups a color to use for plotting
#' @param showFeatureNames boolean wether to show all FeatureNames on the x-axis (default false)
#' @param main plot title
#' @details fill this
#' @import dplyr ggplot2 ggrepel magrittr
#' @export

showWeights <- function(model, view, factor, ntop = 0, ntail = 0, manual = NULL, th = NULL, ntotal= 10, 
                        FeatureGroups = NULL, FeatureGroupsCols=NULL, showFeatureNames = F, main = NULL) {
  
  # Sanity checks
  if (class(model) != "MOFAmodel")
    stop("'model' has to be an instance of MOFAmodel")
  stopifnot(all(view %in% viewNames(model)))  
  
  if(is.null(factorNames)) factorNames(model) <- 1:model@Dimensions$K
  factor <- as.character(factor)
  stopifnot((factor %in% factorNames(model)))

  W <- getExpectations(model,"SW","E")[[view]][,factor]
  
  df_W <- data.frame(loading=W, FeatureName = names(W))
  df_W$imp <- F
  if(ntop>0) df_W$imp[df_W$loading >= sort(df_W$loading, decreasing = T)[ntop]] <- T
  if(ntail>0) df_W$imp[df_W$loading <= sort(df_W$loading, decreasing = F)[ntail]] <- T
  if(!is.null(manual)) df_W$imp[df_W$FeatureName %in% manual] <-T
  if(!is.null(th)) df_W$imp[abs(df_W$loading) >= th] <- T
  if(ntotal > 0) df_W$imp[abs(df_W$loading) >= sort(abs(df_W$loading), decreasing = T)[ntotal]] <- T
  
  df_W %<>% arrange(loading)
  df_W$FeatureName <- factor(df_W$FeatureName, levels=df_W$FeatureName)
  
  #color groups of features
  if(!is.null(FeatureGroups)) {
    #sanity check
    if(length(FeatureGroups) != length(W)) stop("Length of FeatureGroup and number of features in the view do not agree")
    if(!setequal(names(FeatureGroups),as.character(df_W$FeatureName))) stop("Names of FeatureGroup and names of features in the view do not agree")
    
    df_W$groups <- as.factor(FeatureGroups[as.character(df_W$FeatureName)] )
    if(is.null(FeatureGroupsCols)){
      FeatureGroupsCols <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
      if(length(unique(df_W$groups)) > 29) FeatureGroupsCols <- colors(length(unique(df_W$groups)))
      FeatureGroupsCols <- FeatureGroupsCols[1:length(unique(df_W$groups))]
      names(FeatureGroupsCols) <- unique(FeatureGroups)
    }}
  else df_W$groups <- T
  
  if(is.null(main)) main <- paste("Weigths on LF", factor, "on view", view)

  gg_W <- ggplot2::ggplot(df_W, aes(x=FeatureName, y=loading, col= groups)) + geom_point() +
    ggrepel::geom_text_repel(data = filter(df_W, imp), aes(label = FeatureName, col = groups),
                             segment.alpha=0.2, box.padding = unit(0.5, "lines"), show.legend= F) +
    ggtitle(main)
  
  if(!is.null(FeatureGroups)) gg_W <- gg_W + scale_color_manual(values=FeatureGroupsCols)
    else gg_W <- gg_W + guides(col=F) + scale_color_manual(values=c("black", "red"))  +
                geom_point(aes(x=FeatureName, y=loading, col=imp)) 
  
    if(showFeatureNames)
      gg_W <- gg_W + theme(axis.title.x=element_blank(),
                       axis.text.x=element_text(angle = 90, hjust = 1),
                       axis.ticks.x=element_blank(),
                       panel.background =element_rect(fill="white"))
    else
      gg_W <- gg_W + theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank(),
                           panel.background =element_rect(fill="white"))
  
  print(gg_W)
}

