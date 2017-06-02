
#' @title showWeights: plot weights for a certain factor and view
#' @name showWeights
#' @description Function to visualize weights in a dotplot highlighting top positive and neagtive featurs
#' @param model a fitted MOFA model
#' @param view name of view from which to get the corresponding weights
#' @param factor factor to plot weights for
#' @param ntop number of highest positive weights to label
#' @param ntail number of lowest negative weights to label
#' @details fill this
#' @return a dataframe with of d_m rows (feautres) with weights per feature on the sepcified factor and view
#' @import ggplot2 ggrepel
#' @export

showWeights <- function(model, view, factor, ntop = 5, ntail =5, manual = NULL) {
  
  # Sanity checks
  if (class(model) != "MOFAmodel")
    stop("'model' has to be an instance of MOFAmodel")
  stopifnot(all(view %in% viewNames(model)))  
  
  if(is.null(factorNames)) factorNames(model) <- 1:model@Dimensions$K
  stopifnot((factor %in% factorNames(model)))

  W <- getExpectations(model,"SW","E")[[view]][,factor]
  
  df_W <- data.frame(loading=W, FeatureName = names(W))
  df_W$imp <- F
  if(ntop>0) df_W$imp[df_W$loading >= sort(df_W$loading, decreasing = T)[ntop]] <- T
  if(ntail>0) df_W$imp[df_W$loading <= sort(df_W$loading, decreasing = F)[ntail]] <- T
  if(!is.null(manual)) df_W$imp[df_W$FeatureName %in% manual] <-T
  
  df_W %<>% arrange(loading)
  df_W$FeatureName <- factor(df_W$FeatureName, levels=df_W$FeatureName)
  
  gg_W <- ggplot2::ggplot(df_W, aes(x=FeatureName, y=loading, col=imp)) + geom_point()  +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ggrepel::geom_text_repel(data = filter(df_W, imp), aes(label = FeatureName),
                    segment.alpha=0.2, box.padding = unit(0.5, "lines")) +scale_color_manual(values=c("black", "red")) +
    ggtitle(paste("Weigths on LF", factor, "on view", view)) +
    guides(color=F) 
  print(gg_W)
  
  return(df_W)
}

