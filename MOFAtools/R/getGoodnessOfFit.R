
#' @title Get a goodness of fit measure for the model in each view and for each LF
#' @name getGoodnessOfFit
#' @description Method to calculate a a goodness of fit measure for the model in each view and for each LF
#' For gaussian views the coefficient of determination is used as a measure and 
#' is calculated for the whole model as well as for individual factors on each view and feature.
#' For binary views the Brier score is used  as a measure and 
#' is calculated for the whole model as well as for individual factors on each view and feature.
#' Other families are not yet implemented.
#' The resulting measures per factor and overall are plotted as a heatmap and barplot by default.
#' NOTE: The two measure do not measure the same thing and should not be compared!
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views Views to use, default is "all"
#' @param factors Latent variables or factores to use, default is "all"
#' @param plotit boolean, wether to produce a plot (default true)
#' @details fill this
#' @return a list containing list of R2 for all gaussian views and list of BS for all binary views
#' @import pheatmap gridExtra ggplot2
#' @export

getGoodnessOfFit <- function(object, views="all", factors="all", plotit=T) {
  
  # Define views
  if (paste0(views,sep="",collapse="") =="all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)

  # Sort views by likelihood
  if ("ModelOpts" %in% names(getSlots("MOFAmodel"))) {
    FamilyPerView <- object@ModelOpts$likelihood
  } else {
    stop("Likelihoods not defined")
  }
  
  # Only gaussian and bernoulli views included so far
  stopifnot(all(FamilyPerView %in% c("gaussian", "bernoulli")))
  
  # Define factors
  if (paste0(factors,sep="",collapse="") == "all") { 
    factors <- factorNames(object) 
    #old object are not compatible with factro names
    if(is.null(factors)) factors <- 1:ncol(getExpectations(object,"Z","E"))
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  K <- length(factors)

  # Collect relevant expectations
  SW <- getExpectations(object,"SW","E")
  Z <- getExpectations(object,"Z","E")
  # Y <- getExpectations(object,"Y","E")
  Y <- object@TrainData
  
  # List of goodnes of fit measures by liklihood family
  GoodnessFitList <- c()

  # Calculate goodness of fit in Gaussian views using the coefficient of determination
  gaussianViews <- which(FamilyPerView=="gaussian")
  if (length(gaussianViews)>0) {
    # Calculate predictions
    Ypred_m <- lapply(views[gaussianViews], function(m) Z%*%t(SW[[m]])); names(Ypred_m) <- views[gaussianViews]
    # Ypred_mk <- lapply(views[gaussianViews], function(m) sapply(seq_along(factors), function(k) Z[,k]%*%t(SW[[m]][,k]) ) ); names(Ypred_mk) <- views[gaussianViews]
    Ypred_mk <- lapply(views[gaussianViews], function(m) sapply(factors, function(k) Z[,k]%*%t(SW[[m]][,k]) ) ); names(Ypred_mk) <- views[gaussianViews]
    
    # Calculate coefficient of determination
    #   per view
    fvar_m <- sapply(views[gaussianViews], function(m) 1 - sum((Y[[m]]-Ypred_m[[m]])**2, na.rm=T) / sum(sweep(Y[[m]],2,apply(Y[[m]],2,mean,na.rm=T),"-")**2, na.rm=T))
    #   per view and feature
    fvar_md <- lapply(views[gaussianViews], function(m) 1 - colSums((Y[[m]]-Ypred_m[[m]])**2,na.rm=T) / colSums(sweep(Y[[m]],2,apply(Y[[m]],2,mean,na.rm=T),"-")**2,na.rm=T))
    #   per factor and view
    # fvar_mk <- sapply(views[gaussianViews], function(m) sapply( seq_along(factors), function(k) 1 - sum((Y[[m]]-Ypred_mk[[m]][,k])**2) / sum(sweep(Y[[m]],2,apply(Y[[m]],2,mean),"-")**2)))
    fvar_mk <- sapply(views[gaussianViews], function(m) sapply(factors, function(k) 1 - sum((Y[[m]]-Ypred_mk[[m]][,k])**2, na.rm=T) / sum(sweep(Y[[m]],2,apply(Y[[m]],2,mean,na.rm=T),"-")**2, na.rm=T) ))
    #   per factor and view and feature
    # fvar_mdk <- lapply(views[gaussianViews], function(m) sapply( seq_along(factors), function(k) 1 - colSums((Y[[m]]-Ypred_mk[[m]][,k])**2) / colSums(sweep(Y[[m]],2,apply(Y[[m]],2,mean),"-")**2)))
    fvar_mdk <- lapply(views[gaussianViews], function(m) sapply(factors, function(k) 1 - colSums((Y[[m]]-Ypred_mk[[m]][,k])**2,na.rm=T) / colSums(sweep(Y[[m]],2,apply(Y[[m]],2,mean,na.rm=T),"-")**2,na.rm=T)))
    
    # Set names
    names(fvar_m) <- views[gaussianViews]
    names(fvar_md) <- views[gaussianViews]
    names(fvar_mdk) <- views[gaussianViews]
    for(i in names(fvar_md)) names(fvar_md[[i]]) <- colnames(object@TrainData[[i]])
    for(i in names(fvar_mdk)) rownames(fvar_mdk[[i]]) <- colnames(object@TrainData[[i]])
    colnames(fvar_mk) <- views[gaussianViews]
    rownames(fvar_mk) <- factors 
  
    # Heatmap with coefficient of determination (R2) per factor and view
    hm_gaussian <- pheatmap::pheatmap(fvar_mk, main= "", silent=T,
                          color = colorRampPalette(c("white","darkblue"))(100), 
                          cluster_cols = F, cluster_rows = F)
    fvar_m_df <- data.frame(view=names(fvar_m), R2=fvar_m)
    
    # Barplot with coefficient of determination (R2) per view
    bplt_gaussian <- ggplot( fvar_m_df, aes(x=view, y=R2)) + 
      geom_bar(stat="identity", fill="deepskyblue4", width=0.7) + 
      ggtitle("") + xlab("") +
      theme_minimal() +theme(plot.margin = unit(c(1,1,1,1), "cm"))
    
    # Join the two plots
    # gg_gaussian<-gridExtra::grid.arrange(hm_gaussian$gtable, bplt_gaussian, ncol=1, heights=c(K,5))
    gg_gaussian <- gridExtra::arrangeGrob(hm_gaussian$gtable, bplt_gaussian, ncol=1, heights=c(10,5))
    
    # Store results
    GoodnessFitList$gaussian <- list(
      R2Total = fvar_m,
      R2PerFactor = fvar_mk, 
      R2PerFactorAndFeature = fvar_mdk,
      R2PerFeature = fvar_md,
      measure = "R2")
  }
  
  
  # Calculate goodness of fit in Bernoulli views using the Brier score
  bernoulliViews <- which(FamilyPerView=="bernoulli")
  if (length(bernoulliViews)>0) {
      
    # Calculate predictions
    Yprob_m <- lapply(views[bernoulliViews], function(m) 1/(1+exp(-Z%*%t(SW[[m]])))); names(Yprob_m) <- views[bernoulliViews]
    Yprob_mk <- lapply(views[bernoulliViews], function(m) sapply(seq_along(factors), function(k) 1/(1+exp(-Z[,k]%*%t(SW[[m]][,k]))))); names(Yprob_mk) <- views[bernoulliViews]
    
    # Calculate Brier score     
    #   per view
    BS_m <- sapply(views[bernoulliViews], function(m) 1/sum(!is.na(Y[[m]]))*sum((Yprob_m[[m]]-Y[[m]])^2, na.rm=T))
    #   per view and feature
    BS_md <- lapply(views[bernoulliViews], function(m) 1/colSums(!is.na(Y[[m]]))*colSums((Yprob_m[[m]]-Y[[m]])^2, na.rm=T))
    #   per view and factor
    BS_mk <- sapply(views[bernoulliViews], function(m) sapply( seq_along(factors), function(k) 1/sum(!is.na(Y[[m]]))*sum((Yprob_mk[[m]][,k]-Y[[m]])^2, na.rm=T)))
    #   per view, feature and factor
    BS_mdk <- lapply(views[bernoulliViews], function(m) sapply( seq_along(factors), function(k) 1/colSums(!is.na(Y[[m]]))*colSums((Yprob_mk[[m]][,k]-Y[[m]])^2, na.rm=T)))
    
    # Set names
    names(BS_m) <- views[bernoulliViews]
    names(BS_md) <- views[bernoulliViews]
    names(BS_mdk) <- views[bernoulliViews]
    for(i in names(BS_md)) names(BS_md[[i]])<-colnames(object@TrainData[[i]])
    for(i in names(BS_mdk)) rownames(BS_mdk[[i]])<-colnames(object@TrainData[[i]])
    colnames(BS_mk) <- views[bernoulliViews]
    rownames(BS_mk) <- factors 
    
    # Heatmap with Brier score per factor and view
    hm_bernoulli <- pheatmap::pheatmap(BS_mk, main="", silent=T,
                          color=colorRampPalette(c("orange","white"))(100), 
                          cluster_cols=F,cluster_rows=F)
    BS_m_df <- data.frame(view=names(BS_m), BS=BS_m)
    
    # Barplot with total Brier score per view
    bplt_bernoulli<- ggplot2::ggplot( BS_m_df, aes(x=view, y=BS)) + 
      geom_bar(stat="identity", fill="darkgoldenrod", width=0.7) + 
      ggtitle("") + xlab("") +
      theme_minimal() +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    
    # Join the two plots
    # gg_bernoulli <- gridExtra::grid.arrange(hm_bernoulli$gtable, bplt_bernoulli, ncol=1, heights=c(10,5))
    gg_bernoulli <- gridExtra::arrangeGrob(hm_bernoulli$gtable, bplt_bernoulli, ncol=1, heights=c(10,5))
    
    # Store results
    GoodnessFitList$bernoulli <- list(
      BSTotal=BS_m, 
      BSPerFactor=BS_mk,
      BSPerFactorAndFeature=BS_mdk,
      BSPerFeature=BS_md,
      measure="BS")
  }
  
  
  # Assemble plot
  if (plotit) {
    if (length(bernoulliViews)>0 & length(gaussianViews)>0)
      gridExtra::grid.arrange(gg_gaussian, gg_bernoulli, ncol=2, widths=c(length(gaussianViews), length(bernoulliViews)))
    else if (length(bernoulliViews)>0)
      gridExtra::grid.arrange(gg_bernoulli, ncol=1, widths=c(length(bernoulliViews)))
    else if (length(gaussianViews)>0)
      gridExtra::grid.arrange(gg_gaussian, ncol=1, widths=c(length(gaussianViews)))
  }
  
  return(GoodnessFitList)
 
}

