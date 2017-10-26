
#' @title Calculate variance explained by the model
#' @name calculateVarianceExplained
#' @description Method to calculate variance explained by the MOFA model for each view and latent factor, and optionally also for each feature.
#' As a measure of variance explained we adopt the coefficient of determination (R2).
#' For non-gaussian views the calculations are based on the normally-distributed pseudo-data (for more information on the non-gaussian model see Seeger & Bouchard, 2012).
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views Views to use, default is "all"
#' @param factors Latent factores to use, default is "all"
#' @param perFeature boolean, whether to calculate in addition the variance explained (R2) per feature (default FALSE)
#' @param perView boolean, whether to calculate in addition the variance explained (R2) per view, using all factors (default TRUE)
#' @param totalVar calculate variance explained (R2) with respect to the total variance (TRUE) or the residual variance (FALSE? 
#' @param plotit boolean, wether to produce a plot (default True)
#' @details fill this
#' @return a list with matrices with the amount of variation explained per factor and view, and optionally total variance explained per view and variance explained by each feature alone
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export

calculateVarianceExplained <- function(object, views = "all", factors = "all", perFeature = F, perView = F
                                       totalVar = T, plotit = T) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Define views
  if (paste0(views,sep="",collapse="") =="all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)

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
  Y <- getExpectations(object,"Y","E")
  
  # Calculate predictions under the MOFA model using all or a single factor
  Z[is.na(Z)] <- 0 # replace masked values on Z by 0 (do not contribute to predicitons)
  Ypred_m <- lapply(views, function(m) Z%*%t(SW[[m]])); names(Ypred_m) <- views
  Ypred_mk <- lapply(views, function(m) {
                      ltmp <- lapply(factors, function(k) Z[,k]%*%t(SW[[m]][,k]) ); names(ltmp) <- factors; ltmp
                    }); names(Ypred_mk) <- views
  
  # # Mask the predictions (unnecessary as difference below is NA already)
  # for (m in views) { 
  #   Ypred_m[[m]][which(is.na(Y[[m]]))] <- NA 
  #   for (k in factors) {
  #     Ypred_mk[[m]][[k]][which(is.na(Y[[m]]))] <- NA 
  #   }
  # }

  # Calculate prediction under the null model (intercept only)
    #by default the null model is using the intercept LF if present and not the actual mean
    NullModel <- lapply(views, function(m)  {
      # take intercept from samples which are not missing (same for all patients)
      if(object@ModelOpts$learnIntercept==T) apply(Ypred_mk[[m]][["intercept"]],2, function(c) unique(c[!is.na(c)]))
        else apply(Y[[m]],2,mean,na.rm=T)
      })
    names(NullModel) <- views
    
    resNullModel <- lapply(views, function(m) sweep(Y[[m]],2,NullModel[[m]],"-")); names(resNullModel) <- views
    partialresNull <- lapply(views, function(m) sweep(Ypred_m[[m]],2,NullModel[[m]],"-")); names(partialresNull) <- views
    
  # Remove intercept factor if present
  if(object@ModelOpts$learnIntercept==T) factorsNonconst <- factors[-1] else  factorsNonconst <- factors
    
  # Calculate coefficient of determination
    # per view
    fvar_m <- sapply(views, function(m) 1 - sum((Y[[m]]-Ypred_m[[m]])**2, na.rm=T) / sum(resNullModel[[m]]**2, na.rm=T))
     
    # per view and feature
    if (perFeature)
      fvar_md <- lapply(views, function(m) 1 - colSums((Y[[m]]-Ypred_m[[m]])**2,na.rm=T) / colSums(resNullModel[[m]]**2,na.rm=T))
    
    # per factor and view
     if (totalVar) {
       fvar_mk <- sapply(views, function(m) sapply(factorsNonconst, function(k) 1 - sum((resNullModel[[m]]-Ypred_mk[[m]][[k]])**2, na.rm=T) / sum(resNullModel[[m]]**2, na.rm=T) ))
     } else {
       fvar_mk <- sapply(views, function(m) sapply(factorsNonconst, function(k) 1 - sum((partialresNull[[m]]-Ypred_mk[[m]][[k]])**2, na.rm=T) / sum(partialresNull[[m]]**2, na.rm=T) ))
     }
    
    # per factor and view and feature
    if (perFeature)
      fvar_mdk <- lapply(views, function(m) lapply(factorsNonconst, function(k) 1 - colSums((resNullModel[[m]]-Ypred_mk[[m]][[k]])**2,na.rm=T) / colSums(resNullModel[[m]]**2,na.rm=T)))
    
    # Set names
    names(fvar_m) <- views
    if(perFeature){
      names(fvar_md) <- views
      names(fvar_mdk) <- views
      for(i in names(fvar_md)) names(fvar_md[[i]]) <- colnames(object@TrainData[[i]])
      for(i in names(fvar_mdk)) rownames(fvar_mdk[[i]]) <- colnames(object@TrainData[[i]])
    }
    colnames(fvar_mk) <- views
    rownames(fvar_mk) <- factorsNonconst 
    
    # calculate variance explained by view 
    # TO-DO: CHECK
    if (perView) {
      fvar_mk <- sapply(views, function(m) sapply(factorsNonconst, function(l) sum(sapply(factorsNonconst, function(k) cov(Ypred_mk[[m]][[l]], Ypred_mk[[m]][[k]])))))
    }
    
    # Plot the variance explained
    if (plotit) {
      
      # Sort factors
      fvar_mk_df <- reshape2::melt(fvar_mk, varnames=c("factor","view"))
      fvar_mk_df$factor <- factor(fvar_mk_df$factor)
      if (ncol(fvar_mk)>1) {
        hc <- hclust(dist(t(fvar_mk)))
        fvar_mk_df$view <- factor(fvar_mk_df$view, levels = colnames(fvar_mk)[hc$order])
      }
      fvar_mk_df$factor <- factor(fvar_mk_df$factor, levels = factorsNonconst[factor_order])
      
      # Plot 1: grid with the variance explained per factor in each view
      hm <- ggplot(fvar_mk_df, aes(view,factor)) + 
        geom_tile(aes(fill=value), color="black") +
        guides(fill=guide_colorbar("R2")) +
        scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar") +
        ylab("Latent factor") +
        theme(
          # plot.margin = margin(5,5,5,5),
          plot.title = element_text(size=17, hjust=0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.y = element_text(size=15),
          axis.line = element_blank(),
          axis.ticks =  element_blank(),
          panel.background = element_blank()
          )
      
      hm <- hm + ggtitle("Variance explained per factor")  + 
      if (totalVar) {
        guides(fill=guide_colorbar("R2"))
      } else {
        guides(fill=guide_colorbar("Residual R2")) 
      }
        
      # Plot 2: barplot with coefficient of determination (R2) per view
      fvar_m_df <- data.frame(view=names(fvar_m), R2=fvar_m)
      if (ncol(fvar_mk)>1) {
        fvar_m_df$view <- factor(fvar_m_df$view, levels = colnames(fvar_mk)[hc$order])
      }
      
      bplt <- ggplot( fvar_m_df, aes(x=view, y=R2)) + 
        ggtitle("Total variance explained per view") +
        geom_bar(stat="identity", fill="deepskyblue4", width=0.9) +
        xlab("") + ylab("R2") +
        scale_y_continuous(expand=c(0.01,0.01)) +
        theme(
          plot.margin = unit(c(1,2.4,0,0), "cm"),
          panel.background = element_blank(),
          plot.title = element_text(size=17, hjust=0.5),
          axis.ticks.x = element_blank(),
          # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.y = element_text(size=13, color="black"),
          axis.line = element_line(size=rel(1.0), color="black")
            )
      
      # Join the two plots
      # Need to fix alignment using e.g. gtable...
      # gg_R2 <- gridExtra::arrangeGrob(hm, bplt, ncol=1, heights=c(length(factorsNonconst),7) )
      # gg_R2 <- gridExtra::arrangeGrob(bplt, hm, ncol=1, heights=c(1/2,1/2) )
      # gridExtra::grid.arrange(gg_R2)
      p <- plot_grid(bplt, hm, align="v", nrow=2, rel_heights=c(1/3,2/3))
      print(p)
      
      # Plot 3: variance explained per view (TO CHECK)
      if (perView) {
        cols  <- c(RColorBrewer::brewer.pal(9, "Set1"),RColorBrewer::brewer.pal(8, "Dark2"))
        par(mfrow=c(1,2))
        barplot(t(t(fvar_mk)/colSums(fvar_mk)), col = cols, horiz = T, main = "Variance components per view", ncol = 2)
        plot.new()
        legend("center", fill=cols, legend=factorsNonconst)
      }
     
    }
    
  
  # Store results
    R2_list <- list(
      R2Total = fvar_m,
      R2PerFactor = fvar_mk)
    if (perFeature) {
      R2_list$R2PerFactorAndFeature = fvar_mdk
      R2_list$R2PerFeature = fvar_md
    }
    if (perView) {
      R2_list$PerView <- fvar_mk
    }
  
  return(R2_list)
 
}

