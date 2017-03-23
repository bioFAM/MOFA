
#' @title Calculate proportion of variance explained by each latent variable
#' @name CalculateVariance_Views
#' @description Method to calculate the proportion of variance explained by each latent variable. \cr
#' There are several ways to do this:
#' Bayesian approach: we use the variance estimates from the variational distributions of the variables.
#'   var((E[W_k]*E[S_k])*E[Z_k])/var(Y)
#' Non-bayesian approach 1: to measure the proportion of variance explained by factor k we predict the data using only latent variable k
#' and we compute divide the predicted variance by the total variance (across samples and features).
#'    V_expl_k = var(Ypred_k)/var(Y)
#' Non-bayesian approach 2: one could argue that the non-bayesian approach 1 is not reliable because factors could have a lot of meaningless variance associated with it.
#' Therefore, better compute a goodness of fit such as the coefficient of determination:
#'   sum((Y[[m]]-Ypred_m[[m]])**2) / sum((Y[[m]]-mean(unlist(Y[[m]])))**2)
#' @param object a \code{\link{MOFAmodel}} object.
#' @param plot generate a barplot
#' @details fill this
#' @return a matrix with dim (M,K) with the proportion of residual variance explained by each factor. Barplot with the proportion of variance explained by the model for each view
#' @import ggplot2 tidyr
#' @export

CalculateVariance_Views <- function(object, views="all", factors="all", method=NULL) {
  
  # TO-DO: VIEWS, FACTORS, METHOD
  
  # Collect dimensionalities
  N <- object@Dimensions[["N"]]
  M <- object@Dimensions[["M"]]
  K <- object@Dimensions[["K"]]
  
  # Collect relevant expectations
  W <- getExpectations(object,"SW","EW")
  S <- getExpectations(object,"SW","ES")
  SW <- getExpectations(object,"SW","ESW")
  Z <- getExpectations(object,"Z","E")
  Y <- getExpectations(object,"Y","E")
  
  # Calculate observed variance
  obs_var <- sapply(names(Y), function(m) apply(Y[[m]],2,var,na.rm=T))
  
  ## Non-Bayesian approach 1 ##
  # V_expl_k = var(Ypred_k)/var(Y) = var((E[W_k]*E[S_k])*E[Z_k])/var(Y)
  # According to Damien, this is not correct  because the variance might be just bullshit signal
  # fvar_m <- sapply(names(Y), function(m) sum(apply(Y[[m]],2,var) - 1/object@Expectations$Tau[[m]]$E) / sum(apply(Y[[m]],2,var)))
  
  
  ## Non-Bayesian approach 2: coefficient of determination ##
  # This is better because it is dependent on how well our object fits the data
  Ypred_m <- lapply(1:M, function(m) Z%*%t(SW[[m]]))
  fvar_m <- sapply(1:M, function(m) 1 - sum((Y[[m]]-Ypred_m[[m]])**2) / sum((Y[[m]]-mean(unlist(Y[[m]])))**2))
  
  Ypred_mk <- lapply(1:M, function(m) lapply(1:K, function(k) Z[,k]%*%t(SW[[m]][,k]) ) )
  fvar_mk <- sapply(1:M, function(m) sapply(1:K, function(k) 1 - sum((Y[[m]]-Ypred_mk[[m]][[k]])**2) / sum((Y[[m]] - mean(unlist(Y[[m]])))**2) ) )
  
  # Set matrix names
  colnames(fvar_mk) <- viewNames(object)
  rownames(fvar_mk) <- colnames(Z)
  
  # Bar plot with the residual variance for each view
  tmp <- as.data.frame(fvar_mk)
  tmp$K <- as.factor(1:K)
  tmp <- tidyr::gather(tmp, key="view", value="fvar", -K)
  
  p <- ggplot2::ggplot(tmp, ggplot2::aes(x=view,y=fvar,fill=K)) +
    geom_bar(stat="identity", position="stack") +
    ylab("Variance explained") + xlab("") +
    # ylim(c(0,0.5)) + 
    coord_flip() +
    theme(
      # plot.margin=margin(10,10,10,10),
      axis.text.x=element_text(size=rel(1.4), color='black', margin=margin(7,0,0,0)),
      # axis.text.y=element_text(size=rel(1.3), color='black'),
      axis.text.y=element_text(size=rel(1.3), color='black'),
      axis.title.x=element_text(size=rel(1.3), margin=margin(10,0,0,0)),
      axis.title.y=element_blank(),
      # axis.line = element_line(colour="black", size=0.8),
      axis.ticks.x = element_line(colour="black", size=0.5),
      axis.ticks.y = element_blank(),
      legend.position='right',
      legend.title=element_text(size=rel(1.1)),
      legend.text=element_text(size=rel(1.0)),
      legend.key=element_rect(fill='transparent'),
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  return(list(plot=p, values=fvar_mk))
}





CalculateVarianceFeatures <- function(object, view, features, factors="all") {
  
  features <- c("ENSG00000003147","ENSG00000005189","ENSG00000005812")
  view <- "expr_mRNA"
  factors <- "all"
  
  if (factors=="all") {
    factors <- factorNames(object)  
  }
  
  # Saniy checks
  stopifnot(all(features) %in% featureNames(object)[[view]])
  stopifnot(view %in% viewNames(object))
  stopifnot(length(view)==1)
  stopifnot(length(features)<=20)
  stopifnot(all(factors %in% factorNames(object)))
  
  # Collect data and expectations
  Y <- object@TrainData[[view]]
  W <- getExpectations(object,"SW","E")[[view]]; rownames(W) <- colnames(Y)
  Z <- getExpectations(object,"Z","E"); rownames(Z) <- rownames(Y)
  
  
  # Calcualte fraction of variance explained by each factor in each feature
  fvar <- sapply(1:K, function(k) apply(as.matrix(Z[,k]) %*% t(as.matrix(W[features,k])), 2, var) / apply(Y[,features],2,var))
  
  # Prepare data for plotting
  tmp <- as.data.frame(fvar)
  tmp$gene <- rownames(tmp)
  tmp <- tidyr::gather(tmp, key="K", value="fvar", -gene)
  
  p <- ggplot2::ggplot(tmp, ggplot2::aes(x=gene,y=fvar,fill=K)) +
    geom_bar(stat="identity", position="stack") +
    ylab("Variance explained") + xlab("") +
    # ylim(c(0,0.5)) + 
    coord_flip() +
    theme(
      # plot.margin=margin(10,10,10,10),
      axis.text.x=element_text(size=rel(1.4), color='black', margin=margin(7,0,0,0)),
      # axis.text.y=element_text(size=rel(1.3), color='black'),
      axis.text.y=element_text(size=rel(1.3), color='black'),
      axis.title.x=element_text(size=rel(1.3), margin=margin(10,0,0,0)),
      axis.title.y=element_blank(),
      axis.line = element_line(colour="black", size=0.8),
      axis.ticks.x = element_line(colour="black", size=0.5),
      axis.ticks.y = element_blank(),
      legend.position='right',
      legend.title=element_text(size=rel(1.1)),
      legend.text=element_text(size=rel(1.0)),
      legend.key=element_rect(fill='transparent'),
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  print(p)
  
  return(list(plot=p, values=fvar))
}
