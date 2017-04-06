
#' @title Calculate proportion of variance explained by each latent variable in each view
#' @name CalculateVariance_Views
#' @description Method to calculate the proportion of variance explained by each latent variable in each view. \cr
#' There are several ways to do this: \cr
#' Bayesian approach: we use the variance estimates from the variational distributions of the variables.\cr
#'   \eqn{var((E[W_k]*E[S_k])*E[Z_k])/var(Y)}
#' Non-bayesian approach 1: to measure the proportion of variance explained by factor k we predict the data using only latent variable k and we compute divide the predicted variance by the total variance (across samples and features): \cr
#'    \eqn{V_expl_k = var(Ypred_k)/var(Y)}
#' (default) Non-bayesian approach 2: one could argue that the non-bayesian approach 1 is not reliable because factors could have a lot of meaningless variance associated with it. Therefore, better compute a goodness of fit such as the coefficient of determination: \cr
#'   \eqn{sum((Y[[m]]-Ypred_m[[m]])**2) / sum((Y[[m]]-mean(unlist(Y[[m]])))**2)}
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views Views to use, default is "all"
#' @param factors Latent variables or factores to use, default is "all"
#' @param method (ignore for now)
#' @details fill this
#' @return a matrix with the proportion of variance explained by each factor (columns) in each view (rows)
#' @import ggplot2 tidyr
#' @export

CalculateVariance_Views <- function(object, views="all", factors="all", method=NULL) {
  
  # Define views
  if (views=="all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)
  
  # Define factors
  if (factors=="all") { 
    factors <- factorNames(object) 
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  K <- length(factors)
  
  # Collect relevant expectations
  # W <- getExpectations(object,"SW","EW")
  # S <- getExpectations(object,"SW","ES")
  SW <- getExpectations(object,"SW","E")
  Z <- getExpectations(object,"Z","E")
  Y <- getExpectations(object,"Y","E")
  
  # Calculate observed variance
  obs_var <- sapply(names(Y), function(m) apply(Y[[m]],2,var,na.rm=T))
  
  # Non-Bayesian approach 1
  # V_expl_k = var(Ypred_k)/var(Y) = var((E[W_k]*E[S_k])*E[Z_k])/var(Y)
  # According to Damien, this is not correct  because the variance might be just bullshit signal
  fvar_m <- sapply(names(Y), function(m) sum(apply(Y[[m]],2,var) - 1/object@Expectations$Tau[[m]]$E) / sum(apply(Y[[m]],2,var)))
  
  # Non-Bayesian approach 2: coefficient of determination
  # This is better because it is dependent on how well our object fits the data
  Ypred_m <- sapply(views, function(m) Z%*%t(SW[[m]]))
  fvar_m <- sapply(views, function(m) 1 - sum((Y[[m]]-Ypred_m[[m]])**2) / sum((Y[[m]]-mean(unlist(Y[[m]])))**2))
  
  Ypred_mk <- sapply(views, function(m) sapply(factors, function(k) Z[,k]%*%t(SW[[m]][,k]) ) )
  fvar_mk <- sapply(views, function(m) sapply(factors, function(k) 1 - sum((Y[[m]]-Ypred_mk[[m]][,k])**2) / sum((Y[[m]] - mean(unlist(Y[[m]])))**2) ) )
  
  # Set matrix names
  colnames(fvar_mk) <- views
  rownames(fvar_mk) <- factors
  
  # Bar plot with the residual variance for each view and factor
  # tmp <- as.data.frame(fvar_mk)
  # tmp$K <- as.factor(1:K)
  # tmp <- tidyr::gather(tmp, key="view", value="fvar", -K)
  # p <- ggplot2::ggplot(tmp, aes(x=view,y=fvar,fill=K)) +
  #   geom_bar(stat="identity", position="stack") +
  #   ylab("Variance explained") + xlab("") +
  #   # ylim(c(0,0.5)) + 
  #   guides(fill=guide_legend(title="Latent variable", title.position="top", title.hjust=0.5)) +
  #   coord_flip() +
  #   theme(
  #     plot.margin=margin(20,20,20,20),
  #     axis.text.x=element_text(size=rel(1.5), color='black', margin=margin(7,0,0,0)),
  #     axis.text.y=element_text(size=rel(1.5), color='black', margin=margin(0,7,0,0)),
  #     axis.title.x=element_text(size=rel(1.3), margin=margin(10,0,0,0)),
  #     axis.title.y=element_blank(),
  #     axis.line = element_line(colour="black"),
  #     axis.ticks.x = element_line(colour="black", size=0.5),
  #     axis.ticks.y = element_blank(),
  #     legend.position='top',
  #     legend.direction="horizontal",
  #     legend.title=element_text(size=rel(1.4)),
  #     legend.text=element_text(size=rel(1.2)),
  #     legend.key=element_rect(fill='transparent'),
  #     panel.border=element_blank(),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank()
  #   )
  # print(p)
  
  # Barplot with the residual variance for each view
  tmp <- data.frame(view=names(fvar_m), fvar=round(fvar_m,2))
  
  p <- ggplot(tmp, aes(x=view,y=fvar)) +
    geom_bar(stat="identity", fill="steelblue") +
    ylab("Coefficient of determination") +
    scale_y_continuous(expand=c(0,0.01), lim=c(0,1)) +
    geom_text(aes(label=fvar), vjust=1.6, color="white", position = position_dodge(0.9), size=5.5) +
    theme(
      plot.margin=margin(10,10,10,10),
      axis.text.x=element_text(size=rel(1.8), color='black', margin=margin(7,0,0,0), angle=90),
      axis.text.y=element_text(size=rel(1.3), color='black'),
      axis.title.y=element_text(size=rel(1.7), margin=margin(0,15,0,0)),
      axis.title.x=element_blank(),
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(colour="black", size=0.75),
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
  
  # return(list(plot=p, values=fvar_mk))
  return(fvar_mk)
}





# CalculateVarianceFeatures <- function(object, view, features, factors="all") {
# 
#   features <- c("ENSG00000003147","ENSG00000005189","ENSG00000005812")
#   view <- "expr_mRNA"
#   factors <- "all"
# 
#   if (factors=="all") {
#     factors <- factorNames(object)
#   }
# 
#   # Saniy checks
#   stopifnot(all(features) %in% featureNames(object)[[view]])
#   stopifnot(view %in% viewNames(object))
#   stopifnot(length(view)==1)
#   stopifnot(length(features)<=20)
#   stopifnot(all(factors %in% factorNames(object)))
# 
#   # Collect data and expectations
#   Y <- object@TrainData[[view]]
#   W <- getExpectations(object,"SW","E")[[view]]; rownames(W) <- colnames(Y)
#   Z <- getExpectations(object,"Z","E"); rownames(Z) <- rownames(Y)
# 
# 
#   # Calcualte fraction of variance explained by each factor in each feature
#   fvar <- sapply(1:K, function(k) apply(as.matrix(Z[,k]) %*% t(as.matrix(W[features,k])), 2, var) / apply(Y[,features],2,var))
# 
#   # Prepare data for plotting
#   tmp <- as.data.frame(fvar)
#   tmp$gene <- rownames(tmp)
#   tmp <- tidyr::gather(tmp, key="K", value="fvar", -gene)
# 
#   p <- ggplot2::ggplot(tmp, ggplot2::aes(x=gene,y=fvar,fill=K)) +
#     geom_bar(stat="identity", position="stack") +
#     ylab("Variance explained") + xlab("") +
#     # ylim(c(0,0.5)) +
#     coord_flip() +
#     theme(
#       # plot.margin=margin(10,10,10,10),
#       axis.text.x=element_text(size=rel(1.4), color='black', margin=margin(7,0,0,0)),
#       # axis.text.y=element_text(size=rel(1.3), color='black'),
#       axis.text.y=element_text(size=rel(1.3), color='black'),
#       axis.title.x=element_text(size=rel(1.3), margin=margin(10,0,0,0)),
#       axis.title.y=element_blank(),
#       axis.line = element_line(colour="black", size=0.8),
#       axis.ticks.x = element_line(colour="black", size=0.5),
#       axis.ticks.y = element_blank(),
#       legend.position='right',
#       legend.title=element_text(size=rel(1.1)),
#       legend.text=element_text(size=rel(1.0)),
#       legend.key=element_rect(fill='transparent'),
#       panel.border=element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank()
#     )
#   print(p)
# 
#   return(list(plot=p, values=fvar))
# }
