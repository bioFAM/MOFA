
#' @title Calculate proportion of variance explained by each latent variable
#' @name CalculateProportionVariance
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
#' @import ggplot2
#' @export


CalculateProportionVariance <- function(object, plot=F) {
  
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
  if (plot==T) {
    df <- data.frame(view=viewNames(object), fvar=fvar_m)
    p <- ggplot2::ggplot(df, aes(x=view,y=fvar)) +
      geom_bar(stat="identity", fill="blue") +
      ylab("Coefficient of determination") +
      ylim(c(0,1)) + 
      theme(
        plot.margin=margin(10,10,10,10),
        axis.text.x=element_text(size=rel(1.4), color='black', margin=margin(7,0,0,0)),
        axis.text.y=element_text(size=rel(1.3), color='black'),
        axis.title.y=element_text(size=rel(1.7), margin=margin(0,15,0,0)),
        axis.title.x=element_blank(),
        axis.line = element_line(colour="black", size=0.8),
        axis.ticks.y = element_line(colour="black", size=0.5),
        axis.ticks.x = element_blank(),
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
  }
  return(fvar_mk)
}
