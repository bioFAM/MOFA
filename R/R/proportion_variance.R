#################################################
## Calculate variance explained by each factor ##
#################################################

# Option 1: Bayesian approach
# Option 2: non-Bayesian approach
#   var_expl_k = var(Ypred_k)/var(Y) = var((E[W_k]*E[S_k])*E[Z_k])/var(Y)

# Considerations
#  for the sparse model we must only consider the active genes

# library(ggplot2)
# setwd("/Users/ricard/git/scGFA/R/R")
# source("loadModel.R")
# model = loadModel("/tmp/test/asd.hd5")

# Input: 
# Output:
# - Barplot with the proportion of variance explained by the model for each view
# - Returns a matrix with dim (M,K) with the proportion of residual variance explained by each factor
# CalculateProportionResidualVariance <- function(model, active_genes_threshold=0, min_genes=5, plot=F) {
CalculateProportionResidualVariance <- function(model, plot=F) {
  
  # Collect relevant expectations
  W <- lapply(model@Expectations$SW, function(w) w$EW)
  S <- lapply(model@Expectations$SW, function(w) w$ES)
  SW <- lapply(model@Expectations$SW, function(w) w$ESW)
  Z <- model@Expectations$Z$E
  Y <- lapply(model@Expectations$Y, function(y) y$E)
  M <- length(Y)
  K <- ncol(Z)
  
  ## Bayesian approach ##
  tau <- model@Expectations$Tau
  alpha <- model@Expectations$Alpha
  residual_var <- sapply(names(Y), function(m) sum(apply(Y[[m]],2,var,na.rm=T) - 1/tau[[m]]$E))
  # f_residual_var <- sapply(names(Y), function(m) sum(apply(Y[[m]],2,var) - 1/model@Expectations$Tau[[m]]$E) / sum(apply(Y[[m]],2,var)))
  D <- sapply(Y, ncol)
  prvar_mk <- sapply(names(Y), function(m) ((D[m]/alpha[[m]]$E)*apply(Z,2,var)/residual_var[m]))
  
  ## Non-Bayesian approach ##
  
  # # Calculate observed variance of each view
  # obs_var <- lapply(Y, function(y) apply(y,2,var,na.rm=T) )
  # 
  # # Calculate predicions using the model and the corresponding variance
  # Ypred_m <- lapply(1:M, function(m) Z%*%t(SW[[m]]))
  # pred_var <- lapply(Ypred_m, function(y) apply(y,2,var,na.rm=T) )
  # 
  # # Calculate proportion of residual variance (prvar) of each view (m)
  # prvar_m <- sapply(1:M, function(m) sum(pred_var[[m]])/sum(obs_var[[m]]) )
  # 
  # # Calculate fraction of residual variance explained by each factor in each view
  # # DOESNT SUM UP TO ONE
  # Ypred_mk <- lapply(1:M, function(m) lapply(1:K, function(k) Z[,k]%*%t(SW[[m]][,k]) ) )
  # # prvar_mk <- sapply(1:M, function(m) sapply(1:K, function(k) sum(apply(Ypred_mk[[m]][[k]],2,var)) / sum(obs_var[[m]]) ) )
  # prvar_mk <- sapply(1:M, function(m) sapply(1:K, function(k) sum(apply(Ypred_mk[[m]][[k]],2,var)) / sum(pred_var[[m]]) ) )
  # 
  
  # DOESNT SUM UP TO ONE IF I ONLY USE ACTIVE GENES
  # rvar_mk <- sapply(1:M, function(m) sapply(1:K, function(k) {
  #     active_genes <- S[[m]][,k] > active_genes_threshold
  #     if (sum(active_genes) < min_genes) {
  #       return(0)
  #     } else {
  #       sum(apply(Ypred_mk[[m]][[k]][,active_genes],2,var)) / sum(apply(Ypred_m[[m]][,active_genes],2,var)) 
  #     }}))
  # rvar_mk <- sapply(1:M, function(m) sapply(1:K, function(k) 
    # sum(apply(Ypred_mk[[m]][[k]][,active_genes],2,var)) / sum(apply(Ypred_m[[m]][,active_genes],2,var)) ))
  
  # maybe use alpha...
  colnames(prvar_mk) <- names(Y)
  rownames(prvar_mk) <- colnames(Z)
    
  # Bar plot with the residual variance for each view
  if (plot==T) {
    # TO-DO: PUT THE ZERO DOWN
    df <- data.frame(view=names(model@Data), fvar=rvar_m)
    ggplot(df, aes(x=view,y=fvar)) +
      geom_bar(stat="identity", fill="blue") +
      ylab("Fraction of variance explained by the model") +
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
  }
  return(prvar_mk)
}
