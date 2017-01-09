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
  Zvar <- N*model@Parameters$Z$var
  Y <- lapply(model@Expectations$Y, function(y) y$E)
  
  # Collect dimensionalities
  M <- length(Y)
  K <- ncol(Z)
  N <- nrow(Z)
  D <- sapply(Y, ncol)
  
  obs_var <- sapply(names(Y), function(m) apply(Y[[m]],2,var,na.rm=T))
  
  ## Bayesian approach ##
  # tau <- model@Expectations$Tau
  # alpha <- model@Expectations$Alpha
  # residual_var <- sapply(names(Y), function(m) sum(apply(Y[[m]],2,var,na.rm=T) - 1/tau[[m]]$E))
  # # f_residual_var <- sapply(names(Y), function(m) sum(apply(Y[[m]],2,var) - 1/model@Expectations$Tau[[m]]$E) / sum(obs_var[[m]]))
  # 
  # m=1
  # k=8
  # 
  # var_Wdk_0 <- (1-S[[m]])*model@Parameters$SW[[m]]$var_S0
  # var_Wdk_1 <- S[[m]]*model@Parameters$SW[[m]]$var_S1
  # var_Sdk <- S[[m]]*(1-S[[m]])
  # 
  # var_SW_dk <- var_Wdk_0+var_Wdk_1+var_Sdk
  # 
  # prvar_mk <- (sum(var_SW_dk[,k])+sum(Zvar[,k])) / residual_var[m]
  # prvar_mk <- (sum(var_SW_dk[,k])+sum(Zvar[,k])) / sum(obs_var[[m]])
  # 
  # term1 <- sum(var_SW_dk[,k])*sum(Zvar[,k])
  # term2 <- sum(var_SW_dk[,k])*sum(Z[,k]**2)
  # term3 <- sum(Z[,k])*sum(var_SW_dk[,k]**2)
  # prvar_mk <- (term1 + term3 + term3) / sum(obs_var[[m]])
  # 
  ## Non-Bayesian approach 1 ##
  # V_expl_k = var(Ypred_k)/var(Y) = var((E[W_k]*E[S_k])*E[Z_k])/var(Y)
  # According to Damien, this is not correct  because the variance might be just bullshit signal
  

  ## Non-Bayesian approach 2: coefficient of determination ##
  Ypred_m <- lapply(1:M, function(m) Z%*%t(SW[[m]]))
  Ypred_mk <- lapply(1:M, function(m) lapply(1:K, function(k) Z[,k]%*%t(SW[[m]][,k]) ) )
  prvar_mk <- sapply(1:M, function(m) sapply(1:K, function(k) 1 - sum((Y[[m]]-Ypred_mk[[m]][[k]])**2) / sum((Y[[m]] - mean(unlist(Y[[m]])))**2) ) )

  # TO-DO: USE ONLY ACTIVE GENES?
  
  # Set matrix names
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
