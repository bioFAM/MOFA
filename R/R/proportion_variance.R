

CalculateResidualVariance <- function(model) {
  
  # Calculate total variance of each view
  M <- length(Y)
  var <- lapply(model$Y, function(y) apply(y,2,var) )
  residual_var <- sapply(1:M, function(m) sum( total_var[m]-1/model$tau[[m]] ))
  
  # Extract alpha matrix
  alphaMat <- cbind(data.frame(Z=1:ncol(model$alpha),alpha=t(model$alpha))
                    
  pvar <- apply(alphaMat[,-1], 2, function(x) ((1/x)*apply(model$Z,2,var))/sum((1/x)*apply(model$Z,2,var)))
}

#     var = s.var(Y,axis=0)
#     residual_var = (var - 1/tau).sum()
#     for k in xrange(self.dim["K"]):
#         factor_var = (self.dim["D"][m]/self.nodes["alpha"].Q[m].E[k]) * s.var(self.nodes["Z"].Q.E[:,k])
#         factor_pvar[m,k] = factor_var / residual_var
