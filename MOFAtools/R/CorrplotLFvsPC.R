library(corrplot)

#' @title Correlate latent factors to principal components of specified view
#' @name CorrplotLFvsPC
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param viewname4PC name of view to do PCA on
#' @param noPCs number of PCs to compare to (default = 5)
#' @details asd
#' @return Correlation Matrix Latent factors versus Principal components
#' @reference fill this
#' @import corrplot
#' @export

CorrplotLFvsPC<-function(model, viewname4PC, noPCs=5){
  Z<-model@Expectations$Z$E
  singleview<-model@TrainData[[viewname4PC]]
  pc.out<-prcomp(singleview)
  
  corrmatrix<-apply(pc.out$x[,1:noPCs],2, function(pc) {
    apply(Z,2, function(lv){
      cor(pc, lv)
    })
  })
  
  corrplot::corrplot(corrmatrix, order="original", title=viewname4PC,mar = c(1, 1, 3, 1))
  return(corrmatrix)
}


# methods: svd, ppca, bpca
CorrplotLFvsallPC<-function(model, noPCs=5, method="svd"){
  
  # Collect expectations
  Z <- getExpectations(model,"Z","E")
  
  # Perform PCAs
  listPCs <- lapply(viewNames(model), function(m) {
    # pc.out<-prcomp(model@TrainData[[m]])
    pc.out <- pcaMethods::pca(model@TrainData[[m]], method=method, center=TRUE, nPcs=noPCs)
    # tmp <- pc.out$x[,1:noPCs]
    tmp <- pc.out@scores
    colnames(tmp) <- paste(m, colnames(tmp), sep="_")
    tmp
  })
  
  # Calculate correlation matrix between latent factors and PCs
  matPCs <- do.call(cbind,listPCs)
  corrmatrix <- cor(matPCs,Z)
  
  corrplot::corrplot(corrmatrix, order="original", title="", mar = c(1, 1, 3, 1))
  return(corrmatrix)
}
