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

CorrplotLFvsPC<-function(modelobject, viewname4PC, noPCs=5){
  Z<-modelobject@Expectations$Z$E
  singleview<-modelobject@TrainData[[viewname4PC]]
  pc.out<-prcomp(singleview)
  
  corrmatrix<-apply(pc.out$x[,1:noPCs],2, function(pc) {
    apply(Z,2, function(lv){
      cor(pc, lv)
    })
  })
  
  corrplot::corrplot(corrmatrix, order="original", title=viewname4PC,mar = c(1, 1, 3, 1))
  return(corrmatrix)
}

CorrplotLFvsallPC<-function(modelobject, noPCs=5){
  Z<-modelobject@Expectations$Z$E
  listPCs<-lapply(viewNames(modelobject), function(viewname4PC){
  singleview<-modelobject@TrainData[[viewname4PC]]
  pc.out<-prcomp(singleview)
  tmp<-pc.out$x[,1:noPCs]
  colnames(tmp) <- paste(viewname4PC, colnames(tmp), sep="_")
  tmp
  })
  matPCs<-do.call(cbind,listPCs)
  corrmatrix<-apply(matPCs,2, function(pc) {
    apply(Z,2, function(lv){
      cor(pc, lv)
    })
  })
  
  corrplot::corrplot(corrmatrix, order="original", title="LFs vs single-view PCs",mar = c(1, 1, 3, 1))
  return(corrmatrix)
}
