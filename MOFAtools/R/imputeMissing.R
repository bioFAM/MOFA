#' @title Impute missing value from a fitted MOFA model
#' @name imputeMissing
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param viewnms character vector containing the names of views to be imputed (default: all)
#' @details asd
#' @return List of imputed data, each list element corresponding to specified views.
#' @reference fill this
#' @export
#' 

#Making no use of posterior but simply taking expected values of W and Z
# Better: Full posterior of p(unobserved|observed)= integral dN(unobs; ZW, tau-1)*q(Z)*q(W)*q(alpha) d(Z,W,alpha) -->intractable?
imputeMissing <- function(model, viewnms="all"){
  if(viewnms=="all"){
    viewnms = names(model@TrainData)
  }
  #THIS SHOULD BE THE NON_CENTERED DATA IN TrainData object to get values on original scale!
  stopifnot(all(viewnms %in% names(model@TrainData)))
  
  Z<-model@Expectations$Z$E
  W<-model@Expectations$SW
  imputedData<-lapply(viewnms, function(view){
    imputedView<- Z%*% t(W[[view]]$E) 
  })    
  names(imputedData)<-viewnms
  
  #We could add a slot in the GFA model with 'imputed data'
  return(imputedData)
  }
