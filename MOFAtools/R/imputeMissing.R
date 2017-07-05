
####################################################
## Functions perform imputation of missing values ##
####################################################

#' @title Impute missing value from a fitted MOFA model
#' @name imputeMissing
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param viewnms character vector containing the names of views to be imputed (default: all)
#' @param type of imputations returned. By default values are imputed using "inRange". "response" gives mean for gaussian and poisson and probabilities for bernoulli , 
#' "link" gives the linear predictions, "inRange" rounds the fitted values from "terms" for integer-valued distributions to the next integer.
#' @param onlyMissing By default, only values missing in Training Data are replaced by imputed ones. If all predicitons based on MOFA are wanted, this needs to be set to FALSE. 
#' @details asd
#' @return List of imputed data, each list element corresponding to specified views.
#' @references fill this
#' @export
imputeMissing <- function(object, viewnms="all", type = c("inRange","response", "link"), onlyMissing =T, factors = "all"){
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  if(viewnms=="all"){
    viewnms = viewNames(object)
  }
  
  if(factors=="all"){
    factors = factorNames(object)
  } else factors <- c("0", factors)
  
  type = match.arg(type)
  stopifnot(all(viewnms %in% viewNames(object)))
  
  Z <- object@Expectations$Z$E[, factors]
  W <- lapply(object@Expectations$SW, function(list) list$E[,factors])
  
  imputedData<-lapply(sapply(viewnms, grep, viewNames(object)), function(viewidx){
    
    # make imputation based on linear model
    imputedView <- t(Z%*% t(W[[viewidx]]$E)) 
    
    # make predicitons based on underlying model
    if(type!="link"){
    lk <- object@ModelOpts$likelihood[viewidx]
    if(lk == "gaussian") imputedView <- imputedView
      else if (lk == "bernoulli") {imputedView <- (exp(imputedView)/(1+exp(imputedView))); if(type=="inRange") imputedView <- round(imputedView)}
        else if (lk == "poisson") {imputedView <- (exp(imputedView)); if(type=="inRange") imputedView <- round(imputedView)}
          else stop("Liklihood not implemented for imputation")
    }
    # values that have been observed are kept
    if(onlyMissing){
      observed <- which(!is.na(object@TrainData[[viewidx]]), arr.ind = T)
      imputedView[observed] <- object@TrainData[[viewidx]][observed]
    }
    imputedView
  })    

  # re- arrange list in accordance with other data slots in the model
  names(imputedData) <- viewnms
  imputedData <- imputedData[viewNames(object)]

  # save in model slot  
  object@ImputedData <- imputedData
  object
  }
