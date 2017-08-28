
#######################################################
## Functions to perform imputation of missing values ##
#######################################################

#' @title Impute missing values from a fitted MOFA model
#' @name imputeMissing
#' @description This function uses the latent factors and the weights infered from MOFA to impute missing values in the input views.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views vector containing the names of views (character) or index of views (numeric) to be imputed (default: "all")
#' @param factors vector with the factors indices (numeric) or factor names (character) to use (default is "all")
#' @param type type of imputations returned. By default values are imputed using "inRange". "response" gives mean for gaussian and poisson and probabilities for bernoulli , 
#' "link" gives the linear predictions, "inRange" rounds the fitted values from "terms" for integer-valued distributions to the next integer.
#' @param onlyMissing By default, only values missing in Training Data are replaced by imputed ones. If all predicitons based on MOFA are wanted, this needs to be set to FALSE.
#' @details asd
#' @return List of imputed data, each list element corresponding to specified views.
#' @references fill this
#' @export
imputeMissing <- function(object, views="all", factors = "all", type = c("inRange","response", "link"), onlyMissing = T){
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Get views  
  if (views=="all") {
    views = viewNames(object)
  } else {
    stopifnot(all(views%in%viewNames(object)))
  }
  
  # Get factors
  if (factors=="all") {
    factors = factorNames(object)
  } else {
    stopifnot(all(factors%in%factorNames(object)))
    factors <- c("intercept", factors)
  } 
  
  # Get weights
  W <- getWeights(object, views=views, factors=factors)
  
  type = match.arg(type)
  
  # mask passenger factors
  object <- detectPassengers(object)
  Z <- getFactors(object)[,factors]
  Z[is.na(Z)] <- 0 # set missing values in Z to 0 to exclude from imputations
  
  # Impute data
  imputedData <- lapply(sapply(views, grep, viewNames(object)), function(viewidx){
    
    # make imputation based on linear model
    imputedView <- t(Z%*% t(W[[viewidx]])) 
    
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
  names(imputedData) <- views
  imputedData <- imputedData[viewNames(object)]

  # save in model slot  
  object@ImputedData <- imputedData
  object
}
