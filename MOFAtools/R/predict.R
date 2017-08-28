
#######################################################
## Functions to perform predict views ##
#######################################################

#' @title Predict values in a view from a fitted MOFA model
#' @name predict
#' @description This function uses the latent factors and the weights infered from MOFA to predict values in the input views.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views vector containing the names of views (character) or index of views (numeric) to be predicted (default: "all")
#' @param factors vector with the factors indices (numeric) or factor names (character) to use (default is "all")
#' @param type type of prediction returned. By default values are predicted using "inRange". "response" gives mean for gaussian and poisson and probabilities for bernoulli , 
#' "link" gives the linear predictions, "inRange" rounds the fitted values from "terms" for integer-valued distributions to the next integer.
#' @details asd
#' @return List of predicted data, each list element corresponding to specified views.
#' @references fill this
#' @export

predict <- function(object, views="all", factors = "all", type = c("inRange","response", "link")){

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
  W <- getWeights(object, views="all", factors=factors)

  # Get factors
  Z <- getFactors(object)[,factors]
  Z[is.na(Z)] <- 0 # set missing values in Z to 0 to exclude from imputations
 
  # Get type of predictions wanted 
  type = match.arg(type)
 
  # Predict data based on MOFA model
  predictedData <- lapply(sapply(views, grep, viewNames(object)), function(viewidx){
    
    # calculate terms based on linear model
    predictedView <- t(Z%*% t(W[[viewidx]])) 
    
    # make predicitons based on underlying likelihood
    if(type!="link"){
    lk <- object@ModelOpts$likelihood[viewidx]
    if(lk == "gaussian") predictedView <- predictedView
      else if (lk == "bernoulli") {predictedView <- (exp(predictedView)/(1+exp(predictedView))); if(type=="inRange") predictedView <- round(predictedView)}
        else if (lk == "poisson") {predictedView <- (exp(predictedView)); if(type=="inRange") predictedView <- round(predictedView)}
          else stop("Liklihood not implemented for imputation")
    }
    predictedView
  })

  names(predictedData) <- views

  return(predictedData)
}