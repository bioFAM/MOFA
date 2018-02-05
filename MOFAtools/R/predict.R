
######################################
## Functions to perform predictions ##
######################################

#' @title Do predictions using a fitted MOFA model
#' @name predict
#' @description This function uses the latent factors and the weights to do predictions in the input data
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es), default is "all".
#' @param factors character vector with the factor name(s) or numeric vector with the factor index(es), default is "all".
#' @param type type of prediction returned. "response" gives mean for gaussian and poisson, and probabilities for bernoulli , 
#' @param include_intercept logical indicating whether to include the intercept factors for the prediction (default is TRUE)
#' "link" gives the linear predictions, "inRange" rounds the fitted values from "terms" for integer-valued distributions to the next integer. Default is "inRange".
#' @details the denoised and condensed low-dimensional representation of the data captures the main sources of heterogeneity of the data. 
#' These representation can be used to do predictions using the equation Y = WX. This is the key step underlying imputation, see \code{\link{imputeMissing}} and Methods section of the article.
#' @return List with data predictions, each element corresponding to a view.
#' @export

predict <- function(object, views = "all", factors = "all", type = c("inRange","response", "link"), include_intercept = TRUE) {

  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Get views  
  if (is.numeric(views)) {
    stopifnot(all(views<=object@Dimensions$M))
    views <- viewNames(object)[views] 
  } else {
    if (paste0(views,sep="",collapse="") =="all") { 
      views = viewNames(object)
    } else {
      stopifnot(all(views%in%viewNames(object)))
    }
  }
  
  # Get factors
  if (paste0(factors,collapse="") == "all") { 
    factors <- factorNames(object) 
  } else if(is.numeric(factors)) {
      if (object@ModelOpts$learnIntercept == T) factors <- factorNames(object)[factors+1]
      else factors <- factorNames(object)[factors]
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }

  # add intercept factor for prediction
  if(!"intercept" %in% factors & object@ModelOpts$learnIntercept & include_intercept) factors <- c("intercept", factors)  
  if(!include_intercept & "intercept" %in% factors) factors <- factors[factors!="intercept"]
  
  # Get type of predictions wanted 
  type = match.arg(type)
  
  # Collect weights
  W <- getWeights(object, views=views, factors=factors)

  # Collect factors
  Z <- getFactors(object)[,factors]
  Z[is.na(Z)] <- 0 # set missing values in Z to 0 to exclude from imputations
 
  # Predict data based on MOFA model
  # predictedData <- lapply(sapply(views, grep, viewNames(object)), function(viewidx){
  predictedData <- lapply(views, function(i){
    
    # calculate terms based on linear model
    predictedView <- t(Z%*% t(W[[i]])) 
    
    # make predicitons based on underlying likelihood
    if (type!="link") {
      lk <- object@ModelOpts$likelihood[i]
      if (lk == "gaussian") { 
        predictedView <- predictedView 
      }
      else if (lk == "bernoulli") { 
        predictedView <- (exp(predictedView)/(1+exp(predictedView)))
        if (type=="inRange") predictedView <- round(predictedView)
      } else if (lk == "poisson") { 
        # predictedView <- (exp(predictedView))
        predictedView <- log(1 + exp(predictedView))
        if(type=="inRange") predictedView <- round(predictedView)
      }
      else { 
        stop(sprintf("Likelihood %s not implemented for imputation",lk)) 
      }
    }
    predictedView
  })

  names(predictedData) <- views

  return(predictedData)
}