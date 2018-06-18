
######################################
## Functions to perform predictions ##
######################################

#' @title Do predictions using a fitted MOFA model
#' @name predict
#' @description This function uses the factors and the corresponding weights to do data predictions.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param factors character vector with the factor name(s) or numeric vector with the factor index(es). 
#' Default is "all".
#' @param type type of prediction returned, either: 
#' \itemize{
#'  \item{\strong{response}:}{ gives the response vector, the mean for Gaussian and Poisson, and success probabilities for Bernoulli.}
#'  \item{\strong{link}:}{ gives the linear predictions.}
#'  \item{\strong{inRange}:}{ rounds the fitted values of integer-valued distributions (Poisson and Bernoulli) to the next integer.
#'  This is the default option.}
#' }
#' @param include_intercept logical indicating whether to include the optional intercept factor for the prediction.
#' Default is TRUE.
#' @details the denoised and condensed low-dimensional representation of the data captures the main sources of heterogeneity of the data. 
#' These representation can be used to do predictions using the equation Y = WX. 
#' This is the key step underlying imputation, see \code{\link{imputeMissing}} and the Methods section of the article.
#' @return Returns a list with data predictions.
#' @export
#' @examples 
#' # Example on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' # predict drug response data based on all factors
#' predict(MOFA_CLL, view="Drugs")
#' # predict all views based on all factors
#' predict(MOFA_CLL)
#' # predict mutation data based on all factors returning Bernoulli probabilities
#' predict(MOFA_CLL, view="Mutations", type="response")
#' # predict mutation data based on all factors returning binary classes
#' predict(MOFA_CLL, view="Mutations", type="inRange")
#' 
#' # Example on the scMT data
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFA_scMT <- loadModel(filepath)
#' # predict all views based on all factors (default)
#' predict(MOFA_scMT)
#' 
predict <- function(object, views = "all", factors = "all", 
                    type = c("inRange","response", "link"), 
                    include_intercept = TRUE) {

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
      if (object@ModelOptions$learnIntercept == T) factors <- factorNames(object)[factors+1]
      else factors <- factorNames(object)[factors]
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }

  # add intercept factor for prediction
  if(!"intercept" %in% factors & object@ModelOptions$learnIntercept & include_intercept) factors <- c("intercept", factors)  
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
    lks <- object@ModelOptions$likelihood
    names(lks) <- viewNames(object)
    if (type!="link") {
      lk <- lks[i]
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