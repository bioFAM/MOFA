
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
#'  \item{\strong{response}:}{ gives the response vector, the mean for Gaussian and Poisson,
#'   and success probabilities for Bernoulli.}
#'  \item{\strong{link}:}{ gives the linear predictions.}
#'  \item{\strong{inRange}:}{ rounds the fitted values of integer-valued distributions 
#'  (Poisson and Bernoulli) to the next integer.
#'  This is the default option.}
#' }
#' @details Matrix factorization models generate a denoised and condensed low-dimensional representation 
#' of the data which capture the main sources of heterogeneity of the data. 
#' Such representation can be used to do predictions (data reconstruction) and imputation (see \code{\link{impute}}). \cr
#' For mathematical details, see the Methods section of the MOFA article.
#' @return Returns a list with data predictions.
#' @export
#' @examples 
#' library(ggplot2)
#' 
#' # Example on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")
#' MOFA_CLL <- loadModel(filepath)
#' 
#' # predict drug response data using all factors
#' predictedDrugs <- predict(MOFA_CLL, view="Drugs")
#' 
#' # predict all views using all factors (default)
#' predictedAll <- predict(MOFA_CLL)
#' 
#' # predict Mutation data using all factors, returning Bernoulli probabilities
#' predictedMutations <- predict(MOFA_CLL, view="Mutations", type="response")
#' 
#' # predict Mutation data using all factors, returning binary classes
#' predictedMutationsBinary <- predict(MOFA_CLL, view="Mutations", type="inRange")
#' 
#' # Compare the predictions with the true data
#' pred <- as.numeric(predictedAll$Drugs)
#' true <- as.numeric(getTrainData(MOFA_CLL)$Drugs)
#' qplot(pred,true) + geom_hex(bins=100) + coord_equal() + 
#'    geom_abline(intercept=0, slope=1, col="red")
#' 
#' # Example on the scMT data
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFA_scMT <- loadModel(filepath)
#' 
#' # Predict all views using all factors (default)
#' predictedAll <- predict(MOFA_scMT)
#'  
#' # Compare the predictions with the true data
#' view <- "RNA expression"
#' pred <- as.numeric(predictedAll[[view]])
#' true <- as.numeric(getTrainData(MOFA_scMT)[[view]])
#' qplot(pred,true) + geom_hex(bins=100) + coord_equal() + 
#'    geom_abline(intercept=0, slope=1, col="red") 

predict <- function(object, views = "all", factors = "all", 
                    type = c("inRange","response", "link")) {

  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  
  # Get views  
  if (is.numeric(views)) {
    stopifnot(all(views <= getDimensions(object)[["M"]]))
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
      factors <- factorNames(object)[factors]
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }


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
    predictedView <- sweep(t(Z%*% t(W[[i]])),1, -FeatureIntercepts(object)[[i]]) # add intercept row-wise
    
    # make predicitons based on underlying likelihood
    lks <- ModelOptions(object)[["likelihood"]]
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