
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
#' @details asd
#' @return List of imputed data, each list element corresponding to specified views.
#' @references fill this
#' @export
imputeMissing <- function(object, views = "all", factors = "all", type = c("inRange","response", "link")){
  
  type = match.arg(type)

  # sanity checks are perfomred inside the predcit function
  predData <- predict(object, views=views, factors = factors, type = type)

  # Get views  
  if (views=="all") {
    views = viewNames(object)
  } else {
    stopifnot(all(views%in%viewNames(object)))
  }

  #replace NAs with predicted values
  imputedData <- getTrainData(object, views = views)
  imputedData <- lapply(names(imputedData), function(viewnm) {
      view <- imputedData[[viewnm]]
      non_observed <- which(is.na(view), arr.ind = T)
      if(viewnm %in% names(predData)) view[non_observed] <- predData[[viewnm]][non_observed]
      view
  })

  # re- arrange list in accordance with other data slots in the model
  names(imputedData) <- views
  imputedData <- imputedData[viewNames(object)]
  names(imputedData) <- viewNames(object)


  # save in model slot  
  object@ImputedData <- imputedData
  object
}
