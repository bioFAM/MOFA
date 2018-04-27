
################################
## Functions to do subsetting ##
################################

#' @title Subset factors
#' @name subsetFactors
#' @description Method to subset (or sort) factors
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors.
#' @param keep_intercept bool whether intercept is kept when subsetting (default TRUE).
#' @export
subsetFactors <- function(object, factors, keep_intercept=T) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors) <= object@Dimensions[["K"]])

    # Get factors
   if(is.numeric(factors)) {
      if (object@ModelOptions$learnIntercept == T) factors <- factorNames(object)[factors+1]
      else factors <- factorNames(object)[factors]
    }
      else{ stopifnot(all(factors %in% factorNames(object))) }

  if (keep_intercept & object@ModelOptions$learnIntercept == T & !"intercept" %in% factors) {
    factors <- c("intercept", factors)
  }
  
  # Subset relevant slots
  object@Expectations$Z <- object@Expectations$Z[,factors, drop=F]
  object@Expectations$AlphaW <- sapply(object@Expectations$AlphaW, function(x) x[factors], simplify = F, USE.NAMES = T)
  object@Expectations$W <- sapply(object@Expectations$W, function(x) x[,factors, drop=F], simplify = F, USE.NAMES = T)
  object@Expectations$Theta <- sapply(object@Expectations$Theta, function(x) x[factors], simplify = F, USE.NAMES = T)

  # Modify dimensionality
  object@Dimensions[["K"]] <- length(factors)
  
  # Modify factor names
  factorNames(object) <- as.character(factors)
  
  return(object)
}



#' @title Subset samples
#' @name subsetSamples
#' @description Method to subset (or sort) samples
#' @param object a \code{\link{MOFAmodel}} object.
#' @param samples character vector with the sample names, numeric vector with the sample indices or logical vector with the samples to be kept as TRUE.
#' @export
subsetSamples <- function(object, samples) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(samples) <= object@Dimensions[["N"]])
  warning("Removing samples is fine for an exploratory analysis, but we recommend removing them before training!\n")
  
  # Get samples
  if (is.character(samples)) {
    stopifnot(all(samples %in% sampleNames(object)))
  } else {
    samples <- sampleNames(object)[samples]
  }
  
  # Subset relevant slots
  object@Expectations$Z <- object@Expectations$Z[samples,, drop=F]
  object@Expectations$Y <- sapply(object@Expectations$Y, function(x) x[samples,], simplify = F, USE.NAMES = T)
  object@TrainData <- sapply(object@TrainData, function(x) x[,samples], simplify = F, USE.NAMES = T)
  object@InputData <- object@InputData[,samples,] 
  if (length(object@ImputedData)==0) { object@ImputedData <- sapply(object@ImputedData, function(x) x[,samples], simplify = F, USE.NAMES = T)}

  # Modify dimensionality
  object@Dimensions[["N"]] <- length(samples)
  
  # Modify sample names in the MOFAobject
  sampleNames(object) <- samples
  
  return(object)
}


#' @title Subset views
#' @name subsetViews
#' @description Method to subset (or sort) views
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view names, numeric vector with the view indices or logical vector with the view to be kept as TRUE.
#' @export
subsetViews <- function(object, views) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(views) <= object@Dimensions[["N"]])
  warning("Removing views is fine for an exploratory analysis, but we recommend removing them before training!\n")
  
  # Get views
  if (is.character(views)) {
    stopifnot(all(views %in% viewNames(object)))
  } else {
    views <- viewNames(object)[views]
  }
  
  # Subset relevant slots
  object@Expectations$Y <- object@Expectations$Y[views]
  object@Expectations$W <- object@Expectations$W[views]
  object@TrainData <- object@TrainData[views]
  if (length(object@ImputedData)==0) { object@ImputedData <- object@ImputedData[views] }
  
  # Modify dimensionality
  object@Dimensions[["M"]] <- length(views)
  
  # Modify sample names in the MOFAobject
  viewNames(object) <- views
  
  return(object)
}