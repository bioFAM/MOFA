
################################
## Functions to do subsetting ##
################################

#' @title Subset factors
#' @name subsetFactors
#' @description Method to subset (or sort) factors. \cr
#' Some factors might not be interesting for the downstream analysis and the user can choose to remove them.
#' This has no effect on the values of the other factors. For example, this could be done it the model contains factors
#' which are inactive in all views.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors character vector with the factor names (LF1,LF2,...), or numeric vector with the index of the factors.
#' @param keep_intercept logical indicating whether the optional intercept factor should be kept when subsetting (default TRUE).
#' @examples
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' MOFA_CLL_small <- subsetFactors(MOFA_CLL, factors=c(1,2,3))
#' MOFA_CLL_small <- subsetFactors(MOFA_CLL, factors=c("LF1","LF2","LF3"))
#' @export
subsetFactors <- function(object, factors, keep_intercept=TRUE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors) <= object@Dimensions[["K"]])

    # Get factors
   if(is.numeric(factors)) {
      if (object@ModelOptions$learnIntercept == T) factors <- factorNames(object)[factors+1]
      else factors <- factorNames(object)[factors]
    }
      else{ stopifnot(all(factors %in% factorNames(object))) }

  if (keep_intercept & object@ModelOptions$learnIntercept & !"intercept" %in% factors) {
    factors <- c("intercept", factors)
  }
  
  # Subset relevant slots
  object@Expectations$Z <- object@Expectations$Z[,factors, drop=F]
  object@Expectations$Alpha <- sapply(object@Expectations$Alpha, function(x) x[factors], simplify = F, USE.NAMES = T)
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
#' @description Method to subset (or sort) samples. \cr
#' This function can remove samples from the model. For example, you might want to observe the effect of Factor 1
#' on a subset of samples. You can create a new \code{\link{MOFAmodel}} excluding some samples
#' and then visualise the effect of Factor 1 on the remaining ones, for instance via 
#' \code{\link{plotDataHeatmap}} or \code{\link{plotFactorScatter}}. \cr
#' This functionality is only for exploratory purposes. 
#' In the case of outliers, we strongly recommend removing them before training the model.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param samples character vector with the sample names, numeric vector with the sample indices or 
#' logical vector with the samples to be kept as TRUE.
#' @export
#' @examples
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' # Subset samples via character vector
#' MOFA_CLL_small <- subsetSamples(MOFA_CLL, samples=c("H045","H109","H024","H056"))
#' # Subset samples via numeric vector
#' MOFA_CLL_small <- subsetSamples(MOFA_CLL, samples=1:10)
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
  if (length(object@InputData)>0)
    object@InputData <- object@InputData[,samples,]
  if (length(object@ImputedData)==0)
    object@ImputedData <- sapply(object@ImputedData, function(x) x[,samples], simplify = F, USE.NAMES = T)

  # Modify dimensionality
  object@Dimensions[["N"]] <- length(samples)
  
  # Modify sample names in the MOFAobject
  sampleNames(object) <- samples
  
  return(object)
}


#' @title Subset views
#' @name subsetViews
#' @description Method to subset (or sort) views.
#' This function can remove entire views from the model. 
#' For example, you might want to generate the \code{\link{plotVarianceExplained}} plot excluding a particular view. \cr
#' This functionality is only for exploratory purposes. If some view(s) are not of interest we strongly recommend 
#' removing them before training the model.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view names, numeric vector with the view indices or logical vector with the view to be kept as TRUE.
#' @export
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' # Subset views via character vector
#' MOFA_CLL_small <- subsetViews(MOFA_CLL, views=c("Drugs","Methylation"))
#' # Subset views via numeric vector
#' MOFA_CLL_small <- subsetViews(MOFA_CLL, views=2:3)
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