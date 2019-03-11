
################################
## Functions to do subsetting ##
################################

#' @title Subset factors
#' @name subsetFactors
#' @description Method to subset (or sort) factors. \cr
#' Some factors might not be interesting for the downstream analysis and the user can choose to remove them.
#' This has no effect on the values of the other factors.
#' For example, this could be done if the model contains factors
#' which are inactive in all views.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors character vector with the factor names (LF1,LF2,...),
#'  or numeric vector with the index of the factors.
#' @return \code{\link{MOFAmodel}} object with a subset of factors
#' @examples
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")
#' MOFA_CLL <- loadModel(filepath)
#' MOFA_CLL_small <- subsetFactors(MOFA_CLL, factors=c(1,2,3))
#' MOFA_CLL_small
#' MOFA_CLL_small <- subsetFactors(MOFA_CLL, factors=c("LF1","LF2","LF3"))
#' MOFA_CLL_small
#' @export
subsetFactors <- function(object, factors) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors) <= getDimensions(object)[["K"]])

    # Get factors
   if(is.numeric(factors)) {
      factors <- factorNames(object)[factors]
    }
      else{ stopifnot(all(factors %in% factorNames(object))) }

  # Subset relevant slots
  Expectations(object)$Z <- Expectations(object)$Z[,factors, drop=FALSE]
  Expectations(object)$Alpha <- lapply(Expectations(object)$Alpha,
                                      function(x) x[factors])
  Expectations(object)$W <- lapply(Expectations(object)$W,
                                  function(x) x[,factors, drop=FALSE])
  Expectations(object)$Theta <- lapply(Expectations(object)$Theta,
                                      function(x) x[factors])

  # Modify dimensionality
  Dimensions(object)[["K"]] <- length(factors)
  
  # Modify factor names
  factorNames(object) <- as.character(factors)
  
  return(object)
}



#' @title Subset samples
#' @name subsetSamples
#' @description Method to subset (or sort) samples. \cr
#' This function can remove samples from the model. For example,
#' you might want to observe the effect of Factor 1 on a subset of samples.
#' You can create a new \code{\link{MOFAmodel}} excluding some samples
#' and then visualise the effect of Factor 1 on the remaining ones, for instance via 
#' \code{\link{plotDataHeatmap}} or \code{\link{plotFactorScatter}}. \cr
#' This functionality is only for exploratory purposes. 
#' In the case of outliers, we strongly recommend removing them before training the model.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param samples character vector with the sample names, numeric vector with the sample indices or 
#' logical vector with the samples to be kept as TRUE.
#' @return \code{\link{MOFAmodel}} object with a subset of samples
#' @export
#' @examples
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")
#' MOFA_CLL <- loadModel(filepath)
#' # Subset samples via character vector
#' MOFA_CLL_small <- subsetSamples(MOFA_CLL, samples=c("H045","H109","H024","H056"))
#' MOFA_CLL_small
#' # Subset samples via numeric vector
#' MOFA_CLL_small <- subsetSamples(MOFA_CLL, samples=1:10)
#' MOFA_CLL_small
subsetSamples <- function(object, samples) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(samples) <= Dimensions(object)[["N"]])
  # warning("Warning: removing samples is fine for an exploratory analysis...\nbut we recommend removing them before training!\n")
  
  # Get samples
  if (is.character(samples)) {
    stopifnot(all(samples %in% sampleNames(object)))
  } else {
    samples <- sampleNames(object)[samples]
  }
  
  # Subset relevant slots
  Expectations(object)$Z <- Expectations(object)$Z[samples,, drop=FALSE]
  Expectations(object)$Y <- lapply(Expectations(object)$Y, function(x) x[samples,])
  TrainData(object) <- lapply(TrainData(object), function(x) x[,samples])
  if (length(InputData(object))>0)
    InputData(object) <- InputData(object)[,samples,]
  if (length(ImputedData(object))==0)
    ImputedData(object) <- lapply(ImputedData(object), function(x) x[,samples])

  # Modify dimensionality
  Dimensions(object)[["N"]] <- length(samples)
  
  # Modify sample names in the MOFAobject
  sampleNames(object) <- samples
  
  return(object)
}


#' @title Subset views
#' @name subsetViews
#' @description Method to subset (or sort) views.
#' This function can remove entire views from the model. 
#' For example, you might want to generate the \code{\link{plotVarianceExplained}} plot
#'  excluding a particular view. \cr
#' This functionality is only for exploratory purposes.
#' If some view(s) are not of interest we strongly recommend removing them before training the model.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view names, numeric vector with the view indices
#'  or logical vector with the view to be kept as TRUE.
#' @return \code{\link{MOFAmodel}} object with a subset of views
#' @export
#' @examples
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")
#' MOFA_CLL <- loadModel(filepath)
#' # Subset views via character vector
#' MOFA_CLL_small <- subsetViews(MOFA_CLL, views=c("Drugs","Methylation"))
#' MOFA_CLL_small
#' # Subset views via numeric vector
#' MOFA_CLL_small <- subsetViews(MOFA_CLL, views=2:3)
#' MOFA_CLL_small
#' 
subsetViews <- function(object, views) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(views) <= Dimensions(object)[["N"]])
  warning("Removing views is fine for an exploratory analysis,\n
          but we recommend removing them before training!\n")
  
  # Get views
  if (is.character(views)) {
    stopifnot(all(views %in% viewNames(object)))
  } else {
    views <- viewNames(object)[views]
  }
  
  # Subset relevant slots
  Expectations(object)$Y <- Expectations(object)$Y[views]
  Expectations(object)$W <- Expectations(object)$W[views]
  TrainData(object) <- TrainData(object)[views]
  if (length(ImputedData(object))==0) { ImputedData(object) <- ImputedData(object)[views] }
  
  # Modify dimensionality
  Dimensions(object)[["M"]] <- length(views)
  
  # Modify sample names in the MOFAobject
  viewNames(object) <- views
  
  return(object)
}