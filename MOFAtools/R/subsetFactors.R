
################################
## Functions to do subsetting ##
################################

#' @title subsetFactors: subset factors
#' @name subsetFactors
#' @description Method to subset or sort factors
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors vector with the factors indices (numeric) or factor names (character) to subset.
#' @export

subsetFactors <- function(object, factors) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors) <= object@Dimensions[["K"]])
  if (is.character(factors)) stopifnot(all(factors %in% factorNames(object)))
  if (is.numeric(factors)) stopifnot(all(factors %in% 1:object@Dimensions[["K"]]))
  
  # Sort Expectations
  object@Expectations$Z$E <- object@Expectations$Z$E[,factors]
  object@Expectations$AlphaW <- sapply(object@Expectations$AlphaW, function(m) sapply(m, function(x) x[factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  object@Expectations$SW <- sapply(object@Expectations$SW, function(m) sapply(m, function(x) x[,factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  object@Expectations$Theta <- sapply(object@Expectations$Theta, function(m) sapply(m, function(x) x[factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  
  # Sort Parameters
  object@Parameters$Z <- sapply(object@Parameters$Z, function(x) x[,factors], simplify = F, USE.NAMES = T)
  object@Parameters$SW <- sapply(object@Parameters$SW, function(m) sapply(m, function(x) x[,factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  object@Parameters$Theta <- sapply(object@Parameters$Theta, function(m) sapply(m, function(x) x[factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  object@Parameters$AlphaW <- sapply(object@Parameters$AlphaW, function(m) sapply(m, function(x) x[factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  
  # Sort factor names
  factorNames(object) <- factors
  
  # Modify dimensionality
  object@Dimensions[["K"]] <- length(factors)
  
  return(object)
}


#' @title subsetSamples: subset samples
#' @name subsetSamples
#' @description Method to subset (or sort) samples
#' @param object a \code{\link{MOFAmodel}} object.
#' @param samples numeric or character vector with the sample indices (numeric) or sample names (character) to subset (or sort)
#' @export

subsetSamples <- function(object, samples) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(samples) <= object@Dimensions[["N"]])
  if (is.character(samples)) stopifnot(all(samples %in% sampleNames(object)))
  if (is.numeric(samples)) stopifnot(all(samples %in% 1:object@Dimensions[["N"]]))
  
  # Subset Expectations
  stop()
  object@Expectations$Z$E <- object@Expectations$Z$E[samples,]
  # object@Expectations$SW <- sapply(object@Expectations$SW, function(m) sapply(m, function(x) x[,factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  # object@Expectations$Theta <- sapply(object@Expectations$Theta, function(m) sapply(m, function(x) x[,factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  
  # Subset Parameters
  # object@Parameters$Z <- sapply(object@Parameters$Z, function(x) x[,factors], simplify = F, USE.NAMES = T)
  # object@Parameters$SW <- sapply(object@Parameters$SW, function(m) sapply(m, function(x) x[,factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  # object@Parameters$Theta <- sapply(object@Parameters$Theta, function(m) sapply(m, function(x) x[,factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  # object@Parameters$AlphaW <- sapply(object@Parameters$AlphaW, function(m) sapply(m, function(x) x[factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  
  # Modify sample names
  factorNames(object) <- factors
  
  # Modify dimensionality
  object@Dimensions[["K"]] <- length(factors)
  
  return(object)
}

