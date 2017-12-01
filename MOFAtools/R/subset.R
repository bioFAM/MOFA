
################################
## Functions to do subsetting ##
################################

#' @title Subset factors
#' @name subsetFactors
#' @description Method to subset (or sort) factors
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors.
#' @export

subsetFactors <- function(object, factors) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors) <= object@Dimensions[["K"]])
  if (is.character(factors)) stopifnot(all(factors %in% factorNames(object)))
  if (is.numeric(factors)) stopifnot(all(factors %in% 1:object@Dimensions[["K"]]))
  
  # Subset expectations
  object@Expectations$Z$E <- object@Expectations$Z$E[,factors]
  object@Expectations$AlphaW <- sapply(object@Expectations$AlphaW, function(m) sapply(m, function(x) x[factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  object@Expectations$SW <- sapply(object@Expectations$SW, function(m) sapply(m, function(x) x[,factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  object@Expectations$Theta <- sapply(object@Expectations$Theta, function(m) sapply(m, function(x) x[factors], simplify = F, USE.NAMES = T), simplify = F, USE.NAMES = T)
  
  # Subset factor names
  factorNames(object) <- factors
  
  # Modify dimensionality
  object@Dimensions[["K"]] <- length(factors)
  
  return(object)
}

