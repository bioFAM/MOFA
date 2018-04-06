
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
      if (object@ModelOpts$learnIntercept == T) factors <- factorNames(object)[factors+1]
      else factors <- factorNames(object)[factors]
    }
      else{ stopifnot(all(factors %in% factorNames(object))) }

  if(keep_intercept & object@ModelOpts$learnIntercept == T & !"intercept" %in% factors) factors <- c("intercept", factors)
  # Subset expectations
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

