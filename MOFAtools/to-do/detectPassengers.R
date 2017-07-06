
##########################################################################
## Functions to detect passengers that have not validated their tickets ##
##########################################################################


#' @title detectPassengers: 
#' @name detectPassengers
#' @description 
#' @param model a MOFA model
#' @param views all
#' @param factors all
#' @param ... further arguments that can be passed to pheatmap
#' @details fill this
#' @return fill this
#' @export
detectPassengers <- function(model, views = "all", factors = "all") {
  
  # Sanity checks
  if (class(model) != "MOFAmodel") stop("'model' has to be an instance of MOFAmodel")
  stopifnot(all(views %in% viewNames(model)))  
  stopifnot(all(factors %in% factorNames(model)))  
  
  # Define views
  
  # Define factors
  factors <- as.character(factors)
  if (paste0(factors,collapse="")=="all") { 
    factors <- factorNames(model) 
  } else {
    stopifnot(all(factors %in% factorNames(model)))  
  }
  
  # No idea how to do this....
  # Maybe we should do it inside imputation? Algorithm:
  (1) Aim: Imputing N samples in view m
  (2) For each latent variable, we check each sample. If the prediction is a clear outlier, we remove the effect of this latent variable?