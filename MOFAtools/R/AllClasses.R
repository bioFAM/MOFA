
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a Multi-Omics Factor Analysis (MOFA) model
#' @description
#' The \code{MOFAmodel} is an S4 class used to store all relevant data to analyse a MOFA model.
#' @section Slots:
#'  \itemize{
#'    \item{\code{InputData}:}{ the input data before being parsed to Training Data. Either a MultiAssayExperiment object or a list of matrices.}
#'    \item{\code{TrainData}:}{ the parsed data used to fit the MOFA model. A list of matrices.}
#'    \item{\code{ImputedData}:}{ the parsed data with the missing values imputed using the MOFA model. A list of matrices.}
#'    \item{\code{Expectations}:}{ expected values of the different variables of the model. Nested list of two levels, accessible by name subsetting: the first level contains the variable names (Z, SW, Alpha, etc.) and the second level contains the expectation name (E,E2, etc.) }
#'    \item{\code{TrainStats}:}{List with the training statistics, for example evidence lower bound (ELBO), final number of active factors, etc.}
#'    \item{\code{DataOpts}:}{List with the data processing options, for example whether to center or scale the data.}
#'    \item{\code{TrainOpts}:}{List with the training options, for example maximum number of iterations, tolerance for convergence, lower bound frequency calculation, etc.}
#'    \item{\code{ModelOpts}:}{List with the model options, for example data likelihoods, initial number of factors, use of sparsity, covariates, etc.}
#'    \item{\code{Dimensions}:}{List with the relevant dimensionalities of the model. N for the number of samples, D for the number of features of each view, M for the number of views and K for the number of infered latent factors.}
#'    \item{\code{Status}:}{Auxiliary variable indicating whether the model has been trained or not}
#'}
#' @name MOFAmodel
#' @rdname MOFAmodel
#' @aliases MOFAmodel-class
#' @exportClass MOFAmodel
setClass("MOFAmodel", 
         slots=c(InputData = "MultiAssayExperiment",
                 TrainData = "list",
                 ImputedData = "list",
                 Expectations = "list", 
                 TrainStats = "list",
                 TrainOpts = "list",
                 DataOpts = "list",
                 ModelOpts = "list",
                 Dimensions = "list",
                 Status = "character")
)

# Printing method
setMethod("show", "MOFAmodel", function(object) {
  
  if(!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
    stop("Error: Dimensions not defined")
  if(!.hasSlot(object,"Status") | length(object@Status) == 0)
    stop("Error: Status not defined")
  
  if (object@Status == "trained") {
    nfactors <- object@Dimensions[["K"]]
    if (object@ModelOpts$learnIntercept==TRUE) { nfactors <- nfactors-1 }
    cat(sprintf("Trained MOFA model with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of samples: %d \n Number of factors: %d ",
                object@Dimensions[["M"]], paste(viewNames(object),collapse=" "), paste(as.character(object@Dimensions[["D"]]),collapse=" "), object@Dimensions[["N"]], nfactors))
  } else {
    cat(sprintf("Untrained MOFA model with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of samples: %d ",
                object@Dimensions[["M"]], paste(viewNames(object),collapse=" "), paste(as.character(object@Dimensions[["D"]]),collapse=" "), object@Dimensions[["N"]]))
    cat("\n")
  }
})
