
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a MultiOmics Factor Analysis (MOFA) model
#' @description
#' The \code{MOFAmodel} is an S4 class used to store all relevant data to analyse a MOFA model. \cr
#' MOFAmodel contains the following information: \cr
#' - Training statistics: evidence lower bound, number of active factors, likelihood, etc. \cr
#' - Training options: maximum number of iterations, lower bound frequency calculation, etc. \cr
#' - Model options: likelihoods, use of covariates, etc. \cr
#' - Expectations: expected values of the distributions of the different variables. \cr
#' - Parameters: parameters of the distributions of the different variables. \cr
#' - Dimensions: dimensionalities. \cr
#' - Training Data: the data used for fitting the MOFA model. \cr
#' - Input Data: the input data before being parsed to Training Data.
#' @section Slots:
#'  \itemize{
#'    \item{\code{InputData}:}{Optional, MultiAssayExperiment object with all assays to be integrated in the MOFA model}
#'    \item{\code{Expectations}:}{Nested list with numerical matrices of all expectations from all variables. The list has three levels, accessible by name subsetting: the first level contains the variables (SW, Alpha, etc.), the second level contains the expectation name (E,E2, etc.) and third level contains the view. }
#'    \item{\code{Parameters}:}{Nested list with numerical matrices of all parameters from all variables. The list has three levels, accessible by name subsetting: the first level contains the variables (SW, Alpha, etc.), the second level contains the parameter name (mean,variance, etc.) and the third level contains the view. }
#'    \item{\code{TrainStats}:}{List with the training statistics.}
#'    \item{\code{TrainData}:}{List with the data matrices used for training.}
#'    \item{\code{DataOpts}:}{List with the data processing options.}
#'    \item{\code{TrainOpts}:}{List with the training options.}
#'    \item{\code{ModelOpts}:}{List with the model options.}
#'    \item{\code{Dimensions}:}{List with dimensionalities of the model. N for the number of samples, D for the number of features of each view, M for the number of views and K for the number of infered latent variables.}
#'    \item{\code{ImputedData}:}{Optional, List with the data matrices containing imputed data.}
#'    \item{\code{Status}:}{Auxiliary variable indicating whether the model has been trained or not}
#'}
#' @name MOFAmodel
#' @rdname MOFAmodel
#' @aliases MOFAmodel-class
#' @exportClass MOFAmodel
setClass("MOFAmodel", 
         slots=c(InputData="MultiAssayExperiment",
                 Expectations="list", 
                 Parameters="list", 
                 TrainStats="list",
                 DataOpts="list",
                 TrainOpts="list",
                 ModelOpts="list",
                 TrainData="list",
                 ImputedData="list",
                 Dimensions="list",
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
  }
  else {
    cat(sprintf("Untrained MOFA model with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of samples: %d ",
                object@Dimensions[["M"]], paste(viewNames(object),collapse=" "), paste(as.character(object@Dimensions[["D"]]),collapse=" "), object@Dimensions[["N"]]))
  }
})
