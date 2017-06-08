
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a MultiOmics Factor Analysis (MOFA) model
#' @description
#' The \code{MOFAmodel} is an S4 class used by MOFAtools to store data from mult-omic experiments, provide input files for the MOFA python code and hold the results of a MOFA trained model. \cr
#' The model has to be trained using the corresponding python package. \cr\cr
#' MOFAmodel contains the following information: \cr
#' - Training statistics: evidence lower bound, number of active factors, likelihood, etc. \cr
#' - Training statistics: evidence lower bound, number of active factors, likelihood, etc. \cr
#' - Training options: maximum number of iterations, lower bound frequency calculation, etc. \cr
#' - Expectations: expected values of the P and Q distributions. \cr
#' - Parameters: parameters of the P and Q distributions.
#' - Dimensions: dimensionalities.
#' @section Slots:
#'  \itemize{
#'    \item{\code{InputData}:}{MultiAssayExperiment with all assays to be integrated in the MOFA model}
#'    \item{\code{Expectations}:}{Nested list with numerical matrices of all expectations from all variables. The list has three levels, accessible by name subsetting: the first level contains the variables (SW, Alpha, etc.), the second level contains the expectation name (E,E2, etc.) and third level contains the view }
#'    \item{\code{Parameters}:}{Nested list with numerical matrices of all parameters from all variables. The list has three levels, accessible by name subsetting: the first level contains the variables (SW, Alpha, etc.), the second level contains the parameter name (mean,variance, etc.) and the third level contains the view }
#'    \item{\code{TrainStats}:}{List with the training statistics.}
#'    \item{\code{TrainData}:}{List with the data matrices used for training.}
#'    \item{\code{ImputedData}:}{List with the data matrices containg imputed data.}
#'    \item{\code{TrainOpts}:}{List with the training options.}
#'    \item{\code{ModelOpts}:}{List with the model options.}
#'    \item{\code{Dimensions}:}{List with dimensionalities of the model. N for the number of samples, D for the number of features of each view, M for the number of views and K for the number of infered latent variables.}
#'}
#' @name MOFAmodel
#' @rdname MOFAmodel
#' @aliases MOFAmodel-class
#' @import MultiAssayExperiment
#' @exportClass MOFAmodel
setClass("MOFAmodel", 
         slots=c(InputData="MultiAssayExperiment",
                 Expectations="list", 
                 Parameters="list", 
                 TrainStats="list",
                 TrainOpts="list",
                 ModelOpts="list",
                 TrainData="list",
                 ImputedData="list",
                 Dimensions="list",
                 Status = "character")
)

setMethod("show", "MOFAmodel", function(object) {
  if(!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
    stop("Error: Dimensions not defined")
  if(!.hasSlot(object,"Status") | length(object@Status) == 0)
    stop("Error: Status not defined")
  if(object@Status == "trained") { tmp <- "Trained" } else { tmp <- "Untrained" } 
  cat(sprintf("%s MOFA model with the following characteristics: \n views: %d \n viewnames: %s \n number of features: %s \n Number of samples: %d \n Number of latent variables: %d ",
              tmp, object@Dimensions[["M"]], paste(viewNames(object),collapse=" "), paste(as.character(object@Dimensions[["D"]]),collapse=" "), object@Dimensions[["N"]], object@Dimensions[["K"]]))
})
