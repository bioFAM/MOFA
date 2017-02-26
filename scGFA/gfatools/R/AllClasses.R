
###########################################################
## Define a general class to store a scGFA trained model ##
###########################################################

#' @title central object for a trained GFA model
#' @description
#' S4 class and the main class used by gfatools to hold the results of a Group Factor Analysis trained model \cr
#' The model has to be trained using the corresponding python package. \cr
#' GFATrainedModel contains: \cr
#' - Training stats: evidence lower bound, number of active factors, likelihood... \cr
#' - Training options: maximum number of iterations, lower bound frequency calculation... \cr
#' - Expectations: expected values of the Q distributions \cr
#' - Parameters: parameters of the Q distributions
#' @section Slots:
#'  \describe{
#'    \item{\code{Expectations}:}{Nested list with numerical matrices of all expectations from all variables. Three levels: first level define variables (SW, Alpha, etc.), second level defines the expectation name (E,E2, etc.) and third layer is the view (Expression, Methylation, etc. }
#'    \item{\code{Parameters}:}{Nested list with numerical matrices of all parameters from all variables. Three levels: first level define variables (SW, Alpha, etc.), second level defines the parameter name (mean,variance, etc.) and third layer is the view (Expression, Methylation, etc. }
#'    \item{\code{TrainStats}:}{List with training statics such as lower bound, number of active latent variables, etc.}
#'    \item{\code{TrainData}:}{List with the data views used for training.}
#'    \item{\code{TrainOpts}:}{List with training options.}
#'    \item{\code{Dimensions}:}{List with dimensionalities of the model. N for the number of samples, D for the number of features of each view, M for the number of views and K for the number of infered latent variables.}
#'}
#' @name GFATrainedModel
#' @rdname GFATrainedModel
#' @aliases GFATrainedModel-class
#' @exportClass GFATrainedModel
setClass("GFATrainedModel", 
         slots=c(Expectations="list", 
                 Parameters="list", 
                 TrainStats="list",
                 TrainOpts="list",
                 TrainData="list",
                 Dimensions="list")
)

setMethod("show", "GFATrainedModel", function(object) {
  if(!.hasSlot(object,"Dimensions"))
    stop("Error: Dimensions not defined")
  cat(sprintf("scGFAtrained model with the following characteristics: \n views: %d \n number of features: %s \n Number of samples: %d \n Number of latent variables: %d ",
              object@Dimensions[["M"]], stringr::str_c(as.character(object@Dimensions[["D"]]),collapse=","), object@Dimensions[["N"]], object@Dimensions[["K"]]))
})
