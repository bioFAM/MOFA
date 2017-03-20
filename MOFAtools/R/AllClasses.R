
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a MultiOmics Factor Analysis (MOFA) model
#' @description
#' S4 class used by MOFAtools to hold the results of a MultiView Factor Analysis trained model. \cr
#' The model has to be trained using the corresponding python package. \cr\cr
#' MOFAmodel contains the following information: \cr
#' - Training statistics: evidence lower bound, number of active factors, likelihood, etc. \cr
#' - Training options: maximum number of iterations, lower bound frequency calculation, etc. \cr
#' - Expectations: expected values of the P and Q distributions. \cr
#' - Parameters: parameters of the P and Q distributions.
#' - Dimensions: dimensionalities.
#' @section Slots:
#'  \itemize{
#'    \item{\code{Expectations}:}{Nested list with numerical matrices of all expectations from all variables. The list has three levels, accessible by name subsetting: the first level contains the variables (SW, Alpha, etc.), the second level contains the expectation name (E,E2, etc.) and third level contains the view }
#'    \item{\code{Parameters}:}{Nested list with numerical matrices of all parameters from all variables. The list has three levels, accessible by name subsetting: the first level contains the variables (SW, Alpha, etc.), the second level contains the parameter name (mean,variance, etc.) and the third level contains the view }
#'    \item{\code{TrainStats}:}{List with the training statistics.}
#'    \item{\code{TrainData}:}{List with the data views used for training.}
#'    \item{\code{TrainOpts}:}{List with the training options.}
#'    \item{\code{Dimensions}:}{List with dimensionalities of the model. N for the number of samples, D for the number of features of each view, M for the number of views and K for the number of infered latent variables.}
#'}
#' @name MOFAmodel
#' @rdname MOFAmodel
#' @aliases MOFAmodel-class
#' @exportClass MOFAmodel
setClass("MOFAmodel", 
         slots=c(Expectations="list", 
                 Parameters="list", 
                 TrainStats="list",
                 TrainOpts="list",
                 TrainData="list",
                 Dimensions="list")
)

setMethod("show", "MOFAmodel", function(object) {
  if(!.hasSlot(object,"Dimensions"))
    stop("Error: Dimensions not defined")
  cat(sprintf("MOFA trained model with the following characteristics: \n views: %d \n number of features: %s \n Number of samples: %d \n Number of latent variables: %d ",
              object@Dimensions[["M"]], paste(as.character(object@Dimensions[["D"]]),collapse=" "), object@Dimensions[["N"]], object@Dimensions[["K"]]))
})
