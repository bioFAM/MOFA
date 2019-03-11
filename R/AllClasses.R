
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a Multi-Omics Factor Analysis (MOFA) model
#' @description
#' The \code{MOFAmodel} is an S4 class used to store all
#'  relevant data to analyse a MOFA model.

#' @slot InputData the input data before being parsed to Training Data. 
#'    Either a MultiAssayExperiment object or a list of matrices, one per view.
#' @slot TrainData the parsed data used to fit the MOFA model
#'    A list with one matrix per view.
#' @slot ImputedData the parsed data with the missing
#'  values imputed using the MOFA model. 
#'    A list with one matrix per view.
#' @slot Expectations expected values of the different
#'  variables of the model. A list of matrices, one per variable.
#'   The most relevant are "W" for weights and "Z" for factors.
#' @slot TrainStats list with training statistics such as evidence lower bound (ELBO),
#'  number of active factors, etc.
#' @slot DataOptions list with the data processing options such as
#'  whether to center or scale the data.
#' @slot TrainOptions list with the training options such as
#'  maximum number of iterations, tolerance for convergence, etc.
#' @slot ModelOptions list with the model options such as
#'  likelihoods, number of factors, etc.
#' @slot FeatureIntercepts list with the feature-wise intercepts. 
#' Only used internally.
#' @slot Dimensions list with the relevant dimensionalities of the model.
#'  N for the number of samples, M for the number of views, 
#'  D for the number of features of each view and K for the number of infered latent factors.
#' @slot Status Auxiliary variable indicating whether the model has been trained.
#' @name MOFAmodel
#' @rdname MOFAmodel
#' @aliases MOFAmodel-class
#' @exportClass MOFAmodel
setClass("MOFAmodel", slots=c(
  InputData = "MultiAssayExperiment", TrainData = "list", ImputedData = "list",
  Expectations = "list", TrainStats = "list", Dimensions = "list",
  DataOptions = "list", TrainOptions = "list", ModelOptions = "list", FeatureIntercepts = "list",
  Status = "character")
)

setValidity("MOFAmodel", function(object) {
    if(!Status(object) %in% c("trained", "untrained")){
        return("Status(object) needs to be trained or untrained")
    }
    if(Status(object) == "trained"){
        
    if(length(Expectations(object)) == 0)
        return("Status(object) = trained but no expectations present")
        
    if(!identical(sort(c("W","Z","Theta","Tau","Alpha","Y")), sort(names(Expectations(object)))))
        return("Expectation names need to be W, Z, Theta, Tau, Alpha,Y.")
  
    if( !is.matrix(Expectations(object)[["Z"]]) |
        !is.list(Expectations(object)[["W"]]) |
        !(all(vapply(Expectations(object)[["W"]], is.matrix, logical(1)))) |
        !(is.list(Expectations(object)[["Y"]])) |
        ! (all(vapply(Expectations(object)[["Y"]], is.matrix, logical(1)))) |
        ! (is.list(Expectations(object)[["Tau"]])) |
        ! (all(vapply(Expectations(object)[["Tau"]], is.numeric, logical(1)))) |
        ! (is.list(Expectations(object)[["Alpha"]])) |
        ! (all(vapply(Expectations(object)[["Alpha"]], is.numeric, logical(1))))
    ) return("Expectations need to be a list of matrices Z and lists W, Y, Tau and Alpha of matrices")
        
    if(!(length(Expectations(object)[["Alpha"]]) == getDimensions(object)[["M"]])|
        !(length(Expectations(object)[["W"]]) == getDimensions(object)[["M"]]) |
        !(length(Expectations(object)[["Y"]]) == getDimensions(object)[["M"]]) |
        !(length(Expectations(object)[["Theta"]]) == getDimensions(object)[["M"]]) |
        !(sapply(Expectations(object)[["W"]], dim) == rbind(getDimensions(object)[["D"]],
                                                                getDimensions(object)[["K"]])) |
        !(sapply(Expectations(object)[["Y"]], dim) == rbind(getDimensions(object)[["N"]],
                                                                getDimensions(object)[["D"]])) |
        !(sapply(Expectations(object)[["Alpha"]], length) == getDimensions(object)[["K"]]) |
        !(sapply(Expectations(object)[["Theta"]], length) == getDimensions(object)[["K"]]) |
        !(ncol(Expectations(object)[["Z"]]) == getDimensions(object)[["K"]]) |
        !(nrow(Expectations(object)[["Z"]]) == getDimensions(object)[["N"]]))
            return("Dimensions of Expectations do not match model dimensions")

    }
        TRUE
})

# Printing method
setMethod("show", "MOFAmodel", function(object) {
  
  if(!.hasSlot(object,"Dimensions") | length(getDimensions(object)) == 0)
    stop("Error: Dimensions not defined")
  if(!.hasSlot(object,"Status") | length(Status(object)) == 0)
    stop("Error: Status not defined")
  
  if (Status(object) == "trained") {
    # check whether the intercept was learnt (depreciated, included for compatibility with old models)
    if(is.null(ModelOptions(object)[["learnIntercept"]])) {
      learnIntercept <- FALSE
      } else {
        learnIntercept <- ModelOptions(object)[["learnIntercept"]]
      }
    dims <- getDimensions(object)
    nfactors <- dims[["K"]]
    if (learnIntercept) { nfactors <- nfactors-1 }
    cat(sprintf("Trained MOFA model with the following characteristics:
  Number of views: %d \n View names: %s 
  Number of features per view: %s 
  Number of samples: %d 
  Number of factors: %d ",
                dims[["M"]], paste(viewNames(object),collapse=" "),
                paste(as.character(dims[["D"]]),collapse=" "),
                dims[["N"]], nfactors))
  } else {
    dims <- getDimensions(object)
    cat(sprintf("Untrained MOFA model with the following characteristics: 
  Number of views: %d 
  View names: %s 
  Number of features per view: %s
  Number of samples: %d ",
                dims[["M"]], paste(viewNames(object),collapse=" "),
                paste(as.character(dims[["D"]]),collapse=" "),
                dims[["N"]]))
    cat("\n")
  }
})
