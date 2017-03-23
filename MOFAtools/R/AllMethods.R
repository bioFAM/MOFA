
# General function to set names to the different expectations and parameters
.setNames <- function(object, values, dimensionality) {
  stopifnot(names(object@Parameters) == names(object@Expectations))
  nodes <- names(object@Parameters)
  views <- names(object@Dimensions$D)
  
  # Loop over training data
  for (m in views) {
    if (nrow(object@TrainData[[m]]) == dimensionality)
      rownames(object@TrainData[[m]]) <- values
    if (ncol(object@TrainData[[m]]) == dimensionality)
      colnames(object@TrainData[[m]]) <- values
  }
  
  
  # Loop over nodes
  for (node in nodes) {
    
    # Multi-view nodes
    if (setequal(names(object@Expectations[[node]]),views)) {
      
      # Loop over views
      for (m in views) {
        
        # Loop over expectations
        for (expectation in names(object@Expectations[[node]][[m]])) {
          if (class(object@Expectations[[node]][[m]][[expectation]]) == "matrix") {
            if (nrow(object@Expectations[[node]][[m]][[expectation]]) == dimensionality)
              rownames(object@Expectations[[node]][[m]][[expectation]]) <- values
            if (ncol(object@Expectations[[node]][[m]][[expectation]]) == dimensionality)
              colnames(object@Expectations[[node]][[m]][[expectation]]) <- values
          } else if (class(object@Expectations[[node]][[m]][[expectation]]) == "array") {
            if (length(object@Expectations[[node]][[m]][[expectation]]) == dimensionality)
              names(object@Expectations[[node]][[m]][[expectation]]) <- values
          }
        }
        # Loop over parameters
        for (parameter in names(object@Parameters[[node]][[m]])) {
          if (class(object@Parameters[[node]][[m]][[parameter]]) == "matrix") {
            if (nrow(object@Parameters[[node]][[m]][[parameter]]) == dimensionality)
              rownames(object@Parameters[[node]][[m]][[parameter]]) <- values
            if (ncol(object@Parameters[[node]][[m]][[parameter]]) == dimensionality)
              colnames(object@Parameters[[node]][[m]][[parameter]]) <- values
          } else if (class(object@Parameters[[node]][[m]][[parameter]]) == "array") {
            if (length(object@Parameters[[node]][[m]][[parameter]]) == dimensionality)
              names(object@Parameters[[node]][[m]][[parameter]]) <- values
          }
        }
        
      }
      
    # Single-view nodes
    } else {
      
      # Loop over expectations
      for (expectation in names(object@Expectations[[node]])) {
        if (class(object@Expectations[[node]][[expectation]]) == "matrix") {
          if (nrow(object@Expectations[[node]][[expectation]]) == dimensionality)
            rownames(object@Expectations[[node]][[expectation]]) <- values
          if (ncol(object@Expectations[[node]][[expectation]]) == dimensionality)
            colnames(object@Expectations[[node]][[expectation]]) <- values
        } else if (class(object@Expectations[[node]][[expectation]]) == "array") {
          if (length(object@Expectations[[node]][[expectation]]) == dimensionality)
            names(object@Expectations[[node]][[expectation]]) <- values
        }
      }
      
      # Loop over parameters
      for (parameter in names(object@Parameters[[node]])) {
        if (class(object@Parameters[[node]][[parameter]]) == "matrix") {
          if (nrow(object@Parameters[[node]][[parameter]]) == dimensionality)
            rownames(object@Parameters[[node]][[parameter]]) <- values
          if (ncol(object@Parameters[[node]][[parameter]]) == dimensionality)
            colnames(object@Parameters[[node]][[parameter]]) <- values
        } else if (class(object@Parameters[[node]][[parameter]]) == "array") {
          if (length(object@Parameters[[node]][[parameter]]) == dimensionality)
            names(object@Parameters[[node]][[parameter]]) <- values
        }
      }
    }
  }
  
  return(object)
}

###################################
## Set and retrieve factor names ##
###################################

#' @rdname factorNames
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases factorNames,MOFAmodel-method
#' @return character vector with the features names
#' @export
setMethod("factorNames", signature(object="MOFAmodel"), function(object) { colnames(object@Expectations$Z$E) } )

#' @rdname factorNames
#' @param value a character vector of factor names
#' @import methods
#' @export
setReplaceMethod("factorNames", signature(object="MOFAmodel", value="vector"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"Expectations"))
     stop("Before assigning factor names you have to assign expectations")
   if (!methods::.hasSlot(object,"Parameters"))
     stop("Before assigning factor names you have to assign parameters")
   if (methods::.hasSlot(object,"Dimensions"))
     if (!length(value)==object@Dimensions["K"])
       stop("Length of factor names does not match the dimensionality of the latent variable matrix")
   if(!length(value)==ncol(object@Expectations$Z$E)) 
     stop("factor names do not match the number of columns in the latent variable matrix")
    
   object <- .setNames(object, value, object@Dimensions[["K"]])
   object
})



###################################
## Set and retrieve sample names ##
###################################

#' @rdname sampleNames
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases sampleNames,MOFAmodel-method
#' @return character vector with the sample names
#' @export
setMethod("sampleNames", signature(object="MOFAmodel"), function(object) { rownames(object@TrainData[[1]]) } )

#' @rdname sampleNames
#' @param value a character vector of sample names
#' @import methods
#' @export
setReplaceMethod("sampleNames", signature(object="MOFAmodel", value="vector"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"TrainData"))
     stop("Before assigning sample names you have to assign the training data")
   if (!methods::.hasSlot(object,"Expectations"))
     stop("Before assigning sample names you have to assign the expectations")
   if (!methods::.hasSlot(object,"Parameters"))
     stop("Before assigning sample names you have to assign the parameters")
   if (methods::.hasSlot(object,"Dimensions"))
     if (!length(value)==object@Dimensions["N"])
       stop("Length of sample names does not match the dimensionality of the model")
   if(!length(value)==nrow(object@TrainData[[1]])) 
     stop("sample names do not match the dimensionality of the data (rows) ")
   
    object <- .setNames(object, value, object@Dimensions[["N"]])
    object
})

####################################
## Set and retrieve feature names ##
####################################

#' @rdname featureNames
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases featureNames,MOFAmodel-method
#' @return list of character vectors with the feature names for each view
#' @export
setMethod("featureNames", signature(object="MOFAmodel"), function(object) { sapply(object@TrainData,colnames) } )

#' @rdname featureNames
#' @param value list of character vectors with the feature names for each view
#' @import methods
#' @export
setReplaceMethod("featureNames", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
    if (!methods::.hasSlot(object,"TrainData"))
      stop("Before assigning feature names you have to assign the training data")
    if (!methods::.hasSlot(object,"Expectations"))
      stop("Before assigning feature names you have to assign the expectations")
    if (!methods::.hasSlot(object,"Parameters"))
      stop("Before assigning feature names you have to assign the parameters")
    if (methods::.hasSlot(object,"Dimensions"))
      if (!all(sapply(value,length) == object@Dimensions[["D"]]))
        stop("Length of feature names does not match the dimensionality of the model")
    if (!sapply(value,length)==sapply(object@TrainData,ncol))
      stop("feature names do not match the dimensionality of the data (columns)")
    
    for (m in 1:length(object@TrainData))
      object <- .setNames(object, value[[m]], object@Dimensions[["D"]][[m]])
    object
})

#################################
## Set and retrieve view names ##
#################################

#' @rdname viewNames
#' @param object a \code{\link{MOFAmodel}} object.
#' @return character vector with the names for each view
#' @rdname viewNames
#' @export
setMethod("viewNames", signature(object="MOFAmodel"), function(object) { names(object@TrainData) } )


#' @rdname viewNames
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("viewNames<-", signature(object="MOFAmodel", value="character"), 
  function(object,value) {
    if (!methods::.hasSlot(object,"TrainData"))
      stop("Before assigning view names you have to assign the training data")
    if (methods::.hasSlot(object,"Dimensions"))
      if (!length(value) == object@Dimensions["M"])
        stop("Length of view names does not match the dimensionality of the model")
    if (!length(value)==length(object@TrainData))
      stop("view names do not match the number of views in the training data")
    
    # We have to modify this
    multiview_nodes <- c("Alpha","SW","Tau","Theta","Y")
    names(object@TrainData) <- value
    for (node in multiview_nodes) { 
      names(object@Expectations[[node]]) <- value 
      names(object@Parameters[[node]]) <- value 
    }
    
    names(object@Dimensions$D) <- value
    
    return(object)
})


####################################
## Set and retrieve training data ##
####################################

#' @rdname TrainData
#' @param object a \code{\link{MOFAmodel}} object.
#' @return list of numeric matrices that contain the training data
#' @rdname TrainData
#' @export
setMethod("TrainData", signature(object="MOFAmodel"), function(object) { object@TrainData } )

#' @import methods
setMethod(".TrainData<-", signature(object="MOFAmodel", value="list"),
  function(object,value) {
    # N <- unique(sapply(value,nrow))
    # if (length(N) > 1) 
    #   stop("Views do not have the same number of samples (rows)")
    # if (methods::.hasSlot(object,"Dimensions")) {
    #   if (object@Dimensions["M"] != length(value))
    #     if (object@Dimensions["N"] != N)
    #       stop("Number of samples in the data do not match the specified dimensionality of the model")
    #   if (all(object@Dimensions["D"] != sapply(value,ncol)))
    #     stop("Number of features in the data do not match the specified dimensionality of the model")
    # }
    object@TrainData <- value
    object
})

#######################################
## Set and retrieve training options ##
#######################################

#' @rdname TrainOpts
#' @param object a \code{\link{MOFAmodel}} object.
#' @rdname TrainOpts
#' @return list of training options
#' @export
setMethod("TrainOpts", "MOFAmodel", function(object) { object@TrainOpts } )
setMethod(".TrainOpts<-", signature(object="MOFAmodel", value="list"),
          function(object,value) {
            object@TrainOpts <- value
            object
          })

##########################################
## Set and retrieve training statistics ##
##########################################

#' @rdname TrainStats
#' @param object a \code{\link{MOFAmodel}} object.
#' @return list of training statistics
#' @export
setMethod("TrainStats", "MOFAmodel", function(object) { object@TrainStats } )
setMethod(".TrainStats<-", signature(object="MOFAmodel", value="list"),
  function(object,value) {
    object@TrainStats <- value
    object
})

###################################
## Set and retrieve expectations ##
###################################

#' @rdname Expectations
#' @param object a \code{\link{MOFAmodel}} object.
#' @rdname Expectations
#' @return list of expectations
#' @export
setMethod("Expectations", "MOFAmodel", function(object) { object@Expectations } )
setMethod(".Expectations<-", signature(object="MOFAmodel", value="list"),
  function(object,value) {
    object@Expectations <- value
    object
})

#################################
## Set and retrieve parameters ##
#################################

#' @param object a \code{\link{MOFAmodel}} object.
#' @rdname Parameters
#' @return value list of Parameters
#' @rdname Parameters
#' @export
setMethod("Parameters", "MOFAmodel", function(object) { object@Parameters } )
setMethod(".Parameters<-", signature(object="MOFAmodel", value="list"),
  function(object,value) {
    object@Parameters <- value
    object
})