###################################
## Set and retrieve expectations ##
###################################

#' @rdname Expectations
#' @description Function to set and retrieve expectations of model components.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases Expectations,MOFAmodel-method
#' @return list of matrices containing expectations of model components
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' names(Expectations(MOFAobject))


setMethod("Expectations", signature(object="MOFAmodel"),
          function(object) {
              object@Expectations
              })

#' @rdname Expectations
#' @param value a list with matrices for expectations of unobserved model components
#' @import methods
#' @export


setReplaceMethod("Expectations", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@Expectations <- value
   validObject(object)
   object
})

###################################
## Set and retrieve training options ##
###################################

#' @rdname TrainOptions
#' @description Function to set and retrieve training options.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases TrainOptions,MOFAmodel-method
#' @return list of training options
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' TrainOptions(MOFAobject)


setMethod("TrainOptions", signature(object="MOFAmodel"),
          function(object) {
              object@TrainOptions
              })

#' @rdname TrainOptions
#' @param value a list with training options
#' @import methods
#' @export


setReplaceMethod("TrainOptions", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@TrainOptions <- value
   validObject(object)
   object
})


###################################
## Set and retrieve model options ##
###################################

#' @rdname ModelOptions
#' @description Function to set and retrieve model options.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases ModelOptions,MOFAmodel-method
#' @return list of model options
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' ModelOptions(MOFAobject)


setMethod("ModelOptions", signature(object="MOFAmodel"),
          function(object) {
              object@ModelOptions
              })

#' @rdname ModelOptions
#' @param value a list with model options
#' @import methods
#' @export


setReplaceMethod("ModelOptions", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@ModelOptions <- value
   validObject(object)
   object
})



###################################
## Set and retrieve data options ##
###################################

#' @rdname DataOptions
#' @description Function to set and retrieve data options.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases DataOptions,MOFAmodel-method
#' @return list of data options
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' DataOptions(MOFAobject)


setMethod("DataOptions", signature(object="MOFAmodel"),
          function(object) {
              object@DataOptions
              })

#' @rdname DataOptions
#' @param value a list with data options
#' @import methods
#' @export


setReplaceMethod("DataOptions", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@DataOptions <- value
   validObject(object)
   object
})

###################################
## Set and retrieve model dimensions ##
###################################

#' @rdname Dimensions
#' @description Function to set and retrieve model dimensions.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases Dimensions,MOFAmodel-method
#' @return list of dimensions of the model
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' Dimensions(MOFAobject)


setMethod("Dimensions", signature(object="MOFAmodel"),
          function(object) {
              object@Dimensions
              })

#' @rdname Dimensions
#' @param value list of dimensions
#' @import methods
#' @export


setReplaceMethod("Dimensions", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@Dimensions <- value
   validObject(object)
   object
})


###################################
## Set and retrieve input data ##
###################################

#' @rdname InputData
#' @description Function to set and retrieve input data.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases InputData,MOFAmodel-method
#' @return input data as a multi-assay experiment
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' InputData(MOFAobject)


setMethod("InputData", signature(object="MOFAmodel"),
          function(object) {
              object@InputData
              })

#' @rdname InputData
#' @param value a multi-assay experiment
#' @import methods
#' @export


setReplaceMethod("InputData", signature(object="MOFAmodel", value="MultiAssayExperiment"), 
  function(object,value) {
   object@InputData <- value
   validObject(object)
   object
})


###################################
## Set and retrieve taining data ##
###################################

#' @rdname TrainData
#' @description Function to set and retrieve taining data.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases TrainData,MOFAmodel-method
#' @return list of matrices containing the training data 
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' str(TrainData(MOFAobject))


setMethod("TrainData", signature(object="MOFAmodel"),
          function(object) {
              object@TrainData
              })

#' @rdname TrainData
#' @param value a list of matrices containing the training data
#' @import methods
#' @export


setReplaceMethod("TrainData", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@TrainData <- value
   validObject(object)
   object
})


###################################
## Set and retrieve imputed data ##
###################################

#' @rdname ImputedData
#' @description Function to set and retrieve imputed data.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases ImputedData,MOFAmodel-method
#' @return list of matrices containing the imputed data 
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' MOFAobject <- impute(MOFAobject)
#' str(ImputedData(MOFAobject))


setMethod("ImputedData", signature(object="MOFAmodel"),
          function(object) {
              object@ImputedData
              })

#' @rdname ImputedData
#' @param value a list of matrices containing the imputed data
#' @import methods
#' @export


setReplaceMethod("ImputedData", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@ImputedData <- value
   validObject(object)
   object
})



###################################
## Set and retrieve training status ##
###################################

#' @rdname Status
#' @description Function to set and retrieve taining status.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases Status,MOFAmodel-method
#' @return character indicating whether the object is trained or untrained
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' Status(MOFAobject)


setMethod("Status", signature(object="MOFAmodel"),
          function(object) {
              object@Status
              })

#' @rdname Status
#' @param value character indicating whether the object is trained or untrained
#' @import methods
#' @export


setReplaceMethod("Status", signature(object="MOFAmodel", value="character"), 
  function(object,value) {
   object@Status <- value
   validObject(object)
   object
})

###################################
## Set and retrieve training statistics ##
###################################

#' @rdname TrainStats
#' @description Function to set and retrieve taining statistics
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases TrainStats,MOFAmodel-method
#' @return list of different training statistics (EBLO, number of factors)
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' str(TrainStats(MOFAobject))


setMethod("TrainStats", signature(object="MOFAmodel"),
          function(object) {
              object@TrainStats
              })

#' @rdname TrainStats
#' @param value list of different training statistics (EBLO, number of factors)
#' @import methods
#' @export


setReplaceMethod("TrainStats", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@TrainStats <- value
   validObject(object)
   object
})


###################################
## Set and retrieve feature intercepts ##
###################################

#' @rdname FeatureIntercepts
#' @description Function to set and retrieve feature intercepts
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases FeatureIntercepts,MOFAmodel-method
#' @return list of feature intercepts (per view)
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' str(FeatureIntercepts(MOFAobject))


setMethod("FeatureIntercepts", signature(object="MOFAmodel"),
          function(object) {
              object@FeatureIntercepts
              })

#' @rdname FeatureIntercepts
#' @param value list of feature intercepts (per view)
#' @import methods
#' @export


setReplaceMethod("FeatureIntercepts", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
   object@FeatureIntercepts <- value
   validObject(object)
   object
})

###################################
## Set and retrieve factor names ##
###################################

#' @rdname factorNames
#' @description Function to set and retrieve factor names.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases factorNames,MOFAmodel-method
#' @return character vector with the features names
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' factorNames(MOFAobject)
#' factorNames(MOFAobject) <- paste("Factor",1:3,sep="_")
#' factorNames(MOFAobject)

setMethod("factorNames", signature(object="MOFAmodel"), function(object) { colnames(object@Expectations$Z) } )

#' @rdname factorNames
#' @param value a character vector of factor names
#' @import methods
#' @export


setReplaceMethod("factorNames", signature(object="MOFAmodel", value="vector"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"Expectations")  | length(object@Expectations) == 0)
     stop("Before assigning factor names you have to assign expectations")
   if (methods::.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
     if (length(value)!=object@Dimensions["K"])
       stop("Length of factor names does not match the dimensionality of the latent variable matrix")
   if (length(value)!=ncol(object@Expectations$Z)) 
     stop("factor names do not match the number of columns in the latent variable matrix")
    
   object <- .setNames(object, value, object@Dimensions[["K"]])
   object
})



###################################
## Set and retrieve sample names ##
###################################

#' @rdname sampleNames
#' @description Function to set and retrieve sample names.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases sampleNames,MOFAmodel-method
#' @return character vector with the sample names
#' @export
#' @examples
#' data("CLL_data", package = "MOFAdata")
#' MOFAobject  <- createMOFAobject(CLL_data)
#' head(sampleNames(MOFAobject))
setMethod("sampleNames", signature(object="MOFAmodel"), function(object) { colnames(object@TrainData[[1]]) } )

#' @rdname sampleNames
#' @param value a character vector of sample names
#' @import methods
#' @export
setReplaceMethod("sampleNames", signature(object="MOFAmodel", value="vector"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"TrainData") | length(object@TrainData) == 0)
     stop("Before assigning sample names you have to assign the training data")
   if (!methods::.hasSlot(object,"Expectations") | length(object@Expectations) == 0)
     stop("Before assigning sample names you have to assign the expectations")
   if (methods::.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
     if (!length(value)==object@Dimensions["N"])
       stop("Length of sample names does not match the dimensionality of the model")
   if(length(value)!=ncol(object@TrainData[[1]])) 
     stop("sample names do not match the dimensionality of the data (cols) ")
   
    object <- .setNames(object, value, object@Dimensions[["N"]])
    object
})

####################################
## Set and retrieve feature names ##
####################################

#' @rdname featureNames
#' @description Function to set and retrieve feature names.
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases featureNames,MOFAmodel-method
#' @return list of character vectors with the feature names for each view
#' @export
#' @examples
#' data("CLL_data", package = "MOFAdata")
#' MOFAobject  <- createMOFAobject(CLL_data)
#' featureNames(MOFAobject)$Mutations
#' head(featureNames(MOFAobject)$Drugs)

setMethod("featureNames", signature(object="MOFAmodel"), function(object) { 
  tmp <- lapply(object@TrainData,rownames)
  names(tmp) <- viewNames(object); return(tmp) } 
)

#' @rdname featureNames
#' @param value list of character vectors with the feature names for each view
#' @import methods
#' @export
setReplaceMethod("featureNames", signature(object="MOFAmodel", value="list"), 
  function(object,value) {
    if (!methods::.hasSlot(object,"TrainData")  | length(object@TrainData) == 0)
      stop("Before assigning feature names you have to assign the training data")
    if (!methods::.hasSlot(object,"Expectations")  | length(object@Expectations) == 0)
      stop("Before assigning feature names you have to assign the expectations")
    if (methods::.hasSlot(object,"Dimensions")  | length(object@Dimensions) == 0)
      if (!all(vapply(value, length, numeric(1)) == object@Dimensions[["D"]]))
        stop("Length of feature names does not match the dimensionality of the model")
    if (!all(vapply(value, length, numeric(1)) == vapply(object@TrainData, nrow, numeric(1))))
      stop("feature names do not match the dimensionality of the data (columns)")
    
    for (m in seq_along(object@TrainData)) {
      object <- .setNames(object, value[[m]], object@Dimensions[["D"]][m], names(object@Dimensions[["D"]][m]))
    }
    return(object)
})

#################################
## Set and retrieve view names ##
#################################

#' @rdname viewNames
#' @description Function to set and retrieve view names.
#' @param object a \code{\link{MOFAmodel}} object.
#' @return character vector with the names for each view
#' @rdname viewNames
#' @export
#' @examples
#' data("CLL_data", package = "MOFAdata")
#' MOFAobject  <- createMOFAobject(CLL_data)
#' viewNames(MOFAobject)
#' viewNames(MOFAobject) <- c("DrugResponses", viewNames(MOFAobject)[2:4])
#' viewNames(MOFAobject) 

setMethod("viewNames", signature(object="MOFAmodel"), function(object) { names(object@TrainData) } )


#' @rdname viewNames
#' @param value character vector with the names for each view
#' @import methods
#' @export
setReplaceMethod("viewNames", signature(object="MOFAmodel", value="character"), 
  function(object,value) {
    if (!methods::.hasSlot(object,"TrainData") | length(object@TrainData) == 0)
      stop("Before assigning view names you have to assign the training data")
    if (methods::.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
      if (!length(value) == object@Dimensions["M"])
        stop("Length of view names does not match the dimensionality of the model")
    if (length(value)!=length(object@TrainData))
      stop("view names do not match the number of views in the training data")
    
    if (object@Status == "trained"){
      multiview_nodes <- c("Alpha","W","Tau","Theta","Y")
      for (node in multiview_nodes) { 
        names(object@Expectations[[node]]) <- value 
      }
      if (methods::.hasSlot(object,"FeatureIntercepts") & length(object@FeatureIntercepts)>0) {
        names(object@FeatureIntercepts) <- value
      }
    }
    
    names(object@TrainData) <- value
    names(object@Dimensions$D) <- value
    
    return(object)
})