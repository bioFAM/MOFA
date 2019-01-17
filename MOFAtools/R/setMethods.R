
# (Hidden) General function to set names
.setNames <- function(object, values, dimensionality, views="all") {
  nodes <- names(object@Expectations)
  if (paste0(views,collapse="") == "all") { 
    views <- names(object@Dimensions$D) 
  } else {
    stopifnot(all(views%in%names(object@Dimensions$D) ))
  } 
  
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
    if (node != "Z") {
    # if (setequal(names(object@Expectations[[node]]),views)) {
      
      # Loop over views
      for (m in views) {
        
        # Loop over expectations
          if (class(object@Expectations[[node]][[m]]) == "matrix") {
            if (nrow(object@Expectations[[node]][[m]]) == dimensionality)
              rownames(object@Expectations[[node]][[m]]) <- values
            if (ncol(object@Expectations[[node]][[m]]) == dimensionality)
              colnames(object@Expectations[[node]][[m]]) <- values
          } else if (class(object@Expectations[[node]][[m]]) == "array") {
            if (length(object@Expectations[[node]][[m]]) == dimensionality)
              names(object@Expectations[[node]][[m]]) <- values
          }
        
      }
      
    # Single-view nodes
    } else {
      
      # Loop over expectations
        if (class(object@Expectations[[node]]) == "matrix") {
          if (nrow(object@Expectations[[node]]) == dimensionality)
            rownames(object@Expectations[[node]]) <- values
          if (ncol(object@Expectations[[node]]) == dimensionality)
            colnames(object@Expectations[[node]]) <- values
        } else if (class(object@Expectations[[node]]) == "array") {
          if (length(object@Expectations[[node]]) == dimensionality)
            names(object@Expectations[[node]]) <- values
        }
      
    }
  }
  
  return(object)
}

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
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' factorNames(MOFAobject)

setMethod("factorNames", signature(object="MOFAmodel"), function(object) { colnames(object@Expectations$Z) } )

#' @rdname factorNames
#' @param value a character vector of factor names
#' @import methods
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' factorNames(MOFAobject)
#' factorNames(MOFAobject) <- paste("Factor",1:3,sep="_")
#' factorNames(MOFAobject)


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
#' data("CLL_data")
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
#' data("CLL_data")
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
      if (!all(sapply(value,length) == object@Dimensions[["D"]]))
        stop("Length of feature names does not match the dimensionality of the model")
    if (!all(sapply(value,length)==sapply(object@TrainData,nrow)))
      stop("feature names do not match the dimensionality of the data (columns)")
    
    for (m in 1:length(object@TrainData)) {
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
#' data("CLL_data")
#' MOFAobject  <- createMOFAobject(CLL_data)
#' viewNames(MOFAobject)
setMethod("viewNames", signature(object="MOFAmodel"), function(object) { names(object@TrainData) } )


#' @rdname viewNames
#' @param value character vector with the names for each view
#' @import methods
#' @export
#' @examples
#' data("CLL_data")
#' MOFAobject  <- createMOFAobject(CLL_data)
#' viewNames(MOFAobject) 
#' viewNames(MOFAobject) <- c("DrugResponses", viewNames(MOFAobject)[2:4])
setMethod("viewNames<-", signature(object="MOFAmodel", value="character"), 
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
