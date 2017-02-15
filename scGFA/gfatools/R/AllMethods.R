
###################################
## Set and retrieve sample names ##
###################################


#' @param object a \code{\link{GFATrainedModel}} object.
#' @rdname sampleNames
#' @aliases sampleNames,GFATrainedModel-method
#' @author Ricard Argelaguet
#' @return character vector with the sample names
#' @export
setMethod("sampleNames", signature(object="GFATrainedModel"), function(object) { rownames(object@TrainData[[1]]) } )

#' @param value a character vector of sample names
#' @rdname sampleNames
#' @export
setReplaceMethod("sampleNames", signature(object="GFATrainedModel", value="vector"), 
                 function(object,value) {
                   if (!.hasSlot(object,"TrainData"))
                     stop("Before assigning sample names you have to assign the training data")
                   if (.hasSlot(object,"Dimensions"))
                     if (!length(value)==object@Dimensions["N"])
                       stop("Length of sample names does not match the dimensionality of the model")
                   if(!length(value)==nrow(object@TrainData[[1]])) 
                     stop("sample names do not match the dimensionality of the data (rows) ")
                   for (m in 1:length(object@TrainData))
                     rownames(object@TrainData[[m]]) <- value
                 })

####################################
## Set and retrieve feature names ##
####################################

#' @param object a \code{\link{GFATrainedModel}} object.
#' @rdname featureNames
#' @aliases featureNames,GFATrainedModel-method
#' @author Ricard Argelaguet
#' @return list of character vectors with the feature names for each view
#' @export
setMethod("featureNames", signature(object="GFATrainedModel"), function(object) { sapply(object@TrainData,colnames) } )

#' @rdname featureNames
#' @param value list of character vectors with the feature names for each view
#' @export
setReplaceMethod("featureNames", signature(object="GFATrainedModel", value="list"), 
          function(object,value) {
            if (!.hasSlot(object,"TrainData"))
              stop("Before assigning feature names you have to assign the training data")
            if (.hasSlot(object,"Dimensions"))
              if (!all(sapply(value,length) == object@Dimensions["D"]))
                stop("Length of feature names does not match the dimensionality of the model")
            if (!sapply(value,length)==sapply(object@TrainData,ncol))
              stop("feature names do not match the dimensionality of the data (columns)")
            for (m in 1:length(object@TrainData))
              colnames(object@TrainData[[m]]) <- value[[m]]
          })

#################################
## Set and retrieve view names ##
#################################

#' @param object a \code{\link{GFATrainedModel}} object.
#' @rdname viewNames
#' @return character vector with the names for each view
#' @author Ricard Argelaguet
#' @rdname viewNames
#' @export
setMethod("viewNames", signature(object="GFATrainedModel"), function(object) { names(object@TrainData) } )


#' @param value character vector with the names for each view
#' @rdname viewNames
#' @export
setMethod("viewNames<-", signature(object="GFATrainedModel", value="list"), 
          function(object,value) {
            if (!.hasSlot(object,"TrainData"))
              stop("Before assigning view names you have to assign the training data")
            if (.hasSlot(object,"Dimensions"))
              if (!length(value) == object@Dimensions["M"])
                stop("Length of view names does not match the dimensionality of the model")
            if (!length(value)==length(object@TrainData))
              stop("view names do not match the number of views in the training data")
            names(object@TrainData) <- value
          })


####################################
## Set and retrieve training data ##
####################################

#' @param object a \code{\link{GFATrainedModel}} object.
#' @rdname TrainData
#' @return list of numeric matrices that contain the training data
#' @author Ricard Argelaguet
#' @rdname TrainData
#' @export
setMethod("TrainData", signature(object="GFATrainedModel"), function(object) { object@TrainData } )

#' @rdname TrainData
#' @param value list of numeric matrices that contain the training data
#' @export
setMethod("TrainData<-", signature(object="GFATrainedModel", value="list"),
          function(object,value) {
            N <- unique(sapply(value,nrow))
            if (length(N) > 1) 
              stop("Views do not have the same number of samples (rows)")
            if (.hasSlot(object,"Dimensions")) {
              if (object@Dimensions["M"] != length(value))
                if (object@Dimensions["N"] != N)
                  stop("Number of samples in the data do not match the specified dimensionality of the model")
              if (all(object@Dimensions["D"] != sapply(value,ncol)))
                stop("Number of features in the data do not match the specified dimensionality of the model")
            }
            object@TrainData <- value
            object
          })

#######################################
## Set and retrieve training options ##
#######################################

#' @param object a \code{\link{GFATrainedModel}} object.
#' @rdname TrainOpts
#' @return list of training options
#' @author Ricard Argelaguet
#' @rdname TrainOpts
#' @export
setMethod("TrainOpts", "GFATrainedModel", function(object) { object@TrainOpts } )

#' @rdname TrainOpts
#' @param value list of training options
#' @export
setMethod("TrainOpts<-", signature(object="GFATrainedModel", value="list"),
          function(object,value) {
            object@TrainOpts <- value
            object
          })

##########################################
## Set and retrieve training statistics ##
##########################################

#' @param object a \code{\link{GFATrainedModel}} object.
#' @rdname TrainStats
#' @return list of training statistics
#' @author Ricard Argelaguet
#' @rdname TrainStats
#' @export
setMethod("TrainStats<-", signature(object="GFATrainedModel", value="list"),
          function(object,value) {
            object@TrainStats <- value
            object
          })

#' @rdname TrainStats
#' @param value list of training statistics
#' @export
setMethod("TrainStats", "GFATrainedModel", function(object) { object@TrainStats } )

###################################
## Set and retrieve expectations ##
###################################

#' @param object an \code{\link{GFATrainedModel}} object.
#' @rdname Expectations
#' @return list of expectations
#' @author Ricard Argelaguet
#' @rdname Expectations
#' @export
setMethod("Expectations", "GFATrainedModel", function(object) { object@Expectations } )

#' @rdname Expectations
#' @param value list of expectations
#' @export
setMethod("Expectations<-", signature(object="GFATrainedModel", value="list"),
          function(object,value) {
            object@Expectations <- value
            object
          })

#################################
## Set and retrieve parameters ##
#################################

#' @param object an \code{\link{GFATrainedModel}} object.
#' @rdname Parameters
#' @return value list of Parameters
#' @author Ricard Argelaguet
#' @rdname Parameters
#' @export
setMethod("Parameters", "GFATrainedModel", function(object) { object@Parameters } )

#' @rdname Parameters
#' @param value list of Parameters
#' @export
setMethod("Parameters<-", signature(object="GFATrainedModel", value="list"),
          function(object,value) {
            object@Parameters <- value
            object
          })