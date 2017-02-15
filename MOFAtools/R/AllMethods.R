
###################################
## Set and retrieve sample names ##
###################################

#' @rdname sampleNames
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases sampleNames,MOFAmodel-method
#' @author Ricard Argelaguet
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
                   if (methods::.hasSlot(object,"Dimensions"))
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

#' @rdname featureNames
#' @param object a \code{\link{MOFAmodel}} object.
#' @aliases featureNames,MOFAmodel-method
#' @author Ricard Argelaguet
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
            if (methods::.hasSlot(object,"Dimensions"))
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

#' @rdname viewNames
#' @param object a \code{\link{MOFAmodel}} object.
#' @return character vector with the names for each view
#' @author Ricard Argelaguet
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
            names(object@TrainData) <- value
          })


####################################
## Set and retrieve training data ##
####################################

#' @rdname TrainData
#' @param object a \code{\link{MOFAmodel}} object.
#' @return list of numeric matrices that contain the training data
#' @author Ricard Argelaguet
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
#' @author Ricard Argelaguet
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
#' @author Ricard Argelaguet
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
#' @author Ricard Argelaguet
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
#' @author Ricard Argelaguet
#' @rdname Parameters
#' @export
setMethod("Parameters", "MOFAmodel", function(object) { object@Parameters } )
setMethod(".Parameters<-", signature(object="MOFAmodel", value="list"),
          function(object,value) {
            object@Parameters <- value
            object
          })