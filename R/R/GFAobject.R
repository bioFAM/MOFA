
#####################################################################
## Script to define a general class to store a scGFA trained model ##
#####################################################################

# GFATrainedModel contains:
# - Training stats: evidence lower bound, number of active factors, likelihood...
# - Training options: maximum number of iterations, lower bound frequency calculation...
# - Expectations: expected values of the Q distributions
# - Parameters: parameters of the Q distributions

# setClass(Class="GFATrainedModel", slots=c(Expectations="list", Parameters="list", TrainStats="list", TrainOpts="list", TrainData="list"), )
setClass(Class="GFATrainedModel", slots=c(Expectations="list", Parameters="list", TrainStats="list", TrainOpts="list", Data="list"), )

# Define the show method
setMethod("show", "GFATrainedModel", function(object) { cat("scGFAtrained model ") } )

# Define method to return training statistics
setGeneric(name="getTrainStats", def=function(object) { standardGeneric("getTrainStats") } )
setMethod("getTrainStats", "GFATrainedModel", function(object) { return(object@TrainStats) } )

# Define method to return training options
setGeneric(name="getTrainOpts", def=function(object) { standardGeneric("getTrainOpts") })
setMethod("getTrainOpts", "GFATrainedModel", function(object) { return(object@TrainOpts) } )

# Define method to return expectations
setGeneric(name="getExpectations", def=function(object) { standardGeneric("getExpectations") })
setMethod("getExpectations", "GFATrainedModel", function(object) { return(object@Expectations) } )

# Define method to return parameters
setGeneric(name="getParameters", def=function(object) { standardGeneric("getParameters") })
setMethod("getParameters", "GFATrainedModel", function(object) { return(object@Parameters) } )


# showMethods(class="GFATrainedModel")
# getMethod(f="show",signature="GFATrainedModel")

# myGFA <- new("GFATrainedModel", TrainStats=, Expectations=, Parameters=, Nodes=)
