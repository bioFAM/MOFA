
# (Hidden) General function to set names
.setNames <- function(object, values, dimensionality, views="all") {
  nodes <- names(Expectations(object))
  if (paste0(views,collapse="") == "all") { 
    views <- names(getDimensions(object)[["D"]]) 
  } else {
    stopifnot(all(views%in% names(getDimensions(object)[["D"]]) ))
  } 
  
  # Loop over training data
  for (m in views) {
    if (nrow(TrainData(object)[[m]]) == dimensionality)
      rownames(TrainData(object)[[m]]) <- values
    if (ncol(TrainData(object)[[m]]) == dimensionality)
      colnames(TrainData(object)[[m]]) <- values
  }
  
  
  # Loop over nodes
  for (node in nodes) {
    
    # Multi-view nodes
    if (node != "Z") {

      # Loop over views
      for (m in views) {
        
        # Loop over expectations
          if (is(Expectations(object)[[node]][[m]], "matrix")) {
            if (nrow(Expectations(object)[[node]][[m]]) == dimensionality)
              rownames(Expectations(object)[[node]][[m]]) <- values
            if (ncol(Expectations(object)[[node]][[m]]) == dimensionality)
              colnames(Expectations(object)[[node]][[m]]) <- values
          } else if (is(Expectations(object)[[node]][[m]], "array")) {
            if (length(Expectations(object)[[node]][[m]]) == dimensionality)
              names(Expectations(object)[[node]][[m]]) <- values
          }
        
      }
      
    # Single-view nodes
    } else {
      
      # Loop over expectations
        if (is(Expectations(object)[[node]], "matrix")) {
          if (nrow(Expectations(object)[[node]]) == dimensionality)
            rownames(Expectations(object)[[node]]) <- values
          if (ncol(Expectations(object)[[node]]) == dimensionality)
            colnames(Expectations(object)[[node]]) <- values
        } else if (is(Expectations(object)[[node]], "array")) {
          if (length(Expectations(object)[[node]]) == dimensionality)
            names(Expectations(object)[[node]]) <- values
        }
      
    }
  }
  
  return(object)
}
