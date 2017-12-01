subset_augment <- function(mat, pats) {
  pats <- unique(pats)
  mat <- t(mat)
  aug_mat <- matrix(NA, ncol=ncol(mat), nrow=length(pats))
  aug_mat <- mat[match(pats,rownames(mat)),,drop=FALSE]
  rownames(aug_mat) <- pats
  colnames(aug_mat) <- colnames(mat)
  return(t(aug_mat))
}


detectPassengers <- function(object, views = "all", factors = "all", r2_threshold = 0.03) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Define views
  if (paste0(views,sep="",collapse="") =="all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)
  
  # Define factors
  factors <- as.character(factors)
  if (paste0(factors,collapse="")=="all") { 
    factors <- factorNames(object) 
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  
  # Collect relevant data
  Z <- getFactors(object)
  
  # Identify factors unique to a single view by calculating relative R2 per factor
  r2 <- calculateVarianceExplained(object, views = views, factors = factors, plotit = F, totalVar = T)$R2PerFactor
  unique_factors <- names(which(rowSums(r2>=r2_threshold)==1))
  
  # Mask samples that have full missing views
  missing <- sapply(getTrainData(object,views), function(view) sampleNames(object)[apply(view, 2, function(x) all(is.na(x)))] )
  names(missing) <- viewNames(object)
  for (factor in unique_factors) {
    # view <- names(which(r2[factor,]>=r2_threshold))
    view <- colnames(r2[,which(r2[factor,]>=r2_threshold),drop=F])
    missing_samples <- missing[[view]]
    if (length(missing_samples)>0) {
      Z[missing_samples,factor] <- NA
    }
  }
  
  # Replace the latent matrix
  object@Expectations$Z$E <- Z
  
  return(object)
  
}


