
# Function to automatically infer the likelihoods from the data
.inferLikelihoods <- function(object) {
  likelihood <- rep(x="gaussian", times=object@Dimensions$M)
  names(likelihood) <- viewNames(object)
  
  for (view in viewNames(object)) {
    data <- getTrainData(object, view)[[1]]
    # if (all(data %in% c(0,1,NA))) {
    if (length(unique(data[!is.na(data)]))==2) {
      likelihood[view] <- "bernoulli"
    } else if (all(data[!is.na(data)]%%1==0 & data[!is.na(data)]>=0)) {
      likelihood[view] <- "poisson"
    }
  }
  
  return(likelihood)
}

# Function to update old models
.updateOldModel <- function(object) {
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  
  # Update node names
  if ("SW" %in% names(object@Expectations)) {
    # object@ModelOptions$schedule[object@ModelOptions$schedule == "SW"] <- "W" # schedule is depreciated from ModelOptions
    names(object@Expectations)[names(object@Expectations)=="SW"] <- "W"
    colnames(object@TrainStats$elbo_terms)[colnames(object@TrainStats$elbo_terms)=="SW"] <- "W"
  }
  
  if ("AlphaW" %in% names(object@Expectations)) {
    names(object@Expectations)[names(object@Expectations) == "AlphaW"] <- "Alpha"
    colnames(object@TrainStats$elbo_terms)[colnames(object@TrainStats$elbo_terms)=="AlphaW"] <- "Alpha"
  }
  
  # Update expectations
  if (is.list(object@Expectations$Z)) {
    object@Expectations$Z <- object@Expectations$Z$E
    for (view in viewNames(object)) {
      object@Expectations$Alpha[[view]] <- object@Expectations$Alpha[[view]]$E
      object@Expectations$W[[view]] <- object@Expectations$W[[view]]$E
      object@Expectations$Tau[[view]] <- object@Expectations$Tau[[view]]$E
      object@Expectations$Theta[[view]] <- object@Expectations$Theta[[view]]$E
      object@Expectations$Y[[view]] <- object@Expectations$Y[[view]]$E
    }
  }
  
  # update learnMean to learnIntercept
  if ("learnMean" %in% names(object@ModelOptions)) {
    tmp <- names(object@ModelOptions)
    tmp[tmp=="learnMean"] <- "learnIntercept"
    names(object@ModelOptions) <- tmp
  }
  object@ModelOptions$learnIntercept <- as.logical(object@ModelOptions$learnIntercept)
  
  
  return(object)
}

# Function to find factors that act as an intercept term for the samples,
.detectInterceptFactors <- function(object, cor_threshold = 0.75) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  
  # Fetch data
  data <- getTrainData(object)
  factors <- getFactors(object, include_intercept = F)
  
  # Correlate the factors with global means per sample
  r <- lapply(data, function(x) abs(cor(apply(x,2,mean),factors, use="complete.obs")))
  for (i in names(r)) {
    if (any(r[[i]]>cor_threshold)) {
      cat(paste0("Factor ",which(r[[i]]>cor_threshold)," is capturing an intercept effect in ",i,"\n"))
      cat("Intercept factors arise from global differences between the samples, which could be different library size, mean methylation rates, etc.")
    }
  }
}


subset_augment <- function(mat, pats) {
  pats <- unique(pats)
  mat <- t(mat)
  aug_mat <- matrix(NA, ncol=ncol(mat), nrow=length(pats))
  aug_mat <- mat[match(pats,rownames(mat)),,drop=FALSE]
  rownames(aug_mat) <- pats
  colnames(aug_mat) <- colnames(mat)
  return(t(aug_mat))
}


# Function to mask passenger samples.
# Passenger samples n occur when factor k is unique to view m, but sample n is missing view m.
# In such a case, the model has no information on the value of sample n on factor k, and the value should be masked.
.detectPassengers <- function(object, views = "all", factors = "all", r2_threshold = 0.02) {
  
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
  r2 <- calculateVarianceExplained(object, views = views, factors = factors)$R2PerFactor
  unique_factors <- names(which(rowSums(r2>=r2_threshold)==1))
  
  # Mask samples that are unique in the unique factors
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
  object@Expectations$Z <- Z
  
  return(object)
  
}

flip_factor <- function(model, factor){
  model@Expectations$Z[,factor] <- - model@Expectations$Z[,factor]
  for(viewnm in names(model@Expectations$W)) {
    model@Expectations$W[[viewnm]][,factor] <- -model@Expectations$W[[viewnm]][,factor]
  }
return(model)
}


