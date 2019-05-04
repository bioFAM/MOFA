
# Function to automatically infer the likelihoods from the data
.inferLikelihoods <- function(object) {
  likelihood <- rep(x="gaussian", times=(getDimensions(object))$M)
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
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")  
  
  # Update node names
  if ("SW" %in% names(Expectations(object))) {
    names(Expectations(object))[names(Expectations(object))=="SW"] <- "W"
    colnames(TrainStats(object)$elbo_terms)[colnames(TrainStats(object)$elbo_terms)=="SW"] <- "W"
  }
  
  if ("AlphaW" %in% names(Expectations(object))) {
    names(Expectations(object))[names(Expectations(object)) == "AlphaW"] <- "Alpha"
    colnames(TrainStats(object)$elbo_terms)[colnames(TrainStats(object)$elbo_terms)=="AlphaW"] <- "Alpha"
  }
  
  # Update expectations
  if (is.list(Expectations(object)$Z)) {
    Expectations(object)$Z <- Expectations(object)$Z$E
    for (view in viewNames(object)) {
      Expectations(object)$Alpha[[view]] <- Expectations(object)$Alpha[[view]]$E
      Expectations(object)$W[[view]] <- Expectations(object)$W[[view]]$E
      Expectations(object)$Tau[[view]] <- Expectations(object)$Tau[[view]]$E
      Expectations(object)$Theta[[view]] <- Expectations(object)$Theta[[view]]$E
      Expectations(object)$Y[[view]] <- Expectations(object)$Y[[view]]$E
    }
  }
  
  # update model dimensions
  Dimensions(object)[["K"]] <- ncol(Expectations(object)$Z)
  
  # update learnMean to learnIntercept
  if ("learnMean" %in% names(ModelOptions(object))) {
    tmp <- names(ModelOptions(object))
    tmp[tmp=="learnMean"] <- "learnIntercept"
    names(ModelOptions(object)) <- tmp
  }
  
  # Add feature-wise means to the gaussian data to restore uncentered data in TrainData
  if (length(FeatureIntercepts(object))>=1) {
    ModelOptions(object)$learnIntercept <- NULL
    # for (m in seq_along(TrainData(object))) {
    #   if (ModelOptions(object)$likelihood[m] == "gaussian") {
    #     if (max(abs(apply(TrainData(object)[[m]],1, mean, na.rm=TRUE))) > 10^(-5))
    #       print("Warning, gaussian data seems to be uncentered")
    #     TrainData(object)[[m]] <- TrainData(object)[[m]] + as.numeric(FeatureIntercepts(object)[[m]])
    #   }
    # }
  }
  
  # update intercept to new model structure and remove intercept from pseudodata
  if(!is.null(ModelOptions(object)$learnIntercept)){
    ModelOptions(object)$learnIntercept <- as.logical(ModelOptions(object)$learnIntercept)
    if(ModelOptions(object)$learnIntercept){
      nonintercept_idx <- which(!apply(Expectations(object)$Z==1,2,all))
      intercept_idx <- which(apply(Expectations(object)$Z==1,2,all))
      if(length(intercept_idx)!=1) stop("No or multiple intercepts were learn despite using learnIntercept.")
      # save intercepts in FeatureIntercepts slot
      FeatureIntercepts(object)  <- lapply(Expectations(object)$W, function(w) w[,intercept_idx])

      #remove intercept form factors and weights
      Expectations(object)$Z <- Expectations(object)$Z[,nonintercept_idx, drop=FALSE]
      Expectations(object)$Alpha <- lapply(Expectations(object)$Alpha,
                                          function(x) x[nonintercept_idx])
      Expectations(object)$W <- lapply(Expectations(object)$W,
                                      function(x) x[,nonintercept_idx, drop=FALSE])
      Expectations(object)$Theta <- lapply(Expectations(object)$Theta,
                                          function(x) x[nonintercept_idx])
      
      # sweep out intercept from pseudodata
      for(m in seq_along(Expectations(object)$Y)) 
        Expectations(object)$Y[[m]] <- sweep(Expectations(object)$Y[[m]],2, FeatureIntercepts(object)[[m]])
      
    } else FeatureIntercepts(object)  <- lapply(Dimensions(object)$D, function(d) rep(0,d))
    ModelOptions(object)$learnIntercept <- NULL
  }

  # Add DataOptions
  if (length(DataOptions(object))==0)
    DataOptions(object) <- list(scaleViews = FALSE, removeIncompleteSamples = FALSE)
  
  # Remove depreciated and detailed model and training options
  TrainOptions(object) <- list(maxiter = TrainOptions(object)$maxiter,
                              tolerance = TrainOptions(object)$tolerance,
                              DropFactorThreshold = TrainOptions(object)$DropFactorThreshold,
                              verbose = TrainOptions(object)$verbose,
                              seed = TrainOptions(object)$seed)
  
  ModelOptions(object) <- list(likelihood = ModelOptions(object)$likelihood,
                              numFactors = ModelOptions(object)$numFactors,
                              sparsity = ModelOptions(object)$sparsity)
  
  return(object)
}

# Function to find "intercept" factors
.detectInterceptFactors <- function(object, cor_threshold = 0.75) {

  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")

  # Fetch data
  data <- getTrainData(object)
  factors <- getFactors(object)

  # Correlate the factors with global means per sample
  r <- lapply(data, function(x) abs(cor(colSums(x,na.rm=T),factors, use="pairwise.complete.obs")))
  
  token <- 0
  for (i in names(r)) {
    if (any(r[[i]]>cor_threshold)) {
      token <- 1
      message(paste0("Warning: Factor ",which(r[[i]]>cor_threshold)," is strongly correlated with the total expression for each sample in ",i))
    }
  }
  if (token==1)
    message("Such (strong) factors usually appear when count-based assays are not properly normalised by library size.")
  
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
# In such a case, the model has no information on the value of sample n on factor k, 
# and the value should be masked.
.detectPassengers <- function(object, views = "all", factors = "all", r2_threshold = 0.02) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  
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
  missing <- lapply(getTrainData(object,views), function(view) {
    sampleNames(object)[apply(view, 2, function(x) all(is.na(x)))] 
    })
  names(missing) <- viewNames(object)
  for (factor in unique_factors) {
    # view <- names(which(r2[factor,]>=r2_threshold))
    view <- colnames(r2[,which(r2[factor,]>=r2_threshold),drop=FALSE])
    missing_samples <- missing[[view]]
    if (length(missing_samples)>0) {
      Z[missing_samples,factor] <- NA
    }
  }
  
  # Replace the latent matrix
  Expectations(object)$Z <- Z
  
  return(object)
  
}

flip_factor <- function(model, factor){
  Expectations(model)$Z[,factor] <- - Expectations(model)$Z[,factor]
  for(viewnm in names(Expectations(model)$W)) {
    Expectations(model)$W[[viewnm]][,factor] <- -Expectations(model)$W[[viewnm]][,factor]
  }
return(model)
}


