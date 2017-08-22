##########################################################
## Functions to perform Feature Set Enrichment Analysis ##
##########################################################

#' @title Feature Set Enrichment Analysis
#' @name FeatureSetEnrichmentAnalysis 
#' @description Method to perform enrichment analysis for feature sets from a given view in the latent factors  \cr
#' Here we use a slightly modified version of the \link{PCGSE} package.
#' @param model a \code{\link{MOFAmodel}} object.
#' @param view name of the view
#' @param factors either name (character vector) or indices (numeric vector) of the latent variables for which enrichment should be computed (default="all")
#' @param feature.sets Data structure that holds feature set membership information. Must be either a binary membership matrix (rows are feature sets and columns are features) or a list of feature set member indexes (see vignette).
#' @param local.statistic The feature statistic used to quantify the association between each genomic variable and each PC. Must be one of the following: loading (default), cor, z.
#' @param global.statistic The feature set statisic computed from the feature statistics. Must be one of the following: "mean.diff" (default) or "rank.sum".
#' @param transformation Optional transformation to apply to the feature-level statistics. Must be one of the following "none" or "abs.value" (default).
#' @param statistical.test The statistical test used to compute the significance of the feature set statistics under a competitive null hypothesis.
#' Must be one of the following: "parametric" (default), "cor.adj.parametric", "permutation".
#' @param min.size Minimum size of a feature set
#' @param nperm Number of permutations to perform. Only relevant if statistical.test is set to "permutation".
#' @param cores Number of cores to run the permutation analysis in parallel. Only relevant if statistical.test is set to "permutation".
#' @param p.adj.method Method to adjust p-values factor-wise for multiple testing. Can be any method in p.adjust.methods(). Default uses Benjamini-Hochberg procedure.
#' @details fill this
#' @return a list with two matrices, one for the feature set statistics and the other for the pvalues. Rows are feature sets and columns are latent variables.
#' @import foreach doParallel
#' @importFrom stats p.adjust
#' @export

FeatureSetEnrichmentAnalysis <- function(model, view, feature.sets, factors="all", local.statistic=c("loading", "cor", "z"),
                                         transformation=c("abs.value", "none"), global.statistic=c("mean.diff", "rank.sum"),
                                         statistical.test=c("parametric", "cor.adj.parametric", "permutation"),
                                         nperm=100, min.size=10, cores=1, p.adj.method = "BH", alpha=0.1) {
  
  # Collect factors
  if (paste0(factors,sep="",collapse="") == "all") { 
    factors <- factorNames(model)
    if (model@ModelOpts$learnMean == T) { factors <- factors[-1] }
  } else if(!all(factors %in% factorNames(model))) stop("Factors do not match factor names in model")
  
  # Collect observed data
  data <- model@TrainData[[view]]
  data <- t(data)
  
  # Collect relevant expectations
  W <- getExpectations(model,"SW","E")[[view]][,factors, drop=FALSE]
  Z <- getExpectations(model,"Z","E")[,factors, drop=FALSE]
  
  # Check that there is no constant factor
  stopifnot( all(apply(Z,2,var,na.rm=T)>0) )
  
  # To-do: check feature.sets input format
  # to-do: to reduce FDR problems, extract only factors that are active in that view
  
  # turn feature.sets into binary membership matrices if provided as list
  if(class(feature.sets) == "list") {
    features <- Reduce(union, feature.sets)
    feature.sets <- sapply(names(feature.sets), function(nm) {
      tmp <- features %in% feature.sets[[nm]]
      names(tmp) <- features
      tmp
    })
    feature.sets <-t(feature.sets)*1
  }
  
  # Check if some features do not intersect between the feature sets and the observed data and remove them
  features <- intersect(colnames(data),colnames(feature.sets))
  if(length(features) == 0 ) stop("Feautre names in feature.sets do not match feature names in model.")
  data <- data[,features]
  W <- W[features,]
  feature.sets <- feature.sets[,features]
  
  # Filter feature sets with small number of features
  feature.sets <- feature.sets[rowSums(feature.sets)>=min.size,]
  
  # Match test options
  # local.statistic <- match.arg(local.statistic)
  # transformation <- match.arg(transformation)
  # global.statistic <- match.arg(global.statistic)
  # statistical.test <- match.arg(statistical.test)
  
  # Print options
  message("Doing feature Ontology Enrichment Analysis with the following options...")
  message(sprintf("View: %s", view))
  message(sprintf("Latent variables: %s", paste(as.character(factors),collapse=" ")))
  message(sprintf("Number of feature sets: %d", nrow(feature.sets)))
  message(sprintf("Local statistic: %s", local.statistic))
  message(sprintf("Transformation: %s", transformation))
  message(sprintf("Global statistic: %s", global.statistic))
  message(sprintf("Statistical test: %s", statistical.test))
  if (statistical.test=="permutation") {
    message(sprintf("Cores: %d", cores))
    message(sprintf("Number of permutations: %d", nperm))
  }
  
  # use own version for permutation test because of bugs in PCGSE package
  if (statistical.test == "permutation") {
    doParallel::registerDoParallel(cores=cores)
    
    null_dist_tmp <- foreach(rnd=1:nperm) %dopar% {
      perm <- sample(ncol(data))
      # Permute rows of the weight matrix to obtain a null distribution
      W_null <- apply(W, 2, function(w) w[perm])
      rownames(W_null) <- rownames(W); colnames(W_null) <- colnames(W)
      
      data_null <- data[,perm]
      rownames(data_null) <- rownames(data)
      
      # Compute null statistic
      s.null <- pcgse(data=data_null, prcomp.output=list(rotation=W_null, x=Z), pc.indexes=1:length(factors), feature.sets=feature.sets, feature.statistic=local.statistic,
                      transformation=transformation, feature.set.statistic=global.statistic, feature.set.test="parametric", nperm=NA)$statistic
      abs(s.null)
    }
    null_dist <- do.call("rbind", null_dist_tmp)
    colnames(null_dist) <- factors
    
    # Compute true statistics
    s.true <- pcgse(data=data, prcomp.output=list(rotation=W, x=Z), pc.indexes=1:length(factors), feature.sets=feature.sets, feature.statistic=local.statistic,
                    transformation=transformation, feature.set.statistic=global.statistic, feature.set.test="parametric", nperm=NA)$statistic
    colnames(s.true) <- factors
    rownames(s.true) <- rownames(feature.sets)
    
    # directly use permutation to control FDR
    calcfdr <- function(th, idx) {
      fp <- sum(abs(null_dist[,idx]) > th)
      tp <- sum(abs(s.true[,idx]) > th )
      fp/(tp+fp) }
    sigPathways <- lapply(1:length(factors), function(j) {
      fdr <-sapply(seq(0,100,0.01), calcfdr,j)
      Idxsel <- which(fdr<=alpha)[1]
      if(!is.na(Idxsel)) {
        tsel <- seq(0,100,0.01)[Idxsel]
        rownames(s.true)[s.true[,j] > tsel]
      } else NULL
    })
      
    # Compute p.values
    p.values <- matrix(NA, nrow=nrow(s.true), ncol=ncol(s.true));
    rownames(p.values) <- rownames(s.true); colnames(p.values) <- factors
    for (j in factors) {
      p.values[,j] <- sapply(s.true[,j], function(x) mean(abs(null_dist[,j]) > abs(x)) )
    }

  } else {
    p.values <- pcgse(data=data, prcomp.output=list(rotation=W, x=Z), pc.indexes=1:length(factors), feature.sets=feature.sets, feature.statistic=local.statistic,
                      transformation=transformation, feature.set.statistic=global.statistic, feature.set.test=statistical.test, nperm=nperm)$p.values
    colnames(p.values) <- factors
    rownames(p.values) <- rownames(feature.sets)
  }
  
  if(!p.adj.method %in%  p.adjust.methods) stop("p.adj.method needs to be an element of p.adjust.methods")
  adj.p.values <- apply(p.values, 2,function(lfw) p.adjust(lfw, method = p.adj.method))
  
  sigPathways2 <- lapply(factors, function(j) rownames(adj.p.values)[adj.p.values[,j] <= alpha])
  if(statistical.test != "permutation") sigPathways <- NULL
  
  return(list(pval = p.values, pval.adj = adj.p.values, sigPathways=sigPathways, sigPathways2=sigPathways2))
}


########################
## Plotting functions ##
########################


#' @title Line plot of Feature Set Enrichment Analysis results
#' @name LinePlot_FeatureSetEnrichmentAnalysis
#' @description Line plot of the Feature Set Enrichment Analyisis results for a specific latent variable
#' @param p.values output of \link{FeatureSetEnrichmentAnalysis} function. A data frame of p.values where rows are feature sets and columns are latent variables
#' @param factor Factor for which to show wnriched pathways in the lineplot
#' @param threshold p.value threshold to filter out feature sets
#' @param max.pathways maximum number of enriched pathways to display
#' @details fill this
#' @return nothing
#' @import ggplot2
#' @export
LinePlot_FeatureSetEnrichmentAnalysis <- function(p.values, factor, threshold=0.1, max.pathways=25, ...) {
  
  # Sanity checks
  # (...)
  
  # Get data  
  tmp <- as.data.frame(p.values[,factor, drop=F])
  tmp$pathway <- rownames(tmp)
  colnames(tmp) <- c("pvalue")
  
  # Filter out pathways
  tmp <- tmp[tmp$pvalue<=threshold,,drop=F]
  if(nrow(tmp)==0) {
    warning("No siginificant pathways at the specified threshold. For an overview use Heatmap_FeatureSetEnrichmentAnalysis().")
    return()
  }
  
  # If there are too many pathways enriched, just keep the 'max_pathways' more significant
  if (nrow(tmp) > max.pathways)
    tmp <- head(tmp[order(tmp$pvalue),],n=max.pathways)
  
  # Convert pvalues to log scale (add a small regulariser to avoid numerical errors)
  tmp$log <- -log10(tmp$pvalue + 1e-6)
  
  # Annotate significcant pathways
  # tmp$sig <- factor(tmp$pvalue<threshold)
  
  #order according to significance
  tmp$pathway <- factor(tmp$pathway <- rownames(tmp), levels = tmp$pathway[order(tmp$pvalue, decreasing = T)])
  
  p <- ggplot(tmp, aes(x=pathway, y=log)) +
    # ggtitle(paste("Enriched sets in factor", factor)) +
    geom_point(size=5) +
    geom_hline(yintercept=-log10(threshold), linetype="longdash") +
    scale_y_continuous(limits=c(0,7)) +
    scale_color_manual(values=c("black","red")) +
    geom_segment(aes(xend=pathway, yend=0)) +
    ylab("-log pvalue") +
    coord_flip() +
    theme(
      axis.text.y = element_text(size=rel(1.2), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.2), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}

#' @title Heatmap of Feature Set Enrichment Analysis results
#' @name Heatmap_FeatureSetEnrichmentAnalysis
#' @description Heatmap of the Feature Set Enrichment Analyisis results
#' @param p.values output of \link{FeatureSetEnrichmentAnalysis} function. A data frame of p.values where rows are gene sets and columns are latent variables
#' @param threshold p.value threshold to filter out feature sets. If a feature set has a p.value lower than 'threshold'
#' @param ... Parameters to be passed to pheatmap function
#' @details fill this
#' @return vector of factors being enriched for at least one feautre set at the threshold specified 
#' @import pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
Heatmap_FeatureSetEnrichmentAnalysis <- function(p.values, threshold=0.05, log=T, ...) {
  p.values <- p.values[!apply(p.values, 1, function(x) sum(x>=threshold)) == ncol(p.values),]
  if (log) {
    p.values <- -log10(p.values+0.001)
    threshold <- -log10(threshold)
    col <- colorRampPalette(c("lightgrey", "red"))(n=10)
  } else {
    col <- colorRampPalette(c("red", "lightgrey"))(n=10)
  }
  pheatmap::pheatmap(p.values, color = col, ...)
  # sigFactors <- colnames(p.values)[which(apply(p.values, 2, function(x) any(x<=threshold)))]
  # return(sigFactors)
}





######################################################################
## From here downwards it is a modified version of the PCGSE module ##
######################################################################

pcgse = function(data, 
                 prcomp.output, 
                 pc.indexes=1, 
                 feature.sets,
                 feature.statistic="z",
                 transformation="none",
                 feature.set.statistic="mean.diff",
                 feature.set.test="cor.adj.parametric",
                 nperm=9999 # for feature.set.test value of "permutation"
) {
  current.warn = getOption("warn")
  options(warn=-1)
  # if (is.na(data)) {
  #   stop("'data must' be specified!")
  # }  
  if (is.na(feature.sets)) {
    stop("'feature.sets' must be specified!")
  }   
  options(warn=current.warn) 
  if (!(feature.statistic %in% c("loading", "cor", "z"))) {
    stop("feature.statistic must be 'loading', 'cor' or 'z'")
  }  
  if (!(transformation %in% c("none", "abs.value"))) {
    stop("transformation must be 'none' or 'abs.value'")
  }  
  if (!(feature.set.statistic %in% c("mean.diff", "rank.sum"))) {
    stop("feature.set.statistic must be 'mean.diff' or 'rank.sum'")
  }    
  if (!(feature.set.test %in% c("parametric", "cor.adj.parametric", "permutation"))) {
    stop("feature.set.test must be one of 'parametric', 'cor.adj.parametric', 'permutation'")
  }
  if (feature.set.test == "permutation" & feature.statistic == "loading") { 
    stop("feature.statistic cannot be set to 'loading' if feature.set.test is 'permutation'")
  }
  if (!is.matrix(feature.sets) & feature.set.test == "permutation") {
    stop("feature.sets must be specified as a binary membership matrix if feature.set.test is set to 'permutation'") 
  }  
  # if (feature.set.test == "parametric") {
  #   warning("The 'parametric' test option ignores the correlation between feature-level test statistics and therefore has an inflated type I error rate. ",
  #           "This option should only be used for evaluation purposes.")    
  # }  
  # if (feature.set.test == "permutation") {
  #   warning("The 'permutation' test option can be extremely computationally expensive given the required modifications to the safe() function. ",
  #           "For most applications, it is recommended that feature.set.test is set to 'cor.adj.parametric'.")
  # }
  
  # Turn the feature set matrix into list form if feature.set.test is not "permutation"
  feature.set.indexes = feature.sets  
  if (is.matrix(feature.sets)) {
    feature.set.indexes = createVarGroupList(var.groups=feature.sets)  
  }
  
  n = nrow(data)
  p = ncol(data)
  
  # Compute the feature-level statistics.
  feature.statistics = matrix(0, nrow=p, ncol=length(pc.indexes))
  for (i in 1:length(pc.indexes)) {
    pc.index = pc.indexes[i]
    feature.statistics[,i] = computefeatureStatistics(data=data, prcomp.output=prcomp.output, pc.index=pc.index, feature.statistic, transformation)
  }
  
  # Perform the specified feature set test for each feature set on each specified PC using the feature-level statistics
  if (feature.set.test == "parametric" | feature.set.test == "cor.adj.parametric") {
    if (feature.set.statistic == "mean.diff") {
      results = pcgseViaTTest(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes,
                              feature.statistics=feature.statistics, cor.adjustment=(feature.set.test == "cor.adj.parametric"))      
    } else if (feature.set.statistic == "rank.sum") {
      results = pcgseViaWMW(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes,
                            feature.statistics=feature.statistics, cor.adjustment=(feature.set.test == "cor.adj.parametric"))
    }     
  } else if (feature.set.test == "permutation") {
    # results = pcgseViaSAFE(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes, 
    #                        feature.statistic=feature.statistic, transformation=transformation, feature.set.statistic=feature.set.statistic, nperm=nperm)
    results = pcgseViaPermutation(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes, 
                                  feature.statistics=feature.statistics, feature.set.statistic=feature.set.statistic, nperm=nperm)        
  }
  
  return (results) 
}




#-------------------------------------------------------------------------------------------------------------------------------
# Internal methods - General
#-------------------------------------------------------------------------------------------------------------------------------

#
# Turn the annotation matrix into a list of var group indexes for the valid sized var groups
#
createVarGroupList = function(var.groups) {
  var.group.indexes = list()  
  for (i in 1:nrow(var.groups)) {
    member.indexes = which(var.groups[i,]==1)
    var.group.indexes[[i]] = member.indexes    
  }
  names(var.group.indexes) = rownames(var.groups)    
  return (var.group.indexes)
}

# 
# Computes the feature-level statistics
#
computefeatureStatistics = function(data, prcomp.output, pc.index, feature.statistic, transformation) {
  p = ncol(data)
  n = nrow(data)
  feature.statistics = rep(0, p)
  if (feature.statistic == "loading") {
    # get the PC loadings for the selected PCs
    feature.statistics = prcomp.output$rotation[,pc.index]
  } else {
    # compute the Pearson correlation between the selected PCs and the data
    feature.statistics = cor(data, prcomp.output$x[,pc.index], use = "complete.obs") 
    if (feature.statistic == "z") {
      # use Fisher's Z transformation to convert to Z-statisics
      feature.statistics = sapply(feature.statistics, function(x) {
        return (sqrt(n-3)*atanh(x))})      
    }    
  }
  
  # Absolute value transformation of the feature-level statistics if requested
  if (transformation == "abs.value") {
    feature.statistics = sapply(feature.statistics, abs)
  }  
  
  return (feature.statistics)
}

#-------------------------------------------------------------------------------------------------------------------------------
# Internal methods - Enrichment via t-test or correlation-adjusted t-test
#-------------------------------------------------------------------------------------------------------------------------------

pcgseViaTTest = function(data, prcomp.output, pc.indexes, feature.set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets = length(feature.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(feature.set.indexes)
  feature.set.statistics = matrix(T, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(feature.set.statistics) = names(feature.set.indexes)    
  
  for (i in 1:num.feature.sets) {
    indexes.for.feature.set = feature.set.indexes[[i]]
    m1 = length(indexes.for.feature.set)
    not.feature.set.indexes = which(!(1:ncol(data) %in% indexes.for.feature.set))
    m2 = length(not.feature.set.indexes)
    
    if (cor.adjustment) {      
      # compute sample correlation matrix for members of feature set
      cor.mat = cor(data[,indexes.for.feature.set], use = "complete.obs")
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
      # compute the VIF, using CAMERA formula from Wu et al., based on Barry et al.
      vif = 1 + (m1 -1)*mean.cor
    }
    
    for (j in 1:length(pc.indexes)) {
      # get the feature-level statistics for this PC
      pc.feature.stats = feature.statistics[,j]
      # compute the mean difference of the feature-level statistics
      mean.diff = mean(pc.feature.stats[indexes.for.feature.set]) - mean(pc.feature.stats[not.feature.set.indexes])
      # compute the pooled standard deviation
      pooled.sd = sqrt(((m1-1)*var(pc.feature.stats[indexes.for.feature.set]) + (m2-1)*var(pc.feature.stats[not.feature.set.indexes]))/(m1+m2-2))      
      # compute the t-statistic
      if (cor.adjustment) {
        t.stat = mean.diff/(pooled.sd*sqrt(vif/m1 + 1/m2))
        df = n-2
      } else {
        t.stat = mean.diff/(pooled.sd*sqrt(1/m1 + 1/m2))
        df = m1+m2-2
      }
      feature.set.statistics[i,j] = t.stat      
      # compute the p-value via a two-sided test
      lower.p = pt(t.stat, df=df, lower.tail=T)
      upper.p = pt(t.stat, df=df, lower.tail=F)        
      p.values[i,j] = 2*min(lower.p, upper.p)      
    }
  } 
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = feature.set.statistics  
  
  return (results)
}

#-------------------------------------------------------------------------------------------------------------------------------
# Internal methods - Enrichment via Wilcoxon Mann Whitney or correlation-adjusted WMW
#-------------------------------------------------------------------------------------------------------------------------------

pcgseViaWMW = function(data, prcomp.output, pc.indexes, feature.set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets = length(feature.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(feature.set.indexes)
  feature.set.statistics = matrix(T, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(feature.set.statistics) = names(feature.set.indexes)    
  
  for (i in 1:num.feature.sets) {
    indexes.for.feature.set = feature.set.indexes[[i]]
    m1 = length(indexes.for.feature.set)
    not.feature.set.indexes = which(!(1:ncol(data) %in% indexes.for.feature.set))
    m2 = length(not.feature.set.indexes)
    
    if (cor.adjustment) {            
      # compute sample correlation matrix for members of feature set
      cor.mat = cor(data[,indexes.for.feature.set])
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
    }
    
    for (j in 1:length(pc.indexes)) {
      # get the feature-level statistics for this PC
      pc.feature.stats = feature.statistics[,j]
      # compute the rank sum statistic feature-level statistics
      wilcox.results = wilcox.test(x=pc.feature.stats[indexes.for.feature.set], y=pc.feature.stats[not.feature.set.indexes],
                                   alternative="two.sided", exact=F, correct=F)
      rank.sum = wilcox.results$statistic                
      if (cor.adjustment) {
        # Using correlation-adjusted formula from Wu et al.
        var.rank.sum = ((m1*m2)/(2*pi))*(asin(1) + (m2 - 1)*asin(.5) + (m1-1)*(m2-1)*asin(mean.cor/2) +(m1-1)*asin((mean.cor+1)/2))
      } else {        
        var.rank.sum = m1*m2*(m1+m2+1)/12
      }
      z.stat = (rank.sum - (m1*m2)/2)/sqrt(var.rank.sum)
      feature.set.statistics[i,j] = z.stat      
      # compute the p-value via a two-sided z-test
      lower.p = pnorm(z.stat, lower.tail=T)
      upper.p = pnorm(z.stat, lower.tail=F)        
      p.values[i,j] = 2*min(lower.p, upper.p)
    }
  } 
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = feature.set.statistics  
  
  return (results)
}
