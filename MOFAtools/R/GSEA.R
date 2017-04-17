
#' @title Gene Set Enrichment Analysis
#' @name GSEA
#' @description Method to perform gene ontology enrichment for a gene sets in a given view  \cr
#' Here we use the PCGSE package.
#' @param model a \code{\link{MOFAmodel}} object.
#' @param view name (character) or index (integer) of the view of interest
#' @param factor.indexes numeric vector with the indices of the latent variables for which enrichment should be computed.
#' @param gene.sets Data structure that holds gene set membership information. Must be either a binary membership matrix (rows are gene sets and columns are genes) or a list of gene set member indexes.
#' @param local.statistic The gene-level statistic used to quantify the association between each genomic variable and each PC. Must be one of the following: loading (default), cor, z.
#' @param global.statistic The gene set statisic computed from the gene-level statistics. Must be one of the following: "mean_diff" (default) or "rank.sum".
#' @param transformation Optional transformation to apply to the gene-level statistics. Must be one of the following "none" (default) or "abs.value".
#' @param statistical.test The statistical test used to compute the significance of the gene set statistics under a competitive null hypothesis.  
#' Must be one of the following: "parametric", "cor.adj.parametric" (default), "permutation".
#' @param nperm Number of permutations to perform. Only relevant if statistical.test is "permutation".
#' @details fill this
#' @return a list with two matrices, one for statistics and one for pvalues. The rows are the different gene sets and the columns are the latent variables.
#' @import safe
#' @export

GSEA <- function(model, view, factor.indexes, gene.sets, local.statistic="loading",
                  transformation="none", global.statistic="mean.diff", statistical.test="cor.adj.parametric", nperm=1, min_size=10) {
  
  if (paste0(factor.indexes,sep="",collapse="") == "all") { 
    factor.indexes <- 1:model@Dimensions[["K"]]
  }
  
  # Collect expectations
  data <- model@TrainData[[view]]
  W <- getExpectations(model,"SW","E")[[view]]; rownames(W) <- colnames(data)
  Z <- getExpectations(model,"Z","E");          rownames(Z) <- rownames(data)
  
  # Remove genes that are not present in either the data or the annotation file
  genes <- intersect(colnames(data),colnames(gene.sets))
  data <- data[,genes]
  W <- W[genes,]
  gene.sets <- gene.sets[,genes]
  
  # Filter gene sets with small number of genes
  gene.sets <- gene.sets[rowSums(gene.sets)>=min_size,]
  

  cat("Doing Gene Ontology Enrichment Analysis with the following options...\n")
  cat(sprintf("View: %s\n", view))
  cat(sprintf("Latent variables: %s\n", paste(as.character(factor.indexes),collapse=" ")))
  cat(sprintf("Number of gene sets: %d\n", nrow(gene.sets)))
  cat(sprintf("Local statistic: %s\n", local.statistic))
  cat(sprintf("Transformation: %s\n", transformation))
  cat(sprintf("Global statistic: %s\n", global.statistic))
  cat(sprintf("Statistical test: %s\n", statistical.test))
  
  p <- pcgse(data=data, prcomp.output=list(rotation=W, x=Z), pc.indexes=factor.indexes, gene.sets=gene.sets, gene.statistic=local.statistic,
        transformation=transformation, gene.set.statistic=global.statistic, gene.set.test=statistical.test, nperm=nperm)$p.values
  colnames(p) <- factorNames(model)[factor.indexes]
  rownames(p) <- rownames(gene.sets)
  
  return(p)
}



######################################################################################################
######################################################################################################
######################################################################################################

pcgse = function(data, 
                 prcomp.output=NA, 
                 pc.indexes=1, 
                 gene.sets,
                 gene.statistic="z",
                 transformation="none",
                 gene.set.statistic="mean.diff",
                 gene.set.test="cor.adj.parametric",
                 nperm=9999 # for gene.set.test value of "permutation"
) {
  current.warn = getOption("warn")
  options(warn=-1)
  if (is.na(data)) {
    stop("'data must' be specified!")
  }  
  if (is.na(gene.sets)) {
    stop("'gene.sets' must be specified!")
  }   
  options(warn=current.warn) 
  if (!(gene.statistic %in% c("loading", "cor", "z"))) {
    stop("gene.statistic must be 'loading', 'cor' or 'z'")
  }  
  if (!(transformation %in% c("none", "abs.value"))) {
    stop("transformation must be 'none' or 'abs.value'")
  }  
  if (!(gene.set.statistic %in% c("mean.diff", "rank.sum"))) {
    stop("gene.set.statistic must be 'mean.diff' or 'rank.sum'")
  }    
  if (!(gene.set.test %in% c("parametric", "cor.adj.parametric", "permutation"))) {
    stop("gene.set.test must be one of 'parametric', 'cor.adj.parametric', 'permutation'")
  }
  if (gene.set.test == "permutation" & gene.statistic == "loading") { 
    stop("gene.statistic cannot be set to 'loading' if gene.set.test is 'permutation'")
  }
  if (!is.matrix(gene.sets) & gene.set.test == "permutation") {
    stop("gene.sets must be specified as a binary membership matrix if gene.set.test is set to 'permutation'") 
  }  
  if (gene.set.test == "parametric") {
    warning("The 'parametric' test option ignores the correlation between gene-level test statistics and therefore has an inflated type I error rate. ",
            "This option should only be used for evaluation purposes.")    
  }  
  # if (gene.set.test == "permutation") {
  #   warning("The 'permutation' test option can be extremely computationally expensive given the required modifications to the safe() function. ",
  #           "For most applications, it is recommended that gene.set.test is set to 'cor.adj.parametric'.")
  # }
  
  # Turn the gene set matrix into list form if gene.set.test is not "permutation"
  gene.set.indexes = gene.sets  
  if (is.matrix(gene.sets) & gene.set.test != "permutation") {
    gene.set.indexes = createVarGroupList(var.groups=gene.sets)  
  }
  
  # Compute PCA if necessary
  if (is.na(prcomp.output)) {
    # Center and scale so that PCA is performed on correlation matrix
    prcomp.output=prcomp(data, scale=T)
  }
  
  n = nrow(data)
  p = ncol(data)
  
  # If gene.set.test is not "permutation", compute the gene-level statistics.
  if (gene.set.test == "permutation") {
    gene.statistics = NA
  } else {
    gene.statistics = matrix(0, nrow=p, ncol=length(pc.indexes))
    for (i in 1:length(pc.indexes)) {
      pc.index = pc.indexes[i]
      gene.statistics[,i] = computeGeneStatistics(data=data, prcomp.output=prcomp.output, pc.index=pc.index, gene.statistic, transformation)
    }
  }
  
  # Perform the specified gene set test for each gene set on each specified PC using the gene-level statistics
  if (gene.set.test == "parametric" | gene.set.test == "cor.adj.parametric") {
    if (gene.set.statistic == "mean.diff") {
      results = pcgseViaTTest(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, gene.set.indexes=gene.set.indexes,
                              gene.statistics=gene.statistics, cor.adjustment=(gene.set.test == "cor.adj.parametric"))      
    } else if (gene.set.statistic == "rank.sum") {
      results = pcgseViaWMW(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, gene.set.indexes=gene.set.indexes,
                            gene.statistics=gene.statistics, cor.adjustment=(gene.set.test == "cor.adj.parametric"))
    }     
  } else if (gene.set.test == "permutation") {
    results = pcgseViaSAFE(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, gene.set.indexes=gene.set.indexes, 
                           gene.statistic=gene.statistic, transformation=transformation, gene.set.statistic=gene.set.statistic, nperm=nperm)        
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
# Computes the gene-level statistics
#
computeGeneStatistics = function(data, prcomp.output, pc.index, gene.statistic, transformation) {
  p = ncol(data)
  n = nrow(data)
  gene.statistics = rep(0, p)
  if (gene.statistic == "loading") {
    # get the PC loadings for the selected PCs
    gene.statistics = prcomp.output$rotation[,pc.index]
  } else {
    # compute the Pearson correlation between the selected PCs and the data
    gene.statistics = cor(data, prcomp.output$x[,pc.index], use = "complete.obs") 
    if (gene.statistic == "z") {
      # use Fisher's Z transformation to convert to Z-statisics
      gene.statistics = sapply(gene.statistics, function(x) {
        return (sqrt(n-3)*atanh(x))})      
    }    
  }
  
  # Absolute value transformation of the gene-level statistics if requested
  if (transformation == "abs.value") {
    gene.statistics = sapply(gene.statistics, abs)
  }  
  
  return (gene.statistics)
}

#-------------------------------------------------------------------------------------------------------------------------------
# Internal methods - Enrichment via t-test or correlation-adjusted t-test
#-------------------------------------------------------------------------------------------------------------------------------

pcgseViaTTest = function(data, prcomp.output, pc.indexes, gene.set.indexes, gene.statistics, cor.adjustment) {
  
  num.gene.sets = length(gene.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.gene.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(gene.set.indexes)
  gene.set.statistics = matrix(T, nrow=num.gene.sets, ncol=length(pc.indexes))    
  rownames(gene.set.statistics) = names(gene.set.indexes)    
  
  for (i in 1:num.gene.sets) {
    indexes.for.gene.set = gene.set.indexes[[i]]
    m1 = length(indexes.for.gene.set)
    not.gene.set.indexes = which(!(1:ncol(data) %in% indexes.for.gene.set))
    m2 = length(not.gene.set.indexes)
    
    if (cor.adjustment) {      
      # compute sample correlation matrix for members of gene set
      cor.mat = cor(data[,indexes.for.gene.set], use = "complete.obs")
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
      # compute the VIF, using CAMERA formula from Wu et al., based on Barry et al.
      vif = 1 + (m1 -1)*mean.cor
    }
    
    for (j in 1:length(pc.indexes)) {
      # get the gene-level statistics for this PC
      pc.gene.stats = gene.statistics[,j]
      # compute the mean difference of the gene-level statistics
      mean.diff = mean(pc.gene.stats[indexes.for.gene.set]) - mean(pc.gene.stats[not.gene.set.indexes])
      # compute the pooled standard deviation
      pooled.sd = sqrt(((m1-1)*var(pc.gene.stats[indexes.for.gene.set]) + (m2-1)*var(pc.gene.stats[not.gene.set.indexes]))/(m1+m2-2))      
      # compute the t-statistic
      if (cor.adjustment) {
        t.stat = mean.diff/(pooled.sd*sqrt(vif/m1 + 1/m2))
        df = n-2
      } else {
        t.stat = mean.diff/(pooled.sd*sqrt(1/m1 + 1/m2))
        df = m1+m2-2
      }
      gene.set.statistics[i,j] = t.stat      
      # compute the p-value via a two-sided test
      lower.p = pt(t.stat, df=df, lower.tail=T)
      upper.p = pt(t.stat, df=df, lower.tail=F)        
      p.values[i,j] = 2*min(lower.p, upper.p)      
    }
  } 
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = gene.set.statistics  
  
  return (results)
}

#-------------------------------------------------------------------------------------------------------------------------------
# Internal methods - Enrichment via SAFE method that computes the permutation distribution of the mean difference of t-statistics
# associated with a linear regression of the PC on the genomic variables.
#-------------------------------------------------------------------------------------------------------------------------------

setGeneric("local.GeneStatistics", function(X.mat, y.vec,...) {standardGeneric("local.GeneStatistics")})
setMethod("local.GeneStatistics",
          c(X.mat = "matrix", y.vec = "numeric"),
          function(X.mat, y.vec, ...) {
            args.local = list(...)[[1]]
            gene.stat = args.local$gene.statistic  
            trans = args.local$transformation    
            return (function(data=X.mat, vector=y.vec, gene.statistic = gene.stat, transformation = trans,...) {  
              p = length(vector)
              n = ncol(data) 
              gene.statistics = rep(0, p)
              # compute the Pearson correlation between the selected PCs and the data
              gene.statistics = cor(t(data), vector, use = "complete.obs") 
              if (gene.statistic == "z") {
                # use Fisher's Z transformation to convert to Z-statisics
                gene.statistics = sapply(gene.statistics, function(x) {
                  return (sqrt(n-3)*atanh(x))})      
              }    
              # Absolute value transformation of the gene-level statistics if requested
              if (transformation == "abs.value") {
                gene.statistics = sapply(gene.statistics, abs)
              }
              return (gene.statistics)
            })
          }
)

# local.GeneStatistics = function(X.mat, y.vec,...){ 
#   args.local = list(...)[[1]]
#   gene.stat = args.local$gene.statistic  
#   trans = args.local$transformation    
#   return (function(data=X.mat, vector=y.vec, gene.statistic = gene.stat, transformation = trans,...) {  
#     p = length(vector)
#     n = ncol(data) 
#     gene.statistics = rep(0, p)
#     # compute the Pearson correlation between the selected PCs and the data
#     gene.statistics = cor(t(data), vector) 
#     if (gene.statistic == "z") {
#       # use Fisher's Z transformation to convert to Z-statisics
#       gene.statistics = sapply(gene.statistics, function(x) {
#         return (sqrt(n-3)*atanh(x))})      
#     }    
#     # Absolute value transformation of the gene-level statistics if requested
#     if (transformation == "abs.value") {
#       gene.statistics = sapply(gene.statistics, abs)
#     }
#     return (gene.statistics)
#   })
# }

setGeneric("global.StandAveDiff", function(C.mat,local.stats,...) { standardGeneric("global.StandAveDiff") })
setMethod("global.StandAveDiff",
          c(C.mat = "matrix", local.stats = "numeric"),
          function(C.mat,local.stats,...) {
            gene.set.mat = t(as.matrix(C.mat))
            return (function(gene.statistics=local.stats,gene.sets=gene.set.mat,...) {  
              num.gene.sets = nrow(gene.set.mat)
              gene.set.statistics = rep(0, num.gene.sets)
              for (i in 1:num.gene.sets) {
                indexes.for.gene.set = which(gene.set.mat[i,]==1)
                m1 = length(indexes.for.gene.set)
                not.gene.set.indexes = which(!(1:ncol(gene.set.mat) %in% indexes.for.gene.set))
                m2 = length(not.gene.set.indexes)
                # compute the mean difference of the gene-level statistics
                mean.diff = mean(gene.statistics[indexes.for.gene.set]) - mean(gene.statistics[not.gene.set.indexes])
                # compute the pooled standard deviation
                pooled.sd = sqrt(((m1-1)*var(gene.statistics[indexes.for.gene.set]) + (m2-1)*var(gene.statistics[not.gene.set.indexes]))/(m1+m2-2))      
                # compute the t-statistic
                gene.set.statistics[i] = mean.diff/(pooled.sd*sqrt(1/m1 + 1/m2))
              }
              return (gene.set.statistics)
            })
          }
)
          
# global.StandAveDiff = function(C.mat,local.stats,...) { 
  # gene.set.mat = t(as.matrix(C.mat))
  # return (function(gene.statistics=local.stats,gene.sets=gene.set.mat,...) {  
  #   num.gene.sets = nrow(gene.set.mat)
  #   gene.set.statistics = rep(0, num.gene.sets)
  #   for (i in 1:num.gene.sets) {
  #     indexes.for.gene.set = which(gene.set.mat[i,]==1)
  #     m1 = length(indexes.for.gene.set)
  #     not.gene.set.indexes = which(!(1:ncol(gene.set.mat) %in% indexes.for.gene.set))
  #     m2 = length(not.gene.set.indexes)
  #     # compute the mean difference of the gene-level statistics
  #     mean.diff = mean(gene.statistics[indexes.for.gene.set]) - mean(gene.statistics[not.gene.set.indexes])
  #     # compute the pooled standard deviation
  #     pooled.sd = sqrt(((m1-1)*var(gene.statistics[indexes.for.gene.set]) + (m2-1)*var(gene.statistics[not.gene.set.indexes]))/(m1+m2-2))      
  #     # compute the t-statistic
  #     gene.set.statistics[i] = mean.diff/(pooled.sd*sqrt(1/m1 + 1/m2))
  #   }
  #   return (gene.set.statistics)
  # })
# }

pcgseViaSAFE = function(data, prcomp.output, pc.indexes, gene.set.indexes, gene.statistic, transformation, gene.set.statistic, nperm=999) {
  num.gene.sets = nrow(gene.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.gene.sets, ncol=length(pc.indexes))
  rownames(p.values) = names(gene.set.indexes)
  gene.set.statistics = matrix(T, nrow=num.gene.sets, ncol=length(pc.indexes))    
  rownames(gene.set.statistics) = names(gene.set.indexes)  
  
  for (j in 1:length(pc.indexes)) {
    pc.index = pc.indexes[j]
    pc = prcomp.output$x[,pc.index]
    global="AveDiff"
    if (gene.set.statistic == "rank.sum") {
      global="Wilcoxon"
    } else {
      # global="StandAveDiff"
      global="AveDiff"
    }
    safe.results = safe::safe(X.mat=t(data), y.vec=pc, C.mat=t(gene.set.indexes), local="GeneStatistics", global=global, Pi.mat=nperm,
                        args.global=list(one.sided=T), args.local=list(gene.statistic=gene.statistic, transformation=transformation), alpha=1.01, print.it = FALSE)
    p.values[,j] = slot(safe.results, "global.pval")
    gene.set.statistics[,j] = slot(safe.results, "global.stat")    
  }
  
  # Turn one-sided p-values into two-sided p-values
  p.values = apply(p.values, c(1,2), function(x) {
    upper = 1 - x
    lower = x
    return (2*min(upper, lower))
  })
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = gene.set.statistics  
  
  return (results)
}

#-------------------------------------------------------------------------------------------------------------------------------
# Internal methods - Enrichment via Wilcoxon Mann Whitney or correlation-adjusted WMW
#-------------------------------------------------------------------------------------------------------------------------------

pcgseViaWMW = function(data, prcomp.output, pc.indexes, gene.set.indexes, gene.statistics, cor.adjustment) {
  
  num.gene.sets = length(gene.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.gene.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(gene.set.indexes)
  gene.set.statistics = matrix(T, nrow=num.gene.sets, ncol=length(pc.indexes))    
  rownames(gene.set.statistics) = names(gene.set.indexes)    
  
  for (i in 1:num.gene.sets) {
    indexes.for.gene.set = gene.set.indexes[[i]]
    m1 = length(indexes.for.gene.set)
    not.gene.set.indexes = which(!(1:ncol(data) %in% indexes.for.gene.set))
    m2 = length(not.gene.set.indexes)
    
    if (cor.adjustment) {            
      # compute sample correlation matrix for members of gene set
      cor.mat = cor(data[,indexes.for.gene.set])
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
    }
    
    for (j in 1:length(pc.indexes)) {
      # get the gene-level statistics for this PC
      pc.gene.stats = gene.statistics[,j]
      # compute the rank sum statistic gene-level statistics
      wilcox.results = wilcox.test(x=pc.gene.stats[indexes.for.gene.set], y=pc.gene.stats[not.gene.set.indexes],
                                   alternative="two.sided", exact=F, correct=F)
      rank.sum = wilcox.results$statistic                
      if (cor.adjustment) {
        # Using correlation-adjusted formula from Wu et al.
        var.rank.sum = ((m1*m2)/(2*pi))*(asin(1) + (m2 - 1)*asin(.5) + (m1-1)*(m2-1)*asin(mean.cor/2) +(m1-1)*asin((mean.cor+1)/2))
      } else {        
        var.rank.sum = m1*m2*(m1+m2+1)/12
      }
      z.stat = (rank.sum - (m1*m2)/2)/sqrt(var.rank.sum)
      gene.set.statistics[i,j] = z.stat      
      # compute the p-value via a two-sided z-test
      lower.p = pnorm(z.stat, lower.tail=T)
      upper.p = pnorm(z.stat, lower.tail=F)        
      p.values[i,j] = 2*min(lower.p, upper.p)
    }
  } 
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = gene.set.statistics  
  
  return (results)
}




