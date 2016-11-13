#######################################
## Gene set enrichment analysis ##
#######################################

# Gene statistics: the gene-level statistic used to quantify the association between
# each genomic variable and each latent variable (default is "loading")
#   "loading": the weight associated with the genomic variable.

# Transformation: Optional transformation to apply to the gene-level statistics (default is "none)
#   "none": no transformation
#   "abs.value": absolute value

# Gene set statistics: the gene set statisic used to quantify the enrichment of a latent variable for a given gene set  (default is "mean.diff")
#   "mean.diff": standardized difference between the mean of the gene-level statistics for members of the gene set and the mean
#             of the gene-level statistics for genomic variables not in the gene set. 
#   "rank.sum": standardized Wilcoxon rank sum statistic computed from the gene-level statistics for members of the gene set. 

# Gene set tests available: the statistical test used to compute the significance of the gene set statistics under a 
#             competitive null hypothesis (default is "cor.adj.parametric")

#   "parametric": 
#      if the mean difference is being used as the gene set statistic, corresponds to a two-sided, two-sample t-test with equal variances.
#      if the rank sum is being used as the gene set statistic, this corresponds to a two-sided, 
#       two-sample z-test based on the standardized rank sum statistic. 
#       NOTE: both of these tests incorrectly assume the gene-level statistics are iid and should therefore only be used for comparative purposes.

#   "cor.adj.parametric": tests statistical significance of the standardized and correlation-adjusted gene set statistic using a two-sided t-test or z-test.  
#             Similar to the CAMERA method by Wu et al., standardization of either the mean different statistic or rank sum statistic is performed   
#             using a variation inflation factor based on the average pair-wise correlation between the gene-level statistics for members of the gene set. 
#             Per Barry et al., this is approximated by the average correlation between the genomic variables. 
#             Per Wu et al., the significance of the correlation-adjusted t-statistic is tested using a two-sided t-test with n-2 df and 
#             the significance of the correlation-adjusted rank sum z-statistic is tested using a two-sided z-test. 

#   "permutation": tests gene set enrichment via the permutation distribution of the gene set statistic.
#             The permutation distribution is computed via permutation of the sample labels, which, in this case, is equivalent to permutation 
#             of the elements of the target PC. This test is realized using the safe() function from the R safe package.
#             The number of permutations is controlled by the "nperm" parameter. The gene.statistic cannot be set to "loadings" with this option.
#             of the gene set. Per Barry et al., this correlation is approximated by the average correlation between the genomic variables.
#   nperm: Number of permutations to perform. Only relevant if gene.set.test is "permutation".
#

# Inputs:
#   data (numeric matrix or df): empirical data matrix with dim (samples,features)
#   W (numerix matrix): weight matrix with dim (features,factors)
#   Z (numerix matrix): latent variable matrix with dim (samples,factors)
#   factor.indexes (vector of int): index of the factors to be tested
#   gene.sets (binary matrix): gene set membership matrix, dim (sets,features)
#   gene.statistic
#   transformation:
#   gene.set.statistic
#   gene.set.test
#   nperm
gsea = function(data, Z, W, factor.indexes=1, sets, feature.statistic,
                     transformation, set.statistic, test, nperm) {
  
  # precompute some terms
  # D = nrow(W)
  # K = ncol(W)
  
  ## Step 1: calculate gene set statistics ##
  statistics = matrix(0, nrow=D, ncol=length(factor.indexes))
  for (i in 1:length(factor.indexes))
    statistics[,i] = computeFeatureStatistics(data, Z, W, factor.indexes[i], feature.statistic, transformation)
  

  ## Step 2: test statistical significance of gene sets ##
  
  # Parametric test
  if (set.test == "parametric") {
    
    # Two sample T-test using the means
    if (set.statistic == "mean.diff") {
      results = pcgseViaTTest(data=data, pc.indexes=pc.indexes, gene.set.indexes=gene.set.indexes,
                              gene.statistics=gene.statistics, cor.adjustment=(gene.set.test == "cor.adj.parametric"))  
    #   
    # # Wilcoxon rank sum statistic
    # } else if (gene.set.statistic == "rank.sum") {
    #   results = pcgseViaWMW(data=data, pc.indexes=pc.indexes, gene.set.indexes=gene.set.indexes,
    #                         gene.statistics=gene.statistics, cor.adjustment=(gene.set.test == "cor.adj.parametric"))
    # }     
    
  # Non-parametric permutation test
  } else if (gene.set.test == "permutation") {
    stop()
    results = pcgseViaSAFE(data=data, Z=Z, W=W, pc.indexes=pc.indexes, gene.set.indexes=gene.set.indexes, 
                           gene.statistic=gene.statistic, transformation=transformation, gene.set.statistic=gene.set.statistic, nperm=nperm)        
  }
  return (statistics)
}

# Function to compute the feature statistics
# Input:
# Output:  
computeFeatureStatistics = function(data, Z, W, factor.index, statistic, transformation) {
  
  if (statistic == "loading") {
    # Use the weight loadings directly as statistics
    statistics = W[,factor.index]
  } else {
    # compute the Pearson correlation between the selected factor and the data
    statistics = cor(data, Z[,factor.index]) 
    if (statistic == "z") {
      # use Fisher's Z transformation to convert to Z-statisics
      statistics = sapply(statistics, function(x) { return (sqrt(n-3)*atanh(x))})      
    }
  }

  # Absolute value transformation of the statistics
  if (transformation == "abs.value")
    statistics = sapply(statistics, abs)
  
  return (statistics)
}

# Function to compute enrichment via fisher's exact test if given a binary matrix
# Inputs:
# - data:
# - factor.indexes: 
# - set.indexes: list where element i contains all column indexes of features belonging to category i
# - feature.statistics: matrix with dim (features, factors) with the feature-level statistics
# - cor.adjustment: bool indicating whether to adjust by correlation

FisherTest <- function(data, factor.indexes, set.indexes) {
  # fisher.test(rbind(c(1,9),c(11,3)), alternative="less")$p.value
}

# Function to compute enrichment via t-test or correlation-adjusted t-test
# Inputs:
# - data:
# - factor.indexes: 
# - set.indexes: list where element i contains all column indexes of features belonging to category i
# - feature.statistics: matrix with dim (features, factors) with the feature-level statistics
# - cor.adjustment: bool indicating whether to adjust by correlation
# Output:
# - p.values: matrix with dim (ngene_sets, npcs) where element ij contains the 
#     pvalue for pc j being enriched for class i 
# - gene.set.statistics: matrix with dim (ngene_sets,npcs) where element ij contains the 
#     gene set statistic for pc j being enriched for class i 
Ttest <- function(data, factor.indexes, set.indexes, feature.statistics, cor.adjustment) {

  num.gene.sets = length(gene.set.indexes)
  n = nrow(data)
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
