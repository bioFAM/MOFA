##################################
## Gene set enrichment analysis ##
##################################

############
## Theory ##
############

# Gene statistics: the gene-level statistic used to quantify the association between
# each genomic variable and each latent variable (default is "loading")
#   "loading": the weight associated with the genomic variable.

# Transformation: Optional transformation to apply to the gene-level statistics (default is "none")
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
#      if the rank sum is being used as the gene set statistic, this corresponds to a two-sided, two-sample z-test based on the standardized rank sum statistic. 
#      NOTE: both of these tests incorrectly assume the gene-level statistics are iid and should therefore only be used for comparative purposes.

# COMPLETE THIS
#   "cor.adj.parametric": tests statistical significance of the standardized and correlation-adjusted gene set statistic 
#      Similar to the CAMERA method by Wu et al., standardization of either the mean difference statistic or rank sum statistic is performed   
#      using a variation inflation factor based on the average pair-wise correlation between the gene-level statistics for members of the gene set. 
#      Per Barry et al., this is approximated by the average correlation between the genomic variables. 
#      Per Wu et al., the significance of the correlation-adjusted t-statistic is tested using a two-sided t-test with n-2 df and 
#      the significance of the correlation-adjusted rank sum z-statistic is tested using a two-sided z-test. 

# COMPLETE THIS
#   "permutation": tests gene set enrichment via the permutation distribution of the gene set statistic.
#             The permutation distribution is computed via permutation of the sample labels This test is realized using the safe() function from the R safe package.
#             The number of permutations is controlled by the "nperm" parameter. The gene.statistic cannot be set to "loadings" with this option.
#             of the gene set. Per Barry et al., this correlation is approximated by the average correlation between the genomic variables.


#   nperm: Number of permutations to perform. Only relevant if gene.set.test is "permutation".
#

# Simulate a test data set
N = 100
D = 50
K = 10
M = 20

data <- matrix(rnorm(N*D),N,D)
W <- matrix(rnorm(D*K),D,K)
Z <- matrix(rnorm(N*K),N,K)
sets <- matrix(sample(c(0,1),size=M*D,replace=T,prob=c(0.8,0.2)),M,D)
local.statistic="loading"
transformation="abs.value"
global.statistic="mean.diff"
statistical.test="parametric"
nperm=NA
factor.indexes=c(1,2,3)

# Inputs:
#   data (numeric matrix): empirical data matrix with dim (samples,features)
#   W (numerix matrix): weight matrix with dim (features,factors)
#   Z (numerix matrix): latent variable matrix with dim (samples,factors)
#   factor.indexes (vector of int): index of the factors to be tested
#   gene.sets (binary matrix): gene set membership matrix, dim (sets,features)
#   gene.statistic
#   transformation:
#   gene.set.statistic
#   gene.set.test
#   nperm
gsea = function(data, Z, W, factor.indexes=1, sets, local.statistic,
                     transformation, global.statistic, statistical.test, nperm) {
  
  # Turn the set matrix into list form
  sets.index <- apply(sets, 1, function(set) which(set==1) )

  # Precompute some terms
  N = nrow(Z)
  D = nrow(W)
  K = length(factor.indexes)
  M = nrow(sets)
  
  ## Step 1: calculate local statistics ##
  local.stats = matrix(0, nrow=D, ncol=length(factor.indexes))
  for (i in 1:length(factor.indexes))
    local.stats[,i] = computeLocalStatistics(data, Z, W, factor.indexes[i], local.statistic, transformation)
  
  ## Step 2: calculate global statistics ##
  # global.stats = matrix(0, nrow=M, ncol=length(factor.indexes))
  # for (i in 1:length(factor.indexes))
    # global.stats[,i] = computeGlobalStatistics(local.stats, factor.indexes[i], global.statistic, sets.index)
    
  ## Step 2: calculate global statistics and test statistical significance of sets ##
  
  # Parametric test
  if (statistical.test == "parametric") {
    
    # Two sample T-test using the means
    if (set.statistic == "mean.diff") {
      # results = Ttest(data=data, pc.indexes=pc.indexes, gene.set.indexes=gene.set.indexes,
                              # gene.statistics=gene.statistics, cor.adjustment=(gene.set.test == "cor.adj.parametric"))  
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

# Function to compute the local statistics
# Input:
# Output:  
computeLocalStatistics = function(data, Z, W, factor.index, local.statistic, transformation) {
  
  if (local.statistic == "loading") {
    # Use the weight loadings directly as statistics
    stats = W[,factor.index]
  } else {
    # compute the Pearson correlation between the selected factor and the data
    stats = cor(data, Z[,factor.index]) 
    if (statistic == "z") {
      # use Fisher's Z transformation to convert to Z-statisics
      N = nrow(data)
      stats = sapply(stats, function(x) { return (sqrt(N-3)*atanh(x))})      
    }
  }

  # Absolute value transformation of the statistics
  if (transformation == "abs.value")
    stats = sapply(stats, abs)
  
  return (stats)
}


# Function to compute the global statistics
# Input:
# Output:  
# computeGlobalStatistics = function(local.stats, factor.index, global.statistic, sets.index) {
#   
#   tmp = local.stats[,factor.index]
#   if (global.statistic == "mean.diff") {
#       # stats = apply(sets, 1, function(set) mean(tmp[which(set==1)]) - mean(tmp[which(set==0)]) )
#       stats = sapply(sets.index, function(i) mean(tmp[i]) - mean(tmp[!(1:length(tmp) %in% i)]) )
#   }
#   if (global.statistic == "rank.sum") {
#       # stats = apply(sets, 1, function(set) wilcox.test(tmp[which(set==1)],tmp[which(set==0)],"two.sided") )  
#       stats = sapply(sets.index, 1, function(i) wilcox.test(tmp[i],tmp[!(1:length(tmp) %in% i)],"two.sided") )  
#   }
#   return (stats)
# }


# Function to compute enrichment via fisher's exact test if given a binary matrix
# Inputs:
# - data:
# - factor.indexes: 
# - set.indexes: list where element i contains all column indexes of features belonging to category i
# - feature.statistics: matrix with dim (features, factors) with the feature-level statistics
# - cor.adjustment: bool indicating whether to adjust by correlation

# FisherTest <- function(data, factor.indexes, set.indexes) {
#   # fisher.test(rbind(c(1,9),c(11,3)), alternative="less")$p.value
# }

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
Ttest <- function(data, factor.indexes, sets.index, local.stats, global.stats, cor.adjustment) {

  M <- length(sets.index)
  D <- nrow(local.stats)
  K <- length(factor.indexes)
  
  # Define the pvalue and t statistic matrix
  p.values <- matrix(0, nrow=M, ncol=K)  
  test.statistics <- matrix(0, nrow=M, ncol=K)  
  rownames(p.values) <- names(sets.index)
  rownames(test.statistics) <- names(sets.index)
  
  for (m in 1:M) {
    i <- sets.index[[m]]
    j <- which(!(1:ncol(data) %in% i))
    m1 <- length(i)
    m2 <- length(j)
    for (k in 1:length(factor.indexes)) {

      # Compute the mean difference as global statistic
      IS DIM CORRECT?
      mean.cor <- sapply(sets.index, function(i) mean(local.stats[i,k]) - mean(local.stats[!(1:D %in% i)]) )
      
      # Adjust the global statistic taking into account the intra-set correlation between features
      if (cor.adjustment) {
        # compute sample correlation matrix for members of gene set
        cor.mat <- cor(data[,i])
        # compute the mean pair-wise correlation 
        mean.cor <- (sum(cor.mat) - m1)/(m1*(m1-1))    
        # compute the VIF, using CAMERA formula from Wu et al., based on Barry et al.
        vif <- 1 + (m1 -1)*mean.cor
      }
      
      
      mean.diff <- global.stats[m,k]
      
      # compute the pooled standard deviation
      pooled.sd = sqrt(((m1-1)*var(local.stats[i]) + (m2-1)*var(local.stats[j]))/(m1+m2-2))
      
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
  }
  
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
