
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
#' @import PCGSE
#' @export

GSEA <- function(model, view, factor.indexes, gene.sets, local.statistic="loading",
                  transformation="none", global.statistic="mean.diff", statistical.test="cor.adj.parametric", nperm=1, min_size=10) {
  
  if (factor.indexes == "all") {
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
        transformation=transformation, gene.set.statistic=global.statistic, gene.set.test=statistical.test, nperm=nperm)
  
  return(p$p.values)
}
