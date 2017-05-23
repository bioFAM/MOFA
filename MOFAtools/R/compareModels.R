
#' @title Compare the latent factors of different \code{\link{MOFAmodel}} 
#' @name compareModels
#' @description Different objects of \code{\link{MOFAmodel} are compared in terms of correlation between 
#' their latent factors. The correlation is calculated only on those samples which are present in all models 
#' used for a comparison. 
#' @param ModelList a list containing \code{\link{MOFAmodel}} objects.
#' @param comparison tye of comparison either 'pairwise' or 'all'
#' @details asd
#' @return Plots a heatmap of correlation of Latent Factors in all models when 'comparison' is 'all'. 
#' Otherwise, for each pair of models, a seperate heatmap is produced comparing one model againt the other.
#' The corresponding correlation matrix or list or pairwise correlation matrices is returned
#' @reference fill this
#' @import pheatmap
#' @export

compareModels <- function(ModelList, comparison="all") {
  #check inputs
  if(class(ModelList)!="list")
    stop("'ModelList' has to be a list")
  if (!all(sapply(ModelList, function (l) class(l)=="MOFAmodel")))
    stop("List elements of 'ModelList' have to be an instance of MOFAmodel")
  if (!comparison %in% c("all", "pairwise"))
    stop("'comparison' has to be either 'all' or 'pairwise'")
  
  # give generic names if no names present
  if(is.null(names(ModelList))) names(ModelList) <- paste("model", 1: length(ModelList), sep="")
  
  #get latent factors
  LFs <- lapply(seq_along(ModelList), function(modelidx){
    model <- ModelList[[modelidx]]
    Z <- getExpectations(model, 'Z', 'E')
    if(model@ModelOpts$learnMean) Z <- Z[,-1]
    if(is.null(rownames(Z))) rownames(Z) <- rownames(model@TrainData[[1]])
    if(is.null(colnames(Z))) 
      if(model@ModelOpts$learnMean) colnames(Z) <- paste("LF", 2:(ncol(Z)+1), sep="") else colnames(Z) <- paste("LF", 1:ncol(Z), sep="")
    Z
    })
  for(i in seq_along(LFs)) 
    colnames(LFs[[i]]) <- paste(names(ModelList)[i], colnames(LFs[[i]]), sep="_")
  
    if(comparison=="all"){
    #get common samples between models
    commonSamples <- Reduce(intersect,lapply(LFs, rownames))
    if(is.null(commonSamples)) 
      stop("No common samples in all models for comparison")
    
    #subset LFs to common samples
    LFscommon <- Reduce(cbind, lapply(LFs, function(Z) Z[commonSamples,]))

    # calculate correlation
    corLFs <- cor(LFscommon)
    
    #annotation by model
    modelAnnot <- data.frame(model = rep(names(ModelList), times=sapply(LFs, ncol)))
    rownames(modelAnnot) <- colnames(LFscommon)
    
    #plot heatmap
    pheatmap(abs(corLFs), show_rownames = F,
             color=colorRampPalette(c("white", "orange" ,"red"))(100), 
             annotation_col = modelAnnot, main= "Absolute correlation between latent factors")
    
    return(corLFs)
    }
  
  if(comparison=="pairwise"){
    PairWiseCor <- lapply(seq_along(LFs[-length(LFs)]), function(i){
      LFs1<-LFs[[i]]
        sublist <- lapply((i+1):length(LFs), function(j){
          LFs2<-LFs[[j]]
          common_pairwise <- intersect(rownames(LFs1), rownames(LFs2))
          if(is.null(common_pairwise)) {
            warning(paste("No common samples between models",i,"and",j,"- No comparison possible"))
            NA
          }
          else{
          corLFs_pairs <- cor(LFs1[common_pairwise,], LFs2[common_pairwise,])
          pheatmap(abs(corLFs_pairs),color=colorRampPalette(c("white", "orange" ,"red"))(100),
                   main=paste("Absolute correlation between factors in model", i,"and",j))
          corLFs_pairs
          }
        })
        names(sublist) <- names(ModelList)[(i+1):length(LFs)]
    })
    names(PairWiseCor) <- names(ModelList[-length(ModelList)])
    return(PairWiseCor)
  }
}