
################################################
## Functions to compare different MOFA models ##
################################################


#' @title Correlation of the latent factors across different trials
#' @name compareFactors
#' @description Different objects of \code{\link{MOFAmodel}} are compared in terms of correlation between 
#' their latent factors. The correlation is calculated only on those samples which are present in all models.
#' Ideally, the output should look like a block diagonal matrix, suggesting that all detected factors are robust under different initialisations.
#' If not, it suggests that some factors are weak and not captured by all models.
#' @param models a list containing \code{\link{MOFAmodel}} objects.
#' @param comparison tye of comparison, either 'pairwise', i.e. compare one model with another one at a time, or 'all', i.e. calculate correlation between factors from all model. By default, all models are compared.
#' @param ... extra arguments passed to pheatmap
#' @details This function can be helpful to evaluate the robustness of factors across different random initilizations. 
#' Large block of factors from different models in the correlation matrix show consistent factors, while stand-alone factors that are only recovered in a single model instance are less reliable.
#' @return Plots a heatmap of correlation of Latent Factors in all models when 'comparison' is 'all'. 
#' Otherwise, for each pair of models, a seperate heatmap is produced comparing one model againt the other.
#' The corresponding correlation matrix or list or pairwise correlation matrices is returned
#' @references fill this
#' @importFrom stats cor
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export

#' @examples

#' ### Example on simulated data
#' #Simulate Data
#' MOFAobject <- makeExampleData()
#' 
#' # Prepare MOFA
#' TrainOptions <- getDefaultTrainOptions()
#' ModelOptions <- getDefaultModelOptions(MOFAobject)
#' DataOptions <- getDefaultDataOptions()
#' MOFAobject <- prepareMOFA(MOFAobject, DataOptions = DataOptions, ModelOptions = ModelOptions, TrainOptions = TrainOptions)
#' 
#' # Train MOFA multiple times
#' n_inits <- 3 
#' MOFAlist <- lapply(1:n_inits, function(it) runMOFA(MOFAobject, outfile=tempfile()))
#' compareFactors(MOFAlist, comparison="all")
#' compareFactors(MOFAlist, comparison="pairwise")


compareFactors <- function(models, comparison = "all", show_rownames=FALSE, show_colnames=FALSE,...) {
  
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="MOFAmodel")))
    stop("Each element of the the list 'models' has to be an instance of MOFAmodel")
  if (!comparison %in% c("all", "pairwise"))
    stop("'comparison' has to be either 'all' or 'pairwise'")
  
  # give generic names if no names present
  if(is.null(names(models))) names(models) <- paste("model", 1: length(models), sep="")
  
  # get latent factors
  LFs <- lapply(seq_along(models), function(modelidx){
    model <- models[[modelidx]]
    Z <- getFactors(model)
    if (model@ModelOptions$learnIntercept==TRUE) 
      Z <- Z[,-1, drop=FALSE]
    Z
    })
  
  for (i in seq_along(LFs)) 
    colnames(LFs[[i]]) <- paste(names(models)[i], colnames(LFs[[i]]), sep="_")
  
  if (comparison=="all") {
    #get common samples between models
    commonSamples <- Reduce(intersect,lapply(LFs, rownames))
    if(is.null(commonSamples)) 
      stop("No common samples in all models for comparison")
    
    #subset LFs to common samples
    LFscommon <- Reduce(cbind, lapply(LFs, function(Z) {
      Z <- Z[commonSamples,, drop=FALSE]
      nonconst <- apply(Z,2,var, na.rm=TRUE) > 0
      if(sum(nonconst) < ncol(Z)) message("Removing ", sum(!nonconst), " constant factors from the comparison.")
      Z[, nonconst]
    })
      )

    # calculate correlation
    corLFs <- cor(LFscommon, use="complete.obs")
    
    #annotation by model
    # modelAnnot <- data.frame(model = rep(names(models), times=sapply(LFs, ncol)))
    # rownames(modelAnnot) <- colnames(LFscommon)
    
    #plot heatmap
    # if(is.null(main)) main <- "Absolute correlation between latent factors"
    if(length(unique(as.numeric(abs(corLFs))))>1){
    pheatmap(abs(corLFs), show_rownames = show_rownames,show_colnames = show_colnames,
             color = colorRampPalette(c("white",RColorBrewer::brewer.pal(9,name="YlOrRd")))(100),
             # color=colorRampPalette(c("white", "orange" ,"red"))(100), 
             # annotation_col = modelAnnot, main= main , ...)
             ...)
    } else warning("No plot produced as correlations consist of only one value")
    # return(corLFs)
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
          # if(is.null(main)) main <- paste("Absolute correlation between factors in model", i,"and",j)
          nonconst1 <- apply(LFs1[common_pairwise,],2,var, na.rm=TRUE) > 0
          nonconst2 <- apply(LFs2[common_pairwise,],2,var, na.rm=TRUE) > 0
          if(sum(!nonconst1) + sum(!nonconst2)>0) message("Removing ", sum(!nonconst1), " and " ,sum(!nonconst2), " constant factors from the each model for the comparison, respectively.")
          corLFs_pairs <- cor(LFs1[common_pairwise,nonconst1], LFs2[common_pairwise,nonconst2], use="complete.obs")
          if(length(unique(abs(corLFs_pairs)))>1){
          pheatmap(abs(corLFs_pairs),color=colorRampPalette(c("white", "orange" ,"red"))(100), ...)
          } else warning("No plot produced as correlations consist of only one value")
          corLFs_pairs
          }
        })
        names(sublist) <- names(models)[(i+1):length(LFs)]
        sublist
    })
    
    names(PairWiseCor) <- names(models[-length(models)])
    return(NULL)
    #return(PairWiseCor)
  }
}


#' @title Compare different instances of trained \code{\link{MOFAmodel}} 
#' @name compareModels
#' @description Different objects of \code{\link{MOFAmodel}} are compared in terms of the final value of the ELBO statistics. 
#' For model selection the model with the highest ELBO value is selected. The height of the bar indicates the number of inferred factors and the color of the bar the value of the ELBO statistic.
#' @param models a list containing \code{\link{MOFAmodel}} objects.
#' @param show_modelnames boolean, whether to indicate the name of each model instance (names of the list in models) or not
#' @export
#' @examples
#' ### Example on simulated data
#' #Simulate Data
#' MOFAobject <- makeExampleData()
#' 
#' # Prepare MOFA
#' TrainOptions <- getDefaultTrainOptions()
#' ModelOptions <- getDefaultModelOptions(MOFAobject)
#' DataOptions <- getDefaultDataOptions()
#' MOFAobject <- prepareMOFA(MOFAobject, DataOptions = DataOptions, ModelOptions = ModelOptions, TrainOptions = TrainOptions)
#' 
#' # Train MOFA multiple times
#' n_inits <- 3 
#' MOFAlist <- lapply(1:n_inits, function(it) runMOFA(MOFAobject, outfile=tempfile()))
#' compareModels(MOFAlist)

compareModels <- function(models, show_modelnames = FALSE) {
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="MOFAmodel")))
    stop("Each element of the the list 'models' has to be an instance of MOFAmodel")

  elbo_vals <- sapply(models, getELBO)
  n_factors <- sapply(models, function(m) {
    n_fac <- getDimensions(m)$K
    if(m@ModelOptions$learnIntercept) n_fac <- n_fac - 1
    n_fac
    })
  if(is.null(names(models))) names(models) <- paste0("model_", seq_along(models))
  df <- data.frame(ELBO=elbo_vals, model = names(models))
  gg <- ggplot(df, aes(x=model,y=n_factors, fill=ELBO)) + geom_bar(stat="identity")+
    ylab("Number of inferred factors")
  if(show_modelnames) gg <- gg + theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
  else gg <- gg + theme(axis.text.x=element_blank())
  return(gg)
}

#' @title Select the best model from a list of trained \code{\link{MOFAmodel}} objects
#' @name selectModel
#' @description Different trained objects of \code{\link{MOFAmodel}} are compared in terms of the final value of the ELBO statistics 
#' and the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{MOFAmodel}} objects.
#' @param plotit show a plot of the characteristics of the compared  \code{\link{MOFAmodel}} objects (ELBO value and number of inferred factors)?
#' @export
#' @examples
#' ### Example on simulated data
#' #Simulate Data
#' MOFAobject <- makeExampleData()
#' 
#' # Prepare MOFA
#' TrainOptions <- getDefaultTrainOptions()
#' ModelOptions <- getDefaultModelOptions(MOFAobject)
#' DataOptions <- getDefaultDataOptions()
#' MOFAobject <- prepareMOFA(MOFAobject, DataOptions = DataOptions, ModelOptions = ModelOptions, TrainOptions = TrainOptions)
#' 
#' # Train MOFA multiple times
#' n_inits <- 3 
#' MOFAlist <- lapply(1:n_inits, function(it) runMOFA(MOFAobject, outfile=tempfile()))
#' selectModel(MOFAlist)

selectModel <- function(models, plotit =TRUE) {
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="MOFAmodel")))
    stop("Each element of the the list 'models' has to be an instance of MOFAmodel")

  elbo_vals <- sapply(models, getELBO)
  if(plotit) print(compareModels(models))
  models[[which.max(elbo_vals)]]
}
