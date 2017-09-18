# Function to get covariates
# inputs:
# file: location of patmeta object
# pat2include: patients ids to include
# cov2include: covariates to be included
# outdir, directory to save output to


getCovariates<-function(pat2include, cov2include, outdir){
  #Load object
  data("patmeta", package = "PACEdata")
  meta <- patmeta
 
  #subset
  meta <- meta[pat2include,cov2include, drop=F]
  
  #Save
    save(meta, file=file.path(outdir,"covariates.RData"))
    # write.table(meta, file=file.path(outdir,"meta.txt"),
    #             row.names=TRUE, col.names=TRUE, quote=F)
    
    return(meta)
  }
