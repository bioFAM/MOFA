# Function to get preprocessed meth data
# inputs:
# file: location of meth object
# pat2include: patients ids to include
# Perc2include= percentage of top variable CpG sites to be included
# outdir, directory to save output to


getMeth<-function(pat2include, Frac2include=0.3, outdir, includeXYchr=T, methDataFile=NULL){
  #Load object
  data("methData", package = "PACEdata")
  meth <- methData
  
  #Subset to patients of interest
  meth <- meth[, colnames(meth)%in% pat2include]
  meth <- assay(meth) %>% t
  
  #transform to M-values
  meth <- log2((meth+0.001)/(1-meth+0.001))

  #Optional: Filter out sites on the sex chromosomes
  if(!includeXYchr){
  methMeta <- fread(methDataFile)
  CpGOnXY <- filter(methMeta, chr %in% c("chrX", "chrY"))
  meth<-meth[, !colnames(meth)%in% CpGOnXY$cg]
  }
  
  #Filter out CpG sitew with low variance, only keep top Perc2include percent
  methVars<-apply(meth, 2, var)
  meth<- meth[, methVars>quantile(methVars, 1-Frac2include)]

  #Save
    save(meth, file=file.path(outdir,"meth.RData"))
    # write.table(meth, file=file.path(outdir,"meth.txt"),
                # row.names=TRUE, col.names=TRUE, quote=F)
    
    return(meth)
  }
