# Function to get preprocessed drug response data
# inputs:
# file: location of drpar object
# pat2include: patients ids to include
# conc2include: concentrations to include
# badDrugs drugs to be excluded ---- #default drugs that failes QC: NSC 74859, bortezomib.
# targetedDrugs drugs to be included for sure ---- #default drugs "ibrutinib", "idelalisib",  "PRT062607 HCl", "duvelisib", "spebrutinib", "selumetinib", "MK-2206",  "everolimus", "encorafenib"
# conc4targeted concentration for targeted drugs to include ---- default 4,5
# chemoDrugs rugs to be included for sure ---- #default drugs "fludarabine", "doxorubicine",  "nutlin-3"
# conc4targeted concentration for chemotherapeutical drugs to include ---- default 3,4,5
#  
# Threshold parameters: drugs are accepted that for at least `effectNum` samples 
# have a viability effect less than or equal to `effectVal`. On the other hand, the 
# mean viability effect should not be below `viab`.
#
#
#
# outdir, directory to save output to

stripConc <- function (x) 
  vapply(strsplit(x, "_"), function(x) paste(x[-length(x)], collapse = "_"), 
         character(1))

deckel <- function(x, lower = -Inf, upper = +Inf) ifelse(x<lower, lower, ifelse(x>upper, upper, x))

getViab<-function(pat2include, 
                  badDrugs=c( "D_008",  "D_025"), 
                  conc2include = 2:5,
                  targetedDrugs= c("D_002", "D_003", "D_166", "D_082", "D_079", "D_012", "D_030", "D_063", "D_083") , 
                  conc4targeted = c(4,5),
                  chemoDrugs = c("D_006", "D_159", "D_010"),
                  conc4chemo = 3:5,
                  effectNum = 4,
                  effectVal = 0.7,
                  viab = 0.6, 
                  maxval = 1.1,
                  outdir,
                  plotit = T,
                  verbose = T){
  
  #Load object
  data("lpdAll", package="PACEdata")
  dr<-lpdAll[fData(lpdAll)$type=="viab"]
  
  
  #Subset to patients of interest
  dr<-dr[, colnames(dr) %in% pat2include]
  
  #Select drug fulfilling requirements
  
  ##Filter out bad drugs and concentrations
  candDrugs <- rownames(dr)[ !(fData(dr)$id %in% badDrugs) & fData(dr)$subtype %in% conc2include]

  # Drugs to be included for sure
  targetedDrugs2include <- paste(rep(targetedDrugs, each = length(conc4targeted)), conc4targeted, sep="_" )
  chemoDrugs2include <- paste(rep(chemoDrugs, each = length(conc4chemo)), conc4chemo, sep="_" )

  ##Thresholds on viability effect
  overallMean  <- rowMeans(exprs(dr)[candDrugs, ])
  nthStrongest <- apply(exprs(dr)[candDrugs, ], 1, function(x) sort(x)[effectNum])
  eligibleDrugs <- candDrugs[ overallMean >=viab & nthStrongest <= effectVal ] %>%
    union(targetedDrugs2include) %>% union(chemoDrugs2include)
  if(plotit){
  par(mfrow = c(1, 3))
  hist(overallMean,  breaks = 30, col = "pink"); abline(v = viab,      col="blue")
  hist(nthStrongest, breaks = 30, col = "pink"); abline(v =effectVal, col="blue")
  plot(overallMean, nthStrongest)
  abline(h = effectVal, v = viab, col = "blue")
  }  

  if(verbose){
  message("Including p= ", length(eligibleDrugs), " drug response features")
  message("Different drugs d = ",  length(unique(stripConc(eligibleDrugs))))
  message("Concentration per drug:")
  print(table(stripConc(eligibleDrugs)))
  }
  
  #subset
  dr<-dr[eligibleDrugs,, drop=FALSE ]
  
  #cut-off all viability balues above maxval
  drmat <- t(deckel(exprs(dr), lower=0, upper=maxval))

  #Save as view
  save(drmat, file=file.path(outdir,"viab.RData"))
  # write.table(drmat, file=file.path(outdir,"viab.txt"),
  #             row.names=TRUE, col.names=TRUE, quote=F)

  return(drmat)
}
