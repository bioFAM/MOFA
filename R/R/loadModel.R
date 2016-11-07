# Function to load a trained model using python
library(RcppCNPy)
library(stringr)

# Function to load training statistics


# Function to load model expectations
load_npy <- function(infolder) {
  files = str_replace( list.files(path=infolder,pattern="*.npy"), ".npy", "" )
  data = list()
  vars = unique( sapply(files, function(f) strsplit(f, "\\_")[[1]][1]) )
  for (v in vars) {
    v_files = files[str_detect(files,str_c(v,"_"))]
    params = unique( sapply(v_files, function(f) strsplit(f, "\\_")[[1]][2]) )
    data[[v]] <- sapply(params, function(x) NULL)
    for (p in params) {
      p_files = v_files[str_detect(v_files,p)]
      views = as.numeric( unique( sapply(p_files, function(f) strsplit(f, "\\_")[[1]][3]) ) )
      if (all(is.na(views))) {
        data[[v]][[p]] = npyLoad(sprintf("%s/%s_%s.npy",infolder,v,p))
      } else {
        data[[v]][[p]] = vector("list",length=length(views))
        for (m in views) {
          data[[v]][[p]][[m]] = npyLoad(sprintf("%s/%s_%s_%s.npy",infolder,v,p,m))
        }
      }
    }
  }
  return(data)
}

infolder="/Users/ricard/git/scGFA/scMT/tmp/model"
# indir="/tmp/test"
data = load_npy(infolder)

# npyLoad(sprintf("%s/activeK.npy",indir))

meta = read.table("/Users/ricard/git/scGFA/scMT/e_metadata.txt")
