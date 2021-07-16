library(readxl)
library(jsonlite)
library(dplyr)

# Calculate MFI from .json subsetted Sv data and classifying functions - Gering 07/2021

calc_MFI <- function(ctd_depth, delta = 40, scale_max = 1){
  # create 3D array of Sv data - one slice per frequency
  nd = ncol(ctd_depth[[1]]) # number of depth bins
  np = nrow(ctd_depth[[1]]) # number of pings
  nf = length(ctd_depth) # number of frequencies
  Sv_f = array(0, dim=(c(np, nd, nf)))
  
  for (i in 1:nf){
    Sv_f[, , i] = ctd_depth[[i]]
  }
  # linear data
  sv_f = 10^(Sv_f/10)
  
  # frequencies
  f <- unlist(lapply(names(ctd_depth), function(x) strtoi(x)/1000))
  e_f = 1/f
  
  # all unique combinations of indices
  f_idx = combinations(n=nf, r=2, set=F, repeats.allowed=F)
  d_f = cbind(f_idx, 0)
  
  # Distance function
  for (i in 1:length(f_idx[,1])) {
    d_f[i,3] = 1-exp(-abs(f[f_idx[i,1]]-f[f_idx[i,2]])/delta)
  }
  
  #scaled linear data
  D_f = sv_f
  
  if (scale_max == 0) { # scale to local (per frequency) maximum
    for (fq in 1:nf) {
      print(max(sv_f[,,fq], na.rm=TRUE))
      D_f[,,fq] = sv_f[,,fq]/max(sv_f[,,fq], na.rm=TRUE)
    }
  }
  if (scale_max == 1) { # scale to global maximum
    print(max(sv_f, na.rm=TRUE))
    D_f = sv_f/max(sv_f, na.rm=TRUE)
  }
  
  # calculate the MFI values
  MFI = array(0, dim=(c(np, nd)))
  num = array(0, dim=(c(np, nd)))
  den = array(0, dim=(c(np, nd)))
  for (i in 1:(nf-1)) {
    for (l in (i+1):nf) {
      num = num+
        d_f[which(d_f[,1] == i & d_f[,2] == l),3]*D_f[,,i]*D_f[,,l]*e_f[i]*e_f[l]
      den = den+D_f[,,i]*D_f[,,l]*e_f[i]*e_f[l]
    }
  }
  MFI = ((num/den)-0.4)/0.6
  return(MFI)
}

classify_MFI <- function(x){
  bins <- c(-Inf, 0.4, 0.6, 0.8, Inf)  # some values are very slightly smaller that 0
  catagories <- c("Swimbladder Fish", "Small Bubbles", "Zooplankton", "Non-swimbladder Fish")
  class_MFI <- cut(x, breaks = bins, labels = catagories)
  return(class_MFI)
}
