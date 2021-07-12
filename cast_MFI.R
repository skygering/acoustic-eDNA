library(gtools)
library(plot.matrix)
library(readxl)
library(jsonlite)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Calculate MFI from .json subsetted Sv data - Gering 07/2021

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
      D_f[,,fq] = sv_f[,,fq]/max(sv_f[,,fq], na.rm=TRUE)
    }
  }
  if (scale_max == 1) { # scale to global maximum
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

# Color palette
c_palette <- c("blue", "green", "orange", "red")

# Import Json Data
json_data <- jsonlite::fromJSON("/Volumes/GeringSSD/Subset_Sv/segment_subset_CTD14-15.json")

# Empty data frame to hold summary o f MFI information 
df <- data.frame(double(), double(), double(), double(), double(), double(), integer(), integer())
colnames(df) <- names(summary(0)) %>% append("Cast") %>% append("Depth")

# loops through each csv in the json data
for (n_ctd in seq_along(json_data)){
  ctd_name <- names(json_data)[[n_ctd]]
  print(ctd_name)
  ctd <- json_data[[n_ctd]]
  # for each depth in the CTD cast
  MFI_bars <- data.frame(Classification = character(), Count = double(), Depth = integer())
  for (n_depth in seq_along(ctd)){
    depth_name <- names(ctd)[[n_depth]]
    print(depth_name)
    # list of data arrays for each frequency at this depth
    depth <- ctd[[n_depth]]
    # calculate MFI
    MFI = calc_MFI(depth)
    # plot MFI (uncomment jpeg and dev.off to save)
    jpeg(paste0("/Volumes/GeringSSD/MFI/", ctd_name, "_", depth_name, "_MFI.jpg"), width = 900, height = 800)
    plot(t(MFI), breaks=c(0, 0.4, 0.6, 0.8, 1.0), 
         col=c_palette, 
         main= paste("CTD:", ctd_name, "Depth: ", depth_name, "- MFI"))
    dev.off()
    
    # histogram
    jpeg(paste0("/Volumes/GeringSSD/MFI/", ctd_name, "_", depth_name, "_MFI_hist.jpg"), width = 900, height = 800)
    hist(MFI, col = "blue", main = paste("CTD:", ctd_name, "Depth: ", depth_name, "- Histogram"))
    dev.off()
    
    # binning
    bins <- c(0, 0.4, 0.6, 0.8, 1.0)
    catagories <- c("Swimbladder Fish", "Small Bubbles", "Zooplankton", "Non-swimbladder Fish")
    temp_bars <- as.data.frame(table(cut(MFI, breaks = bins, labels = catagories))) %>% mutate(Depth = strtoi(depth_name))
    colnames(temp_bars) <- c("Classification", "Count", "Depth")
    MFI_bars <- rbind(MFI_bars, temp_bars)
    
    # add summary information to df dataframe
    MFI_summary <- summary(as.vector(MFI))
    temp_df <- data.frame(x=t(matrix(MFI_summary)))
    colnames(temp_df) <- names(MFI_summary)
    temp_df <- temp_df %>% mutate(Depth = strtoi(depth_name)) %>% mutate(Cast = strtoi(sub("^ctd0+", "", ctd_name)))
    df <- rbind(df, temp_df)
  }
  
  jpeg(paste0("/Volumes/GeringSSD/MFI/", ctd_name, "_MFI_classifications.jpg"), width = 900, height = 800)
  par(mar=c(4,10,4,4))
  barplot(height = MFI_bars$Count, names = MFI_bars$Classification, density = MFI_bars$Depth, angle = MFI_bars$Depth,
          horiz = T, las=2, col=c_palette, 
          main = paste("CTD:", ctd_name, "- Classifications"))
  abline(h=c(4.9 , 9.7) , col="black", lwd = 2.5)
  
  legend("topright", legend = unique(MFI_bars$Depth), 
         density = unique(MFI_bars$Depth), angle = unique(MFI_bars$Depth), 
         bty = "n", pt.cex = 2, cex = 2, horiz = FALSE)
  dev.off()
}
# save df
write.csv(df,"/Volumes/GeringSSD/MFI/MFI_subsets_CTD14-15.csv", row.names = FALSE)

