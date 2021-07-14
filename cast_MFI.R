library(gtools)
library(plot.matrix)
library(readxl)
library(jsonlite)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

source("acoustic-eDNA/classification_functions.R")

# MFI and classification visuals from - Gering 07/2021

# Color palette
c_palette <- c("darkorchid4", "green", "orange", "red")

# Import Json Data
json_data <- jsonlite::fromJSON("/Volumes/GeringSSD/Subset_Sv/segment_subset_CTD14-15.json")

# Empty data frame to hold summary of MFI information 
summary_df <- data.frame(double(), double(), double(), double(), double(), double(), integer(), integer())
colnames(summary_df) <- names(summary(0)) %>% append("Cast") %>% append("Depth")

# loops through each csv in the json data
for (n_ctd in seq_along(json_data)){
  ctd_name <- names(json_data)[[n_ctd]]
  print(ctd_name)
  ctd <- json_data[[n_ctd]]
  # for each depth in the CTD cast
  MFI_classes <- data.frame(Classification = character(), Count = double(), Depth = integer())
  for (n_depth in seq_along(ctd)){
    depth_name <- names(ctd)[[n_depth]]
    print(depth_name)
    
    # list of data arrays for each frequency at this depth
    depth <- ctd[[n_depth]]
    # calculate MFI
    MFI = calc_MFI(depth)
    
    # plot MFI (uncomment jpeg and dev.off to save)
    jpeg(paste0("/Volumes/GeringSSD/MFI/", ctd_name, "_", depth_name, "_MFI.jpg"), width = 900, height = 800)
    plot(t(MFI), breaks=c(-0.02, 0.4, 0.6, 0.8, 1.0), # some values are very slightly smaller that 0
         col=c_palette, 
         main= paste("CTD:", ctd_name, "Depth: ", depth_name, "- MFI"))
    dev.off()
    
    # histogram
    #jpeg(paste0("/Volumes/GeringSSD/MFI/", ctd_name, "_", depth_name, "_MFI_hist.jpg"), width = 900, height = 800)
    hist(MFI, col = "blue", main = paste("CTD:", ctd_name, "Depth: ", depth_name, "- Histogram"))
    #dev.off()
    
    # binning
    depth_class <- as.data.frame(table(classify_MFI(MFI))) %>% mutate(Depth = strtoi(depth_name))
    colnames(depth_class) <- c("Classification", "Count", "Depth")
    MFI_classes <- rbind(MFI_classes, depth_class)
    
    # add summary information to df dataframe
    MFI_summary <- summary(as.vector(MFI))
    temp_df <- data.frame(x=t(matrix(MFI_summary)))
    colnames(temp_df) <- names(MFI_summary)
    temp_df <- temp_df %>% mutate(Depth = strtoi(depth_name)) %>% mutate(Cast = strtoi(sub("^ctd0+", "", ctd_name)))
    summary_df <- rbind(summary_df, temp_df)
  }
  
  #jpeg(paste0("/Volumes/GeringSSD/MFI/", ctd_name, "_MFI_classifications.jpg"), width = 900, height = 800)
  par(mar=c(4,10,4,4))
  barplot(height = MFI_classes$Count, names = MFI_classes$Classification, density = MFI_classes$Depth, angle = MFI_classes$Depth,
          horiz = T, las=2, col=c_palette,
          main = paste("CTD:", ctd_name, "- Classifications"))
  abline(h=c(4.9 , 9.7) , col="black", lwd = 2.5)
  
  legend("topright", legend = unique(MFI_classes$Depth), 
         density = unique(MFI_classes$Depth), angle = unique(MFI_classes$Depth), 
         bty = "n", pt.cex = 2, cex = 2, horiz = FALSE)
  #dev.off()
}
# save df
#write.csv(df,"/Volumes/GeringSSD/MFI/MFI_subsets_CTD14-15.csv", row.names = FALSE)

