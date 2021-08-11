# Scripts that reads in subset JSON file and subset_bounds JSON file (both created by python code) and calculates MFI and ABC values
# Script also creates plot of MFI values, a historgram of MFI values for each cast, 
# and violin plot of log(ABC)  as well as calculating Kendall's Tau significance test
# Developer: Skylar Gering 08/2021

library(readxl)
library(jsonlite)
library(gtools)
library(plot.matrix)
library(ggplot2)
library(RColorBrewer)
library(colorspace)

source("acoustic-eDNA/Rcode/classification_functions.R")

# Color palette - matches colors used in Python code, but lightened
c_palette <- c(lighten("#006164", 0.2), lighten("#57c4ad", 0.2), lighten("#eda247", 0.2), lighten("#db4325", 0.2))

# Import Subset Json Data
subset_data <- jsonlite::fromJSON("/Volumes/GeringSSD/GU201905_Output/subsets.json")
subset_bounds <- jsonlite::fromJSON("/Volumes/GeringSSD/GU201905_Output/subset_bounds.json")
edna_data <- read_excel("/Volumes/GeringSSD/GU1905_eDNA/Copy_DNA_extracts_Qubit.xlsx", sheet = "R.Station.vs.DNA")

# Empty data frame to hold summary of MFI information - columns are Min., 1st Qu., Median, Mean, 3rd Qu., Max., Cast, Depth
summary_df <- data.frame(double(), double(), double(), double(), double(), double(), integer(), integer(), double())
colnames(summary_df) <- names(summary(0)) %>% append("Cast") %>% append("Depth") %>% append("ABC")

# loops through each cast in the json data
for (n_cast in seq_along(subset_data)){ # remove [1] - this is just to check functionality
  cast_name <- names(subset_data)[[n_cast]]
  print(paste("Cast:", cast_name))
  cast <- subset_data[[n_cast]]
  MFI_classes <- data.frame(Classification = character(), Count = double(), Depth = integer()) # MFI classifications for this cast

  # for each depth in the CTD cast
  for (n_depth in seq_along(cast)){
    depth_name <- names(cast)[[n_depth]]
    print(paste("Depth:", depth_name))
    
    # list of data arrays for each frequency at this depth
    depth_subsets <- cast[[n_depth]]
    depth_bounds <- subset_bounds[[n_cast]][[n_depth]][["depth"]]
    Sv38 <- depth_subsets[["38000"]]
    
    # Calculate MFI
    MFI = calc_MFI(depth_subsets, bad_fq=c("200000")) # do not use 200kHz in calculations
    
    if (!is.null(MFI)){
      # Calculate ABC for 38kHz - masked for fish
      mask_Sv <- mask_MFI(MFI, Sv38, list(c(0, 0.4), c(0.8, 1))) # fish with and without swimbladder ranges
      abc_val <- calc_ABC(mask_Sv, depth_bounds)
      
      # plot MFI (uncomment jpeg and dev.off to save)
      plot(t(MFI), breaks=c(-0.02, 0.4, 0.6, 0.8, 1.0), # some values are very slightly smaller that 0
           col=c_palette, 
           main= paste("CTD:", cast_name, "Depth: ", depth_name, "- MFI"))
      
      # histogram of MFI values
      hist(MFI, col = "blue", main = paste("CTD:", cast_name, "Depth: ", depth_name, "- Histogram"))
      
      # binning into 4 MFI classifications
      depth_class <- as.data.frame(table(classify_MFI(MFI))) %>% mutate(Depth = strtoi(depth_name))
      colnames(depth_class) <- c("Classification", "Count", "Depth")
      MFI_classes <- rbind(MFI_classes, depth_class) # save into MFI_classes data frame
      
      # add summary information to df dataframe
      MFI_summary <- summary(as.vector(MFI))
      temp_df <- data.frame(x=t(matrix(MFI_summary)))
      colnames(temp_df) <- names(MFI_summary)
      temp_df <- temp_df %>% mutate(Depth = strtoi(depth_name)) %>% mutate(Cast = strtoi(sub("^ctd0+", "", cast_name))) %>%
        mutate(ABC = abc_val)
      summary_df <- rbind(summary_df, temp_df)
    }
    
  }
  if (nrow(MFI_classes)> 0){
    
    # graph bar classifications
    par(mar=c(4,10,4,4))
    barplot(height = MFI_classes$Count, names = MFI_classes$Classification, density = MFI_classes$Depth, angle = MFI_classes$Depth,
            horiz = T, las=2, col=c_palette,
            main = paste("CTD:", cast_name, "- Classifications"))
    
    legend("bottomleft", legend = unique(MFI_classes$Depth), 
           density = unique(MFI_classes$Depth), angle = unique(MFI_classes$Depth), 
           bty=1:2, pt.cex = 1, cex = 1, horiz = FALSE, title = "Depths", bg='lightblue', inset = 0.055)
    
    abline_list = Map(function(x) x*4.9, seq(length(names(cast))-1))
    abline(h=abline_list, col="black", lwd = 2.5)
  }
}

edna_df <- read_eDNA(edna_data)
abc_df <- do.call(rbind, lapply(unique(edna_df$Cast), function(x) approx_depth_mutate(edna_df[edna_df$Cast == x, ], summary_df[summary_df$Cast == x,])))
df <- merge(edna_df, abc_df) # merge data frames by identical eDNA depth column

# Violin plot of log(ABC) values vs X12S Band Intensity for all casts and depths
abc_edna_violin(df)

# Scatter plot of log(ABC) values vs X12S Band Intensity for all casts and depths
abc_edna_scatter(df)

# Simple correlation testing
cor.test(as.numeric(df$X.12S.rank), log(df$ABC), method="kendall")
