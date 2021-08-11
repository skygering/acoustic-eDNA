# Script that reads in both eDNA CSV and ABC CSV (created by python code base) and merges the two by eDNA depth
# Script also creates two plots (violin and scatter plot) as well as calculating Kendall's Tau significance test
# Developer: Skylar Gering 08/2021
library(csv)
library(jsonlite)
library(readxl)
source("acoustic-eDNA/Rcode/classification_functions.R")

# Read in eDNA CSV (has 27 columns: Cruise, SampleID, Station, Station.Number, Lat, Long, Cast, Depth	Date, Time, Sampling.Depth
# Sampling.Depth.Meter, Sampling.Depth.Type, Gear, Filtration.Volume, Filt.Time, DNA.Conc, X.12S.Yuan.reconciled, X.12S.MAK-reconciled
# X.12S.final, X.18S.PCR, X.18S.T.MAK, Agreement.18S, X.18S.GE.MAK, X.16S.PCR, X.16S.MAK, Agreement.16S - may need to adjust script
# if your file has different columns) and ABC JSON file created using mfi_masks_abc.py in acoustic-eDNA.
edna_data <- read_excel("/Volumes/GeringSSD/GU1905_eDNA/Copy_DNA_extracts_Qubit.xlsx", sheet = "R.Station.vs.DNA")
abc_data <- jsonlite::fromJSON("/Volumes/GeringSSD/GU201905_Output/abc.json")

# Select specific columns from eDNA CSV and create new column called X.12s.rank which gives a number rank from 0-4 for eDNA 12s band intensity
edna_df <- read_eDNA(edna_data)

# Unpack all casts ABC data into one data frame
abc_df <- do.call(rbind, lapply(seq_along(abc_data), function(x) abc_df(abc_data, x)))

# Rank eDNA by quantile (note that this is just for given samples)
abc_df$abc_rank <- with(abc_df, cut(ABC, quantile(ABC, seq(0,1,0.2)), labels = 0:4, include.lowest = TRUE))

# For each cast in ABC data frame, add matching eDNA depth for each sample
abc_df <- do.call(rbind, lapply(unique(edna_df$Cast), function(x) approx_depth_mutate(edna_df[edna_df$Cast == x, ], abc_df[abc_df$Cast == x,])))
df <- merge(edna_df, abc_df) # merge data frames by identical eDNA depth column

df <- df[with(df, order(df$Cast, df$eDNA.Depth)),] # sort data frame by cast and then depth within casts

write.csv(df, "/Volumes/GeringSSD/GU201905_Output/edna_abc.csv")

# Violin plot of log(ABC) values vs X12S Band Intensity for all casts and depths
abc_edna_violin(df)

# Scatter plot of log(ABC) values vs X12S Band Intensity for all casts and depths
#abc_edna_scatter(df)

# Simple correlation testing
cor.test(as.numeric(df$X.12S.rank), log(df$ABC), method="kendall")

