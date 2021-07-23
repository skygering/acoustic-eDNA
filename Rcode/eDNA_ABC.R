library(jsonlite)
library(tidyr)
library(plyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(viridis)
library(zoo)

isclose <- function(a, b, abs_tol){
  return(abs(a-b) <= abs_tol)
}

abc_df <- function(data, cast){
  df <- do.call(rbind, data[[cast]])
  if(!is.null(dim(df))){
    df <- cbind(rownames(df), data.frame(df, row.names=NULL)) %>% mutate(cast)
    colnames(df) <- list("Depth", "ABC", "Cast")
  }
  return(df)
}

approx_merge <- function(edna, abc){
  return (mutate(abc, eDNA.Depth = unlist(lapply(abc$Depth, function(x) edna[isclose(strtoi(x), edna$eDNA.Depth, 5), ]$eDNA.Depth))))
}


edna_data <- read_excel("/Volumes/GeringSSD/GU1905_eDNA/Copy_DNA_extracts_Qubit.xlsx", sheet = "R.Station.vs.DNA")
abc_data <- jsonlite::fromJSON("/Volumes/GeringSSD/GU201905_Output/abc.json")

edna_data %>% 
  select(Cast, Lat, Long, Filtration.Volume, Sampling.Depth.Meter, Sampling.Depth.Type, X.12S.final) %>%
  mutate(X.12S.rank = case_when(
    endsWith(X.12S.final, "No") ~ 0,
    endsWith(X.12S.final, "EL") ~ 1,
    endsWith(X.12S.final, "L") ~ 2,
    endsWith(X.12S.final, "M") ~ 3, 
    endsWith(X.12S.final, "H") ~4)) %>% na.omit() %>%
  transform(X.12S.final=factor(X.12S.final,levels=c("No", "EL", "L", "M", "H"))) %>%
  subset(Sampling.Depth.Meter > 5 & Filtration.Volume == 2) -> edna_df


edna_df <- edna_df %>% select(-Filtration.Volume, -Sampling.Depth.Type)
colnames(edna_df)[which(names(edna_df) == "Sampling.Depth.Meter")] <- "eDNA.Depth"

abc_df <- do.call(rbind, lapply(seq_along(abc_data), function(x) abc_df(abc_data, x)))

abc_df$abc_rank <- with(abc_df, cut(ABC, quantile(ABC, seq(0,1,0.2)), labels = 0:4, include.lowest = TRUE))

abc_df <- do.call(rbind, lapply(unique(edna_df$Cast), function(x) approx_merge(edna_df[edna_df$Cast == x, ], abc_df[abc_df$Cast == x,])))
df <- merge(edna_df, abc_df)


df[with(df, order(df$Cast, df$eDNA.Depth)),]

ggplot(df, aes(X.12S.final, as.numeric(levels(abc_rank))[abc_rank], fill = X.12S.final)) + 
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  geom_jitter(width = 0.005) + ggtitle("Fish Area Backscattering Coefficent Quantile Rank vs 12S PCR eDNA Rank") +
  xlab("12S PCR eDNA Rank") + ylab("Fish ABC Quantile Rank") + theme(legend.position = "none") 
