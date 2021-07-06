library(jsonlite)
library(tidyr)
library(plyr)
library(dplyr)
library(readxl)
library(ggplot2)

ctd_mean_sd <- function(ctd, n_ctd){
  n_depth = nrow(ctd)
  n_fq = ncol(ctd)
  cast <- unlist(rep(n_ctd, times = n_depth * n_fq))
  fq <- unlist(rep(colnames(ctd), times = n_depth))
  depth <- unlist(rep(rownames(ctd), each = n_fq))
  avg <- unlist(lapply(t(ctd), function(x) log10(mean(10^(x/10)))*10))
  sdev <- unlist(lapply(t(ctd), sd)) # does this need to go thought the whole linear process?
  df <- data.frame(cast, fq, depth, avg, sdev)
  return(df)
}

ctd_points <- function(df, edna_data, json_data, n_ctd){
  print("Ctd:")
  print(n_ctd)
  
  edna_data %>% subset(Cast == n_ctd) -> edna_df
  print(unlist(edna_df$Sampling.Depth.Meter))
  
  do.call(rbind, json_data[[n_ctd]])  -> cast_data
  
  if (length(cast_data) != 0 & nrow(edna_df) != 0){
    ctd_df <- cast_data %>% ctd_mean_sd(n_ctd)
    ctd_df <- mutate(ctd_df, ranks = lapply(ctd_df$depth, function(x) unlist(edna_df[isclose(edna_df$Sampling.Depth.Meter, strtoi(x), 5), 
                                                                             grep("^X.12S.rank$", colnames(edna_df))]))) %>% unnest_wider(ranks)
    n_rank = c(6:ncol(ctd_df))
    
    df <- do.call("rbind", lapply(n_rank, function(x) {if(length(n_rank) > 1) ctd_df[,-setdiff(n_rank, x)] else ctd_df} %>% rename_at(6, ~"rank")))
  }
  
  return(df)
}

isclose <- function(a, b, abs_tol){
  return(abs(a-b) <= abs_tol)
}

n_ctd = 1:24

df <- data.frame(cast = double(), fq = character(), 
                 depth = character(), avg = double(), 
                 sdev = double(), rank =double())

edna_data <- read_excel("/Volumes/GeringSSD/Copy_DNA_extracts_Qubit.xlsx", sheet = "R.Station.vs.DNA")
json_data <- jsonlite::fromJSON("/Volumes/GeringSSD/GU201905_CTD_JSON/segment_subset.json")

edna_data %>% 
  select(Cast, Lat, Long, Sampling.Depth.Meter, Sampling.Depth.Type, X.12S.final) %>%
  mutate(X.12S.rank = case_when(
    endsWith(X.12S.final, "No") ~ 0,
    endsWith(X.12S.final, "EL") ~ 1,
    endsWith(X.12S.final, "L") ~ 2,
    endsWith(X.12S.final, "M") ~ 3, 
    endsWith(X.12S.final, "H") ~4)) %>% na.omit() %>% 
  subset(Sampling.Depth.Meter > 5) -> edna_data

df <- do.call("rbind", lapply(n_ctd, function(x) rbind(df, ctd_points(df, edna_data, json_data, x))))
  
ggplot(df[df$fq != "200000",], aes(x = rank, y = avg)) + geom_point(aes(shape = fq, color = cast)) + stat_smooth()

