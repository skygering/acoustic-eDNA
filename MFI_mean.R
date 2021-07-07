library(jsonlite)
library(tidyr)
library(plyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)

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
json_data <- jsonlite::fromJSON("/Volumes/GeringSSD/GU201905_CTD_JSON/segment_subset_1m.json")

edna_data %>% 
  select(Cast, Lat, Long, Filtration.Volume, Sampling.Depth.Meter, Sampling.Depth.Type, X.12S.final) %>%
  mutate(X.12S.rank = case_when(
    endsWith(X.12S.final, "No") ~ 0,
    endsWith(X.12S.final, "EL") ~ 1,
    endsWith(X.12S.final, "L") ~ 2,
    endsWith(X.12S.final, "M") ~ 3, 
    endsWith(X.12S.final, "H") ~4)) %>% na.omit() %>%
  transform(X.12S.final=factor(X.12S.final,levels=c("No", "EL", "L", "M", "H"))) %>%
  subset(Sampling.Depth.Meter > 5) -> edna_data

df <- do.call("rbind", lapply(n_ctd, function(x) rbind(df, ctd_points(df, edna_data, json_data, x))))

# all points with trend line
ggplot(df[df$fq != "200000",], aes(x = rank, y = avg))  + geom_point(aes(shape = fq, color = cast)) + stat_smooth()
# box plots for all points
ggplot(na.omit(df[df$fq != "200000",]), aes(x = factor(rank), y = avg, fill = rank))  + geom_boxplot()

cast1 <- df[df$cast ==1,]
# this looks awful ...
f<- ggplot(cast1, aes(x = rank, y = avg, ymin = avg - sdev, ymax = avg + sdev)) + geom_errorbar(width = 0.2) + geom_point(size = 1.5)

# by-hand cast14
do.call(rbind, json_data[[14]])  %>% ctd_mean_sd(14) -> cast14
avg_38k <- cast14[cast14$depth == 31 & cast14$fq == "38000", ]$avg
depth31 <- cast14[cast14$depth == 31, ] %>% mutate(fq_r = avg/avg_38k)
avg_38k <- cast14[cast14$depth == 56 & cast14$fq == "38000", ]$avg
depth56 <- cast14[cast14$depth == 56, ] %>% mutate(fq_r = avg/avg_38k)
cast14 <- rbind(depth31, depth56)
ggplot(cast14, aes(x = strtoi(fq), y = fq_r, color = depth)) + geom_line() + ggtitle("Cast 14")

do.call(rbind, json_data[[15]])  %>% ctd_mean_sd(15) -> cast15
avg_38k <- cast15[cast15$depth == 42 & cast15$fq == "38000", ]$avg
depth42 <- cast15[cast15$depth == 42, ] %>% mutate(fq_r = avg/avg_38k)
avg_38k <- cast15[cast15$depth == 22 & cast15$fq == "38000", ]$avg
depth22 <- cast15[cast15$depth == 22, ] %>% mutate(fq_r = avg/avg_38k)
avg_38k <- cast15[cast15$depth == 9 & cast15$fq == "38000", ]$avg
depth9 <- cast15[cast15$depth == 9, ] %>% mutate(fq_r = avg/avg_38k)
cast15 <- rbind(depth42, depth22, depth9)
ggplot(cast15, aes(x = strtoi(fq), y = fq_r, color = depth)) + geom_line() + ggtitle("Cast 15")


# Samples were data was taken with all 3 water bottle sizes
edna_all_sample <- do.call("rbind", lapply(unique(edna_data$Cast), 
                                    function(x) do.call("rbind", lapply(unique(edna_data[edna_data$Cast == x,]$Sampling.Depth.Meter),
                                    function(y) { if (identical(unique(edna_data[edna_data$Cast == x & edna_data$Sampling.Depth.Meter == y,]$Filtration.Volume), c(1,2,3)))
                                                    edna_data[edna_data$Cast == x & edna_data$Sampling.Depth.Meter == y,]
                                                  else edna_data[0,]}))))

# Bar Plots
plot_vol_bars <- function(data, volume, y_max){
  return (ggplot(data[data$Filtration.Volume == volume, ], aes(x = X.12S.final, fill=X.12S.final)) +  geom_bar() + scale_x_discrete(drop=F) +
                 scale_fill_brewer(palette = "Set2") + coord_cartesian(ylim=c(0,y_max)))
}
figure <- ggarrange(plot_vol_bars(edna_all_sample, 1, 18), plot_vol_bars(edna_all_sample, 2, 18), plot_vol_bars(edna_all_sample, 3, 18),
                    labels = c("1L", "2L", "3L"),
                    ncol = 3, nrow = 1, common.legend = TRUE, legend = "right", hjust = -0.1)
annotate_figure(figure, top = text_grob("Counts of X.12S for 3 Water Bottle Sizes"))

