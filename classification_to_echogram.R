source("acoustic-eDNA/classification_functions.R")

json_data <- jsonlite::fromJSON("/Volumes/GeringSSD/Subset_Sv/segment_subset_CTD14-15.json")

ctd14_56 <- json_data$ctd014$`56`

MFI <- ctd14_56 %>% calc_MFI()

sb_fish <- apply(MFI, c(1,2), function(x) {if (between(x, -Inf, 0.4))  x  else 0})



sb_fish_1800 <- ctd14_56$`18000` * sb_fish

