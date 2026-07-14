# test full join for GPS-PMEN15 file

test <- read.csv("raw_data/monocle-metadata_GPSC14.csv") %>% 
  dplyr::full_join(
    read.csv("raw_data/pmen15.info.csv")
    ,
    by = c("Lane_id" = "ID")
  ) %>% 
  dplyr::filter(is.na(Published)) %>% 
  glimpse()


