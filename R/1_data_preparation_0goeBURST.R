library(stringr)
library(dplyr)
library(readr)

x <- readLines("raw_data/phyloviz_goeBURST_output.txt")
# annoying R dimension requirement
cc <- c()
st <- c()
rdimension_req <- NA
for(i in seq_along(x)){
  
  if(stringr::str_detect(x[i], "^CC")){
    rdimension_req <- stringr::str_extract(x[i], "\\d+")
  }
  
  if(stringr::str_detect(x[i], "^ST")){
    st <- c(st,
            stringr::str_extract(x[i], "(?<=ST )\\d+"))
    cc <- c(cc,rdimension_req)
  }
}

goeresult <- tibble(
  ST = st,
  CC = cc
)

write.csv(goeresult, "raw_data/phyloviz_goeBURST_output_table.csv",
          row.names = FALSE)
