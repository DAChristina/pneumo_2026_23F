# things to do:
# 1. report NT back to Harry


# 2. compile all 23F

# 3. construct script for running gubbins as submitted job
test <- combine_gpsc14 %>% 
  group_by(In_silico_ST, Country) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  dplyr::left_join(
    combine_gpsc14 %>% 
      group_by(Country) %>% 
      summarise(n_country = n()) %>% 
      ungroup()
    ,
    by = "Country"
  ) %>% 
  mutate(
    prop = n/n_country*100
  ) %>% 
  filter(prop > 10) %>% 
  view()

length(unique(test$In_silico_ST))

length(unique(combine_gpsc14$In_silico_ST))


# age stratification due to weird Age_years column 