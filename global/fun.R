pct_summary <- function(data, group_var, fill_var, top_n_fill){
  fill_sym  <- rlang::sym(fill_var)
  group_sym <- rlang::sym(group_var)
  
  # rare = "Other"
  data <- data %>%
    dplyr::mutate(!!fill_sym := forcats::fct_lump_n(
      factor(!!fill_sym),
      n = top_n_fill,
      other_level = "Other"))
  
  data %>%
    dplyr::count(!!group_sym, !!fill_sym) %>%
    dplyr::ungroup() %>% 
    dplyr::left_join(
      data %>%
        dplyr::count(!!group_sym) %>%
        dplyr::ungroup() %>% 
        dplyr::rename(n_group = n)
      ,
      by = group_var
    ) %>% 
    dplyr::mutate(pct = n/n_group*100,
                  report = paste0(n, "/", n_group,
                                  "(", round(pct, 1), "%)")
                  ) %>%
    dplyr::ungroup()
}
