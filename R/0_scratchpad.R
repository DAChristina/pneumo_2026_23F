# test node
ggtree(tre_gubbins) + 
  geom_tiplab(size = 2) +
  geom_label2(aes(subset=!isTip, label=node),
              size=2, color="darkred", alpha=0.5)

collapsed_tree <- tre_gubbins %>% 
  collapse(node = 447) %>% 
  collapse(node = 628) %>% 
  collapse(node = 676) %>% 
  collapse(node = 863)

# analyse subtree:
subtree1 <- ape::extract.clade(tre_gubbins,
                               node = 650)
ggtree(subtree1) + 
  geom_tiplab(size = 2) +
  geom_label2(aes(subset=!isTip, label=node),
              size=2, color="darkred", alpha=0.5)

df_subtree_collapsed <- dplyr::bind_rows(
  data.frame(node = 447,
             selected_id  = ape::extract.clade(tre_gubbins, 447)$tip.label),
  data.frame(node = 628,
             selected_id  = ape::extract.clade(tre_gubbins, 628)$tip.label),
  data.frame(node = 676,
             selected_id  = ape::extract.clade(tre_gubbins, 676)$tip.label),
  data.frame(node = 863,
             selected_id  = ape::extract.clade(tre_gubbins, 863)$tip.label),
) %>% 
  dplyr::left_join(
    dplyr::left_join(
      combine_gpsc14
      ,
      combine_gpsc14_simplified %>% 
        dplyr::rename_all(~ paste0("simplified_", .))
      ,
      by = c("tre_gubbins.tip.label" = "simplified_id")
    )
    ,
    by = c("selected_id" = "tre_gubbins.tip.label")
  ) %>%
  # view() %>%
  glimpse()
