library(tidyverse)
library(ggtree)
library(ggtreeExtra)
# source("global/fun.R")

combine_gpsc14 <- read.csv("inputs/all_gpsc14_data.csv") %>% 
  dplyr::mutate(id = gsub("\\.contigs", "",
                          gsub("#", "_", id)),
                # simplify other columns
                Vaccine_period = case_when(
                  stringr::str_detect(Vaccine_period, "Pre")  ~ "Pre-PCV",
                  stringr::str_detect(Vaccine_period, "Post") ~ "Post-PCV",
                  TRUE ~ "Unknown"
                  ),
                # fix MLST, CC & Age group
                In_silico_ST = case_when(
                  stringr::str_detect(In_silico_ST, "-") |
                    stringr::str_detect(In_silico_ST, stringr::fixed("*")) ~ "Not recognisable",
                  TRUE ~ In_silico_ST
                ),
                CC = case_when(
                  is.na(CC) |
                    stringr::str_detect(CC, "-") |
                    stringr::str_detect(CC, stringr::fixed("*")) ~ "Not recognisable",
                  TRUE ~ as.character(CC)
                  ),
                Age_group = ifelse(is.na(Age_group), "Unknown", Age_group)
                
                ) %>% 
  glimpse()

tre_gubbins <- ape::read.tree("outputs/gubbins/gpsc14_.final_tree.tre")
tre_gubbins$tip.label <- gsub("\\.contigs", "", gsub("#", "_", tre_gubbins$tip.label))

# focused on raxml tree: rearrange label coz' ggtree link label to row_names
combine_gpsc14 <- 
  dplyr::left_join(
    data.frame(tre_gubbins$tip.label)
    ,
    combine_gpsc14,
    by = c("tre_gubbins.tip.label" = "id")
  )

rownames(combine_gpsc14) <- tre_gubbins$tip.label

# compile data for microreact upload

write.csv(combine_gpsc14, "inputs/microreact_id_adjusted_all_gpsc14_data.csv")

# simplified csv
prop_serotype <- combine_gpsc14 %>% 
  dplyr::group_by(In_silico_serotype) %>% 
  dplyr::summarise(n_serotype = n(), .groups = "drop") %>% 
  dplyr::mutate(
    percent_serotype = round(n_serotype/nrow(combine_gpsc14)*100, 1)
  ) %>% 
  dplyr::arrange(desc(percent_serotype)) %>% 
  dplyr::mutate(
    # filter_serotype = ifelse(percent_serotype > 0.8, 1, 0),
    In_silico_serotype_simplified = ifelse(percent_serotype > 1,
                                           In_silico_serotype,
                                           "Other"),
  ) %>% 
  glimpse()

prop_CC <- combine_gpsc14 %>% 
  dplyr::group_by(CC) %>% 
  dplyr::summarise(n_CC = n(), .groups = "drop") %>% 
  dplyr::mutate(
    percent_CC = round(n_CC/nrow(combine_gpsc14)*100, 1)
  ) %>% 
  dplyr::arrange(desc(percent_CC)) %>% 
  dplyr::mutate(
    # filter_CC = ifelse(percent_CC > 1, 1, 0),
    CC_simplified = CC #ifelse(percent_CC > 1,
                                     #CC,
                                     #"Other"),
  ) %>% 
  glimpse()

prop_countries <- combine_gpsc14 %>% 
  dplyr::group_by(Country) %>% 
  dplyr::summarise(n_country = n(), .groups = "drop") %>% 
  dplyr::mutate(
    percent_country = round(n_country/nrow(combine_gpsc14)*100, 1)
  ) %>% 
  dplyr::arrange(desc(percent_country)) %>% 
  dplyr::mutate(
    # filter_CC = ifelse(percent_CC > 1, 1, 0),
    Country_simplified = ifelse(percent_country > 1.5, # include Malawi
                                Country,
                                "Other"),
  ) %>%
  glimpse()


combine_gpsc14_simplified <- combine_gpsc14 %>% 
  dplyr::left_join(
    prop_serotype %>% 
      dplyr::select(contains("In_silico_"))
    ,
    by = "In_silico_serotype"
  ) %>% 
  dplyr::left_join(
    prop_CC %>% 
      dplyr::select(c(CC, CC_simplified))
    ,
    by = "CC"
  ) %>% 
  dplyr::left_join(
    prop_countries %>% 
      dplyr::select(contains("Country"))
    ,
    by = "Country"
  ) %>% 
  # define rownames again
  dplyr::right_join(
    data.frame(tre_gubbins$tip.label)
    ,
    by = "tre_gubbins.tip.label"
  ) %>% 
  dplyr::transmute(
    id = tre_gubbins.tip.label,
    Continent = Continent,
    
    # adjust Indonesia on top
    Country = ifelse(Country_simplified == "Indonesia", " Indonesia",
                     Country_simplified),
    Serotype = In_silico_serotype_simplified,
    CC = CC_simplified,
    Vaccine = case_when(
      stringr::str_detect(Vaccine_period, "Pre")  ~ "Pre-PCV",
      stringr::str_detect(Vaccine_period, "Post") ~ "Post-PCV",
      TRUE ~ "Unknown"
    )
  ) %>% 
  glimpse()

rownames(combine_gpsc14_simplified) <- tre_gubbins$tip.label

write.csv(combine_gpsc14_simplified %>% 
            dplyr::select(-id),
          "inputs/microreact_id_adjusted_all_gpsc14_data_simplified.csv")

write.table(combine_gpsc14_simplified,
            "inputs/microreact_id_adjusted_all_gpsc14_data_simplified.tsv",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

ape::write.tree(tre_gubbins, "inputs/microreact_id_adjusted_gpsc14_.final_tree.tre")


# test node
ggtree(tre_gubbins) + 
  geom_tiplab(size = 2) +
  geom_label2(aes(subset=!isTip, label=node),
              size=2, color="darkred", alpha=0.5)

# analyse subtree:
subtree1 <- ape::extract.clade(tre_gubbins,
                              node = 446)

df_subtree <- dplyr::bind_rows(
  # 719 ST2059-dominant (South Africa)
  data.frame(node = 719,
             selected_id  = ape::extract.clade(tre_gubbins, 719)$tip.label),
  # 844 sg6 (South Africa)
  data.frame(node = 844,
             selected_id  = ape::extract.clade(tre_gubbins, 844)$tip.label),
  # 605 (mix country)/ split to 606 clade (IDN) & 628 (mixed countries)
  data.frame(node = 606,
             selected_id  = ape::extract.clade(tre_gubbins, 606)$tip.label),
  data.frame(node = 628,
             selected_id  = ape::extract.clade(tre_gubbins, 628)$tip.label),
  # 839 test IDN-Cambodia clade
  data.frame(node = 839,
             selected_id  = ape::extract.clade(tre_gubbins, 839)$tip.label),
  # 863 (peruvian clade)
  data.frame(node = 863,
             selected_id  = ape::extract.clade(tre_gubbins, 863)$tip.label),
  # 448 ST6279-dominant (South Africa)
  data.frame(node = 448,
             selected_id  = ape::extract.clade(tre_gubbins, 448)$tip.label),
  # 566 ST242-dominant (China/other country)
  data.frame(node = 566,
             selected_id  = ape::extract.clade(tre_gubbins, 566)$tip.label),
  
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

# proportion calculation per node
prop_calc1 <- df_subtree %>% 
  dplyr::group_by(node,
                  Country,
                  simplified_Serotype,
                  # simplified_CC
  ) %>% 
  dplyr::summarise(n_all = n(), .groups = "drop") %>% 
  dplyr::left_join(
    df_subtree %>% 
      dplyr::group_by(node) %>% 
      dplyr::summarise(n_node = n(), .groups = "drop")
    ,
    by = "node"
  ) %>% 
  dplyr::mutate(
    percent = round(n_all/n_node*100, 2),
    report = paste0(n_all, "/", n_node, " (", percent, "%)")
  ) %>% 
  arrange(desc(percent)) %>% 
  # view() %>% 
  glimpse()

prop_calc2 <- df_subtree %>% 
  dplyr::group_by(node,
                  Country,
                  # simplified_Serotype,
                  simplified_CC
  ) %>% 
  dplyr::summarise(n_all = n(), .groups = "drop") %>% 
  dplyr::left_join(
    df_subtree %>% 
      dplyr::group_by(node) %>% 
      dplyr::summarise(n_node = n(), .groups = "drop")
    ,
    by = "node"
  ) %>% 
  dplyr::mutate(
    percent = round(n_all/n_node*100, 2),
    report = paste0(n_all, "/", n_node, " (", percent, "%)")
  ) %>% 
  arrange(desc(percent)) %>% 
  # view() %>% 
  glimpse()


show_gubbins <- ggtree(tre_gubbins,
                  layout = "rectangular",
                  open.angle=30,
                  size=0.25,
                  # aes(colour=Clade)
) %<+% 
  combine_gpsc14 +
  # geom_tiplab(size = 2) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # # 23F GPSC14
  # geom_hilight(node=638, fill="#A63603", alpha=0.5) +
  # geom_cladelab(
  #   data = data.frame(
  #     node = 638,
  #     name = "23F\n(GPSC14-ST242)\ndominant"
  #   ),
  #   mapping = aes(
  #     node = node,
  #     label = name
  #   ),
  #   align = TRUE,
  #   offset = .23,
  #   offset.text = .045,
  #   hjust = "center",
  #   barsize = .2,
  #   fontsize = 3,
  #   angle = "auto",
  #   horizontal = FALSE
  # ) +
  theme(
    legend.position = "none",
    # plot.margin = grid::unit(c(-15, -15, -15, -15), "mm")
  )
show_gubbins

# gen tree #####################################################################
tree_gen_gubbins <- show_gubbins %<+%
  combine_gpsc14 +
  # serotype
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$In_silico_serotype),
    width=10,
    offset=0.05
  ) +
  scale_fill_viridis_d(
    name = "Serotype",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 7, order = 1)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # CC
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$CC),
    width=10
    # offset=0.02
  ) +
  scale_fill_viridis_d(
    name = "Clonal complex",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 7, order = 2)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # Vaccination status
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$Vaccine_period),
    width=10
    # offset=0.02
  ) +
  scale_fill_viridis_d(
    name = "Vaccination status",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 7, order = 3)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # Continent
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$Continent),
    width=10
    # offset=0.02
  ) +
  scale_fill_viridis_d(
    name = "Continent",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 7, order = 4)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # Country
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$Country),
    width=10
    # offset=0.02
  ) +
  scale_fill_viridis_d(
    name = "Country",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 5, order = 5)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +

  # ADDITIONAL INFORMATION
  # Manifestation
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$Manifestation),
    width=10
    # offset=0.02
  ) +
  scale_fill_viridis_d(
    name = "Manifestation",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 7, order = 6)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm"),
    legend.position = "right"
  ) #+

  # # Age groups
  # ggnewscale::new_scale_fill() +
  # ggtreeExtra::geom_fruit(
  #   geom=geom_tile,
  #   mapping=aes(fill=combine_gpsc14$Age_group),
  #   width=10
  #   # offset=0.02
  # ) +
  # scale_fill_viridis_d(
  #   name = "Age group",
  #   option = "C",
  #   direction = -1,
  #   guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
  #                        ncol = 7, order = 7)
  # ) +
  # theme(
  #   legend.title=element_text(size=12),
  #   legend.text=element_text(size=9),
  #   legend.spacing.y = unit(0.02, "cm"),
  #   legend.position = "right"
  # )

# png("pictures/phylo_gubbins_1epiTree.png",
#     width = 30, height = 20, units = "cm", res = 800)
tree_gen_gubbins
# dev.off()
