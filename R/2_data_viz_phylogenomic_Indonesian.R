library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggrepel)
# source("global/fun.R")

# reload all designated df & tree (renamed)
tre_gubbins <- ape::read.tree("inputs/microreact_id_adjusted_gpsc14_.final_tree.tre")

combine_gpsc14 <- read.csv("inputs/microreact_id_adjusted_all_gpsc14_data.csv") %>% 
  glimpse()
combine_gpsc14_simplified <- read.csv("inputs/microreact_id_adjusted_all_gpsc14_data_simplified.csv") %>% 
  glimpse()

# test node
ggtree(tre_gubbins) + 
  geom_tiplab(size = 2) +
  geom_label2(aes(subset=!isTip, label=node),
              size=2, color="darkred", alpha=0.5)


# using collapse tree instead
tree <- ggtree(tre_gubbins) + 
  geom_tiplab(size = 2) #+
  # geom_label2(aes(subset=!isTip, label=node),
  #             size=2, color="darkred", alpha=0.5)

collapsed_tree <- tree %>% 
  collapse(node = 447) %>% 
  collapse(node = 628) %>% 
  collapse(node = 676) %>% 
  collapse(node = 863)

# label in 606 & 672
node_anno <- tibble(
  node = c(606, 672),
  labels = c("Clade I", "Clade II")
)

pointed_tree <- collapsed_tree %<+% 
  node_anno +
  geom_point2(aes(subset = (node == 447)),
              shape = 19,
              size = 5,
              color = "black") +
  geom_point2(aes(subset = (node == 628)),
              shape = 19,
              size = 5,
              color = "black") +
  geom_point2(aes(subset = (node == 676)),
              shape = 19,
              size = 5,
              color = "black") +
  geom_point2(aes(subset = (node == 863)),
              shape = 19,
              size = 5,
              color = "black") +
  # label in 606 & 672
  geom_hilight(node=606, fill="pink", alpha=0.5) +
  geom_hilight(node=672, fill="pink", alpha=0.5) +
  geom_text(
    aes(x = branch,
        label = labels),
    size  = 3,
    color = "darkred",
    vjust = -0.5,
    hjust = 0.8,
    fontface  = "bold",
    na.rm = TRUE
  )
pointed_tree


################################################################################
show_collapsed <- pointed_tree %<+% 
  combine_gpsc14 +
  # geom_tiplab(size = 2) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  theme(
    legend.position = "none",
    # plot.margin = grid::unit(c(-15, -15, -15, -15), "mm")
  )
show_collapsed

# gen tree #####################################################################
tree_gen_gubbins <- show_collapsed %<+%
  combine_gpsc14 +
  # serotype
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$In_silico_serotype),
    width=10,
    offset=0.12
    # axis.params = list(
    #   axis = "x",
    #   text = "Serotype",
    #   text.angle = -90,
    #   hjust = 1,
    #   vjust = 1,
    #   text.size = 2
    # )
  ) +
  scale_fill_viridis_d(
    name = "Serotype",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 2, order = 1)
  ) +
  # CC
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$CC),
    width=10
    # offset=0.02,
    
    # axis.params = list(
    #   axis = "x",
    #   text = "CC",
    #   text.angle = -90,
    #   hjust = 1,
    #   vjust = 1,
    #   text.size = 2
    # )
  ) +
  scale_fill_viridis_d(
    name = "Clonal complex",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 2, order = 2)
  ) +
  # Vaccination status
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$Vaccine_period),
    width=10
    # offset=0.02,
    # axis.params = list(
    #   axis = "x",
    #   text = "Vaccination",
    #   text.angle = -90,
    #   hjust = 1,
    #   vjust = 1,
    #   text.size = 2
    # )
  ) +
  scale_fill_viridis_d(
    name = "Vaccination status",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 2, order = 3)
  ) +
  # Continent
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$Continent),
    width=10
    # offset=0.02,
    # axis.params = list(
    #   axis = "x",
    #   text = "Continent",
    #   text.angle = -90,
    #   hjust = 1,
    #   vjust = 1,
    #   text.size = 2
    # )
  ) +
  scale_fill_viridis_d(
    name = "Continent",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 2, order = 4)
  ) +
  # Country
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=combine_gpsc14$Country),
    width=7
    # offset=0.02,
    # axis.params = list(
    #   axis = "x",
    #   text = "Country",
    #   text.angle = -90,
    #   hjust = 1,
    #   vjust = 1,
    #   text.size = 2
    # )
  ) +
  scale_fill_viridis_d(
    name = "Country",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 2, order = 5)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm"),
    legend.position = "right"
  )

# adjust coordinates for legends
built <- ggplot_build(tree_gen_gubbins)
for(i in seq_along(built$data)){
  d <- built$data[[i]]
  if(all(c("x", "xmin", "xmax") %in% names(d))){
    cat(paste(round(unique(d$x), 2), collapse = ", "), "\n")
  }
}

built$layout$panel_params[[1]]$y.range
x_range <- built$layout$panel_params[[1]]$x.range
x_range
y_max <- built$layout$panel_params[[1]]$y.range[2]
y_max

png("pictures/gpsc14_chosen_tree_IDN.png",
    width = 30, height = 20, units = "cm", res = 800)
tree_gen_gubbins <- tree_gen_gubbins +
  coord_cartesian(clip = "off",
                  xlim = c(x_range[1], x_range[2]),
                  ylim = c(0, y_max + 5)
                  ) +
  annotate("text",
           x = c(192.86, 198.02, 203.19, (208.36+1), (213.52+3)),
           y = y_max + 0.5,
           label = c("Serotype", "CC", "Vaccination", "Continent", "Country"),
           angle = 90,
           hjust = 0,
           vjust = -1,
           size = 3)
tree_gen_gubbins
dev.off()

