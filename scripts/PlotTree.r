# load necessary R libraries
library(ggtree)
library(tidyverse)
library(phytools) # to determine the maximum tree height and add midpoint root

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]
title <- args[3]
offset <- as.numeric(args[4])
width <- as.numeric(args[5])
height <- as.numeric(args[6])
outgroup <- args[7]
watermark <- args[8]

## load tree file and root midpoint
tree <- read.tree(input)

if (outgroup == "no") {
  tree <- midpoint.root(tree)
} else {
  outgroup_samples <- unlist(str_split(outgroup, ","))

  # Check if all outgroup samples are in the tree's leaves (tip labels)
  if (all(outgroup_samples %in% tree$tip.label)) {
    tree <- root(tree, outgroup = outgroup_samples)
  } else {
    cat("Outgroup samples not found in tree's tip labels. Midpoint rooting will be applied.\n")
    tree <- midpoint.root(tree)
  }
}
## caluculate tree height (on x-axis)
Xmax <- max(nodeHeights(tree))

## only retain Bootstrapping Support > 75%
tree$node.label[as.numeric(tree$node.label) < 75] <- NA

if (watermark == "YES") {
  ## plot tree
  PLOT.tree <- ggtree(tree,
    layout = "roundrect"
  ) +
    ggtitle(title) +
    theme_tree2() +
    theme_bw() +
    ggplot2::xlim(
      0,
      Xmax + Xmax * offset
    ) +
    xlab("av. subst./site") +
    geom_nodelab(
      hjust = 1.25,
      vjust = -0.75,
      size = 3,
      color = "#21614a"
    ) +
    annotate(
      geom = "text",
      x = (Xmax + Xmax * offset) / 2,
      y = 1 + length(tree$node.label) / 2,
      label = "PREVIEW",
      color = "#ff8e3e",
      angle = 45,
      fontface = "bold",
      size = 25,
      alpha = 0.5,
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_cartesian(clip = "off") +
    options(ragg.max_dim = 100000) +
    geom_tiplab()
} else {
  ## plot tree
  PLOT.tree <- ggtree(tree,
    layout = "roundrect"
  ) +
    ggtitle(title) +
    theme_tree2() +
    theme_bw() +
    ggplot2::xlim(
      0,
      Xmax + Xmax * offset
    ) +
    xlab("av. subst./site") +
    geom_nodelab(
      hjust = 1.25,
      vjust = -0.75,
      size = 3,
      color = "#21614a"
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_cartesian(clip = "off") +
    options(ragg.max_dim = 100000) +
    geom_tiplab()
}
PNG <- paste0(output, ".png")
PDF <- paste0(output, ".pdf")
## export tree
ggsave(
  filename = PDF,
  PLOT.tree,
  width = width,
  height = height, limitsize = FALSE
)
ggsave(
  filename = PNG,
  PLOT.tree,
  width = width,
  height = height, limitsize = FALSE
)
