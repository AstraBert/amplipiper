# load necessary R libraries
library("ggtree")
library("tidyverse")
library(phytools) # to determine the maximum tree height and add midpoint root
library(treeio)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]
title <- args[3]
offset <- as.numeric(args[4])
width <- as.numeric(args[5])
height <- as.numeric(args[6])
outgroup <- args[7]
watermark <- args[8]

## read tree
tree <- read.tree(input)

## root, if possible
if (outgroup != "no") {
  outgroup_samples <- unlist(str_split(outgroup, ","))

  # Check if all outgroup samples are in the tree's leaves (tip labels)
  if (all(outgroup_samples %in% tree$tip.label)) {
    tree <- root(tree, outgroup = outgroup_samples)
  } else {
    cat("Outgroup samples not found in tree's tip labels. No rooting performed.\n")
  }
}
## caluculate tree height (on x-axis)
tree$edge.length[tree$edge.length == "NaN"] <- 0
Xmax <- max(nodeHeights(tree))

## only retain Bootstrapping Support > 75%
tree$node.label[as.numeric(tree$node.label) < 0.75] <- NA

## plot tree
if (watermark == "YES") {
  PLOT.tree <- ggtree(tree,
    layout = "roundrect"
  ) +
    # ggtitle(title) +
    theme_tree2() +
    ggtitle("ASTRAL") +
    theme_bw() +
    ggplot2::xlim(
      0,
      Xmax + Xmax * offset
    ) +
    xlab("coalescent units") +
    geom_nodelab(
      hjust = 1.25,
      vjust = -0.75,
      size = 3,
      color = "#21614a"
    ) +
    annotate(
      geom = "text",
      x = (Xmax + Xmax * 0.2) / 2,
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
    # options(ragg.max_dim = 100000) +
    geom_tiplab()
} else {
  PLOT.tree <- ggtree(tree,
    layout = "roundrect"
  ) +
    # ggtitle(title) +
    theme_tree2() +
    theme_bw() +
    ggplot2::xlim(
      0,
      Xmax + Xmax * 0.2
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
    # options(ragg.max_dim = 100000) +
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
