# Load required libraries
library(tidyverse)

## read from commandline
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]

## read input file
DATA <- read.csv(input,
    header = TRUE,
    na.string = "NA"
)

## make heatmap
PLOT <- ggplot(DATA, aes(y = ID, x = Locus, fill = as.factor(NoAlleles))) +
    geom_tile(color = "black") +
    theme_bw() +
    ylab("Samples") +
    xlab("Loci") +
    scale_fill_manual(name = "Number of \nhaplotypes", values = c("#7ac1dc", "#ff8e3e", "#21614a", "#C3E6DA")) +
    # guides(fill = guide_legend(title = "Number of \nhaplotypes")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## make output path
OUTPUT <- paste0(input, ".png")

## save file as PNG
ggsave(PLOT,
    file = OUTPUT,
    width = length(levels(as.factor(DATA$Locus))),
    height = length(levels(as.factor(DATA$ID))) / 4
)

## make heatmap
PLOT <- ggplot(DATA, aes(y = ID, x = Locus, fill = as.factor(ExpectedPloidy))) +
    geom_tile(color = "black") +
    theme_bw() +
    ylab("Samples") +
    xlab("Loci") +
    scale_fill_manual(name = "Modelled\nPloidy", values = c("#7ac1dc", "#ff8e3e", "#21614a", "#C3E6DA")) +
    # guides(fill = guide_legend(title = "Modelled Ploidy")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## make output path
OUTPUT <- paste0(input, "_ExpectedPloidy.png")

## save file as PNG
ggsave(PLOT,
    file = OUTPUT,
    width = length(levels(as.factor(DATA$Locus))),
    height = length(levels(as.factor(DATA$ID))) / 4
)
