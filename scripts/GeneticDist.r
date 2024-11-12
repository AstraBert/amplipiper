# Load required libraries
library(seqinr)
library(ape)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]

setwd(input)

# Function to read alignment file and plot distance matrix
read_and_plot <- function(file_path, max_dist) {
    # Read the alignment file
    alignment <- read.alignment(file_path, format = "fasta")

    # Calculate genetic distances
    distances <- dist.alignment(alignment)

    # Convert distances to a distance matrix
    distance_matrix <- as.matrix(distances)

    # Convert to upper triangle matrix
    distance_matrix[lower.tri(distance_matrix)] <- NA

    # Convert to data frame
    distance_df <- melt(distance_matrix, na.rm = TRUE)

    # Extract ID from file name
    file_name <- basename(file_path)
    id <- str_extract(file_name, ".*(?=_aln.fasta)")

    # Plot the distance matrix as a heatmap using ggplot2
    p <- ggplot(data = distance_df, aes(Var2, Var1, fill = value)) +
        geom_tile() +
        scale_fill_gradient(low = "#7ac1dc", high = "#ff8e3e", limits = c(0, max_dist)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = paste("Pairwise Genetic Distances -", id), x = "Sequence", y = "Sequence")

    return(p)
}

# Get a list of all files with the extension "_aln.fasta" in the current directory and its subdirectories
file_paths <- list.files(path = ".", pattern = "_aln.fasta$", recursive = TRUE, full.names = TRUE)

# Calculate the maximum distance across all files
max_dist <- 0
for (file_path in file_paths) {
    # Read the alignment file
    alignment <- read.alignment(file_path, format = "fasta")

    # Calculate genetic distances
    distances <- dist.alignment(alignment)

    # Find the maximum distance
    max_dist <- max(max_dist, max(distances))
}

# Calculate the optimal dimensions for the grid
num_files <- length(file_paths)
num_cols <- ceiling(sqrt(num_files))
num_rows <- ceiling(num_files / num_cols)

# Create a list to store plots
plots <- list()

# Read each alignment file and plot the distance matrix
for (file_path in file_paths) {
    p <- read_and_plot(file_path, max_dist)
    plots[[file_path]] <- p
}

# Save the plot as a PNG file
ggsave("distance_matrices.png", 
    arrangeGrob(grobs = plots, ncol = num_cols),
    width = num_cols * 7, 
    height = num_rows * 6
)
