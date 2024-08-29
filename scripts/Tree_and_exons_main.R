library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)


# Source the exon_mapper.R script to load the necessary functions
source("~/jaman/Myo2_phylogeny/scripts/exon_mapper.R")

# Source the treePlotterWithExonMaps.R script to load the tree plotting functions
source("~/jaman/Myo2_phylogeny/scripts/treePlotterWithExonMaps.R")

# Path to the Newick tree file
tree_file <- "~/jaman/Myo2_phylogeny/results_MyoI_motors/tree_results/MyoI_motors_hits_aligned.fasta.contree"

# Extract gene IDs from the tree file in a single line
gene_ids <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", read.tree(tree_file_name)$tip.label)

# Pass the gene_ids to the exon mapper functions
exon_data <- map_exons_to_genes(gene_ids)



# Plot the tree
treePlotterWithExonMaps(tree_file, exon_data, outgroupID, label_size = 2, bar_height = .23, y_shift = 2)
treePlotter(tree_file, outgroupID, label_size = 2)


# Extract gene IDs from a branch of the tree
tree_file <- "~/path_to_tree/treefile.newick" # Specify the tree file and its full path
node_to_collect <- 132 # This will be printed on the tree nodes from treePlotter
save_gene_IDs_from_node(tree_file, node_to_collect, format = "txt")