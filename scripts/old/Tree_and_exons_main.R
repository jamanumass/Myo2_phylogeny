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

# Specify the tree file and outgroup
tree_file <- "/home/jaman_umass_edu/jaman/Myo2_phylogeny/results/results_Amoe_MyoII_motors/tree_inference_result/Amoe_MyoII_motors_input_clean_out.fas.contree" 
outgroupID <- "XP_641363.1_myosin_IA_heavy_chain_Dictyostelium_discoideum_AX4"

## Exon data - generate
# Extract gene IDs from the tree file in a single line
gene_ids <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", read.tree(tree_file_name)$tip.label)
# Pass the gene_ids to the exon mapper functions
exon_data <- map_exons_to_genes(gene_ids)

# Plot the tree with Exon maps
treePlotterWithExonMaps(tree_file, exon_data, outgroupID, label_size = 2, bar_height = .23, y_shift = 2)
treePlotter(tree_file, outgroupID, label_size = 2)

# Source the treePlotter function from the external script
source("~/jaman/Myo2_phylogeny/scripts/treePlotter.R")

outgroupID <- "Bnat_34821_MyoIV"

treePlotter(tree_file, outgroupID, label_size = 1)



# Extract gene IDs from a branch of the tree
node_to_collect <- 291 # This will be printed on the tree nodes from treePlotter
save_gene_IDs_from_node(tree_file, node_to_collect, format = "txt")

