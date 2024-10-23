library(ggtree)
library(ape)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)


# Source the tools
source("~/jaman/Myo2_phylogeny/scripts/treePlotterTools.R")

# Set the outgroup
outgroupID <- "Bnat_34821_MyoIV"

# Set the run and tree results directory
run_name <- "pedros_1f_my_2s"
results_dir <- "/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results"
run_dir <- file.path(results_dir, run_name)
hmm_search_dir <- file.path(run_dir, "hmm_search_results")
tree_dir <- file.path(run_dir, "tree_inference_result")
rendering_dir <- file.path(run_dir, "tree_render_results")

dir.create(rendering_dir)


tree_file <- list.files(path = tree_dir, pattern = "\\.contree$", full.names = TRUE)

treePlotter(tree_file, outgroupID, label_size = 1)

# Extract gene IDs from a node
node_to_collect <- 228 # This will be printed on the tree nodes from treePlotter
output_file_IDs <- file.path(rendering_dir, sprintf("%s.node%d.IDs.txt", run_name, node_to_collect))
save_gene_IDs_from_node(tree_file, node_to_collect, output_file_IDs)

# Extract the full sequences from a node
fasta_file <- list.files(path = hmm_search_dir, pattern = "\\_collected_seqs.fa$", full.names = TRUE) #should have all genes in the tree
output_fasta_file <- file.path(rendering_dir, sprintf("%s.node%d.seqs.fasta", run_name, node_to_collect))
extract_sequences_from_node(  tree_file, node_to_collect, fasta_file, output_fasta_file)

source("~/jaman/Myo2_phylogeny/scripts/treePlotterTools.R")

#count sequences to see if you got the expected amount
count_sequences(output_fasta_file)

