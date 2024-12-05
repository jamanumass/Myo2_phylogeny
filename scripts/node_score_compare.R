# This script will show the HMM search scores for the genes in a node (branch) relative to the genes outside that node. 
# It is useful for moving from a preliminary search + tree (which collects good hits and some worse ones) towards a search that excludes most genes outside your target clade while including most genuine hits.

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ape)


# Set the run and tree results directory
run_name <- "results_Myo2_Pedro"
node <- 92  # Example node to collect from (can be modified)

# Define the universal variables
results_dir <- "/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results"
run_dir <- file.path(results_dir, run_name)
hmm_search_dir <- file.path(run_dir, "hmm_search_results")
tree_dir <- file.path(run_dir, "tree_inference_result")

# Find the tree file dynamically from the tree directory
tree_files <- list.files(path = tree_dir, pattern = "\\.contree$", full.names = TRUE)
if (length(tree_files) != 1) {
  stop("Error: There should be exactly one tree file in the directory.")
}
tree_file <- tree_files[1]  # Use the first (and only) tree file



# Function to read scores from a hmmsearch result file and map them to specific genes
read_hmm_scores <- function(hmm_file, gene_ids_selected, gene_ids_outside) {
  # Read the lines of the file
  lines <- readLines(hmm_file)
  
  # Skip lines that start with '#' (comments) and extract only non-comment lines
  data_lines <- lines[!grepl("^#", lines)]
  
  # If no data lines exist, return empty scores
  if (length(data_lines) == 0) {
    return(list(selected = numeric(0), non_selected = numeric(0)))
  }
  
  # Read the table from the processed lines, using `fill = TRUE` to handle rows with varying column numbers
  hmm_data <- read.table(text = paste(data_lines, collapse = "\n"), header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  
  # Columns for the full sequence score and gene ID
  col_gene_id <- 1  # First column: target name (gene ID)
  col_score <- 6    # Sixth column: full sequence score
  
  # Extract the gene IDs and scores
  hmm_gene_ids <- as.character(hmm_data[, col_gene_id])
  hmm_scores <- as.numeric(hmm_data[, col_score])
  
  # Map scores to their respective genes
  selected_scores <- hmm_scores[hmm_gene_ids %in% gene_ids_selected]
  non_selected_scores <- hmm_scores[hmm_gene_ids %in% gene_ids_outside]
  
  # Return lists of genes with their corresponding scores
  return(list(
    selected = selected_scores, 
    non_selected = non_selected_scores
  ))
}

# Function to calculate and print statistics
print_stats <- function(scores_list) {
  if (length(scores_list) > 0) {
    scores <- as.numeric(scores_list)
    message("Number of genes: ", length(scores))
    message("Min: ", min(scores, na.rm = TRUE))
    message("Max: ", max(scores, na.rm = TRUE))
    message("Median: ", median(scores, na.rm = TRUE))
    message("Mean: ", mean(scores, na.rm = TRUE))
  } else {
    message("No scores available.")
  }
}

# Main function to extract and compare HMM search scores and print only the stats
get_hmm_scores_for_node <- function(tree_file, node, hmm_dir) {
  
  # Load the tree and extract gene IDs from the entire tree and from the node
  tree <- read.tree(tree_file)
  all_gene_ids <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", tree$tip.label)
  genes_in_node <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", extract.clade(tree, node)$tip.label)
  genes_outside_node <- setdiff(all_gene_ids, genes_in_node)
  
  # Debug: Print the selected gene IDs for verification
  message("Selected node gene IDs: ", length(genes_in_node))
  message("Genes outside node: ", length(genes_outside_node))
  
  # Initialize lists to store all selected and non-selected scores
  all_selected_scores <- numeric(0)
  all_non_selected_scores <- numeric(0)
  
  # Get the directories that start with 'search_results_'
  subdirs <- list.dirs(hmm_dir, recursive = TRUE, full.names = TRUE)
  search_dirs <- subdirs[grepl("search_results_", basename(subdirs))]
  
  # Loop through each directory
  for (search_dir in search_dirs) {
    
    # Find the hmmsearch result files in each directory
    hmm_files <- list.files(search_dir, pattern = "hmmsearch_out.txt", full.names = TRUE)
    
    # Loop through each hmmsearch result file
    for (hmm_file in hmm_files) {
      # Read and extract the scores from the hmmsearch file
      hmm_scores <- read_hmm_scores(hmm_file, genes_in_node, genes_outside_node)
      
      # Append the scores to the corresponding lists
      all_selected_scores <- c(all_selected_scores, hmm_scores$selected)
      all_non_selected_scores <- c(all_non_selected_scores, hmm_scores$non_selected)
    }
  }
  
  # Print the stats for both groups
  message("\nStatistics for genes in the node:")
  print_stats(all_selected_scores)
  
  message("\nStatistics for genes outside the node:")
  print_stats(all_non_selected_scores)
}

# Call the function with the defined variables
get_hmm_scores_for_node(tree_file, node, hmm_search_dir)
