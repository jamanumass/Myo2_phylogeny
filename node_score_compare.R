# This script will show the HMM search scores for the genes in a node (branch) relative to the genes outside that node. 
# It is useful for moving from a preliminary search + tree (which collects good hits and some worse ones) towards a search that excludes most genes outside your target clade while including most genuine hits.

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Set the run and tree results directory
run_name <- "results_myo1s_2"
node <- 223  # Example node to collect from (can be modified)

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

# Function to read scores from a hmmsearch result file
read_hmm_scores <- function(hmm_file, gene_ids_selected) {
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
  
  # Identify selected and non-selected gene scores
  selected_scores <- hmm_scores[hmm_gene_ids %in% gene_ids_selected]
  non_selected_scores <- hmm_scores[!hmm_gene_ids %in% gene_ids_selected]
  
  return(list(selected = selected_scores, non_selected = non_selected_scores))
}

# Main function to extract and compare HMM search scores
compare_hmm_scores <- function(tree_file, node, hmm_dir) {
  
  # Load the tree
  tree <- read.tree(tree_file)
  
  # Extract gene IDs from the selected node
  collected_from_node_raw <- extract.clade(tree, node)$tip.label
  gene_ids_selected <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", collected_from_node_raw)
  
  # Debug: Print the selected gene IDs
  message("Selected node gene IDs: ", paste(gene_ids_selected, collapse = ", "))
  
  # Get the directories that start with 'search_results_'
  subdirs <- list.dirs(hmm_dir, recursive = TRUE, full.names = TRUE)
  search_dirs <- subdirs[grepl("search_results_", basename(subdirs))]
  
  # Initialize lists to store all selected and non-selected scores
  all_selected_scores <- numeric(0)
  all_non_selected_scores <- numeric(0)
  
  # Loop through each directory
  for (search_dir in search_dirs) {
    
    # Find the hmmsearch result files in each directory
    hmm_files <- list.files(search_dir, pattern = "hmmsearch_out.txt", full.names = TRUE)
    
    # Loop through each hmmsearch result file
    for (hmm_file in hmm_files) {
      # Read and extract the scores from the hmmsearch file
      hmm_scores <- read_hmm_scores(hmm_file, gene_ids_selected)
      
      # Append the scores to the corresponding lists
      all_selected_scores <- c(all_selected_scores, hmm_scores$selected)
      all_non_selected_scores <- c(all_non_selected_scores, hmm_scores$non_selected)
    }
  }
  
  # Combine selected and non-selected scores into one data frame
  if (length(all_selected_scores) > 0 | length(all_non_selected_scores) > 0) {
    selected_df <- data.frame(score = all_selected_scores, category = "Selected")
    non_selected_df <- data.frame(score = all_non_selected_scores, category = "Non-selected")
    combined_df <- rbind(selected_df, non_selected_df)
    
    # Plot histograms for both selected and non-selected gene scores using ggplot2
    ggplot(combined_df, aes(x = score, fill = category)) +
      geom_histogram(binwidth = 10, position = "identity", alpha = 0.6, color = "black") +
      scale_fill_manual(values = c("Selected" = "blue", "Non-selected" = "red")) +
      labs(title = "Histogram of Selected and Non-selected Gene Scores", x = "Score", y = "Frequency") +
      theme_minimal()
  } else {
    message("No scores found for selected or non-selected node genes.")
  }
}

# Call the main function with the defined variables
compare_hmm_scores(tree_file, node, hmm_search_dir)
