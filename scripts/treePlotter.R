treePlotter <- function(tree_file, outgroupID = NULL, label_size = 1) {
  # Read the Newick tree
  tree <- read.tree(tree_file)
  
  # Set the outgroup if provided
  if (!is.null(outgroupID)) {
    node_to_root <- which(grepl(outgroupID, tree$tip.label))
    tree <- root(tree, node_to_root, resolve.root = TRUE)
  }
  
  # Plot the tree using ggtree
  p_tree <- ggtree(tree, layout = "rectangular") + 
    geom_tiplab(aes(label = label), 
                align = TRUE, linesize = 0.5, size = label_size, offset = 1) +  # Increase label size and adjust horizontal justification
    theme_tree2() +   # Apply a tree theme
    geom_text2(aes(label = node), hjust = -0.3, vjust = 0.5, size = label_size * 0.7, color = "blue") + # Apply node labels
    geom_text2(aes(label = ifelse(isTip, NA, label)), hjust = -0.3, vjust = -0.5, size = label_size * 0.7, color = "red") + # Apply node bootstrap support
    xlim(NA, 14) # Adjust plot margins, increase right margin
  
  # Print the tree plot
  print(p_tree)
}

# Function to extract the gene IDs from a tree branch
save_gene_IDs_from_node <- function(tree_file, node_to_collect, format = "txt") {
  # Load the tree from the specified file
  tree <- read.tree(tree_file)
  
  # Extract the tip labels from the specified node
  collected_from_node_raw <- extract.clade(tree, node_to_collect)$tip.label
  
  # Process the labels to extract the gene IDs
  collected_from_node_IDs <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", collected_from_node_raw)
  
  # Get the number of gene IDs collected
  num_genes_collected <- length(collected_from_node_IDs)
  
  # Print the number of genes collected to the console
  message("Number of genes collected: ", num_genes_collected)
  
  # Get the directory of the tree file and prepare the output file name
  tree_file_name <- basename(tree_file)
  tree_file_name <- tools::file_path_sans_ext(tree_file_name)
  output_dir <- dirname(tree_file)
  output_file <- file.path(output_dir, paste0(tree_file_name, "_collected_gene_IDs_node_", node_to_collect, ".", format))
  
  # Save the collected gene IDs to the specified file
  if (format == "txt") {
    write.table(collected_from_node_IDs, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else if (format == "csv") {
    write.csv(collected_from_node_IDs, file = output_file, row.names = FALSE, quote = FALSE)
  } else {
    stop("Unsupported file format. Please use 'txt' or 'csv'.")
  }
  
  message("Gene IDs saved to: ", output_file)
}