
treePlotter <- function(tree_file, outgroupID = NULL, label_size = .5 ) {
  # Read the Newick tree
  tree <- read.tree(tree_file)
  
  # Set the outgroup if provided
  if (!is.null(outgroupID)) {
    node_to_root <- which(grepl(outgroupID, tree$tip.label))
    tree <- root(tree, node_to_root, resolve.root = TRUE)
  }
  
  # Plot the tree using ggtree
  p_tree <- ggtree(tree, layout = "rectangular", size = .1) + 
    geom_tiplab(aes(label = label), 
                align = FALSE ,linetype = "dotted", linesize = 0.2, size = label_size, offset = 0.0) +  # Increase label size and adjust horizontal justification
    #theme_tree2() +   # Apply a tree theme
    geom_text2(aes(label = node), hjust = -0.1, vjust = 0.6, size = label_size * 0.7, color = "blue") +  # Node number labels
    geom_text2(aes(label = ifelse(isTip, NA, label)), hjust = -0.1, vjust = -0.6, size = label_size * 0.7, color = "red")  # Bootstrap support
   #xlim(NA, 3) # Adjust plot margins, increase right margin
  
  # Print the tree plot
  print(p_tree)
}

# Function to extract gene IDs from a tree node and save them to a specified file
save_gene_IDs_from_node <- function(tree_file, node_to_collect, output_file_IDs, format = "txt") {
  
  # Validate tree file
  if (!file.exists(tree_file)) {
    stop("Error: Tree file does not exist.")
  }
  
  # Load the tree from the specified file
  tree <- read.tree(tree_file)
  
  # Extract the tip labels (gene IDs) from the specified node
  collected_from_node_raw <- extract.clade(tree, node_to_collect)$tip.label
  
  # Process the labels to extract only the gene IDs
  collected_from_node_IDs <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", collected_from_node_raw)
  
  # Get the number of gene IDs collected
  num_gene_IDs_collected <- length(collected_from_node_IDs)
  message("Number of genes collected: ", num_gene_IDs_collected)
  
  # Save the collected gene IDs in the specified format
  tryCatch({
    if (format == "txt") {
      write.table(collected_from_node_IDs, file = output_file_IDs, row.names = FALSE, col.names = FALSE, quote = FALSE)
    } else if (format == "csv") {
      write.csv(collected_from_node_IDs, file = output_file_IDs, row.names = FALSE, quote = FALSE)
    } else {
      stop("Unsupported file format. Please use 'txt' or 'csv'.")
    }
  }, error = function(e) {
    stop("Error: Could not save the file. Ensure the directory is writable.")
  })
  
  # Final message after saving
  message("Gene IDs successfully saved to: ", output_file_IDs)
}

# Function to count the number of sequences in a FASTA file
count_sequences <- function(fasta_file) {
  # Read the entire file
  fasta_lines <- readLines(fasta_file)
  
  # Count the number of lines that start with ">"
  num_sequences <- sum(startsWith(fasta_lines, ">"))
  
  return(num_sequences)
}


extract_seqs_from_tree_node <- function(tree_file, fasta_file, node_to_collect, output_fasta_file) {
  
  # Validate input files
  if (!file.exists(tree_file)) stop("Error: Tree file does not exist.")
  if (!file.exists(fasta_file)) stop("Error: FASTA file does not exist.")
  
  # Load tree and extract gene IDs from the specified node
  tree <- read.tree(tree_file)
  collected_from_node_IDs <- extract.clade(tree, node_to_collect)$tip.label
  
  # Clean FASTA file using external script
  temp_fasta_file <- tempfile(fileext = ".fasta")
  cleaning_script <- "/home/jaman_umass_edu/jaman/general_tools/clean_gene_names.sh"
  if (system(paste(cleaning_script, fasta_file, temp_fasta_file)) != 0) {
    stop("Error: Cleaning script failed.")
  }
  
  # Read cleaned FASTA file and extract matching sequences
  fasta_lines <- readLines(temp_fasta_file)
  sequences <- list()
  current_seq <- NULL
  capture_sequence <- FALSE
  
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      if (!is.null(current_seq)) sequences <- c(sequences, list(current_seq))
      gene_id <- substring(line, 2)
      capture_sequence <- gene_id %in% collected_from_node_IDs
      current_seq <- if (capture_sequence) list(header = line, sequence = "") else NULL
    } else if (capture_sequence) {
      current_seq$sequence <- paste0(current_seq$sequence, line)
    }
  }
  if (!is.null(current_seq)) sequences <- c(sequences, list(current_seq))
  
  # Stop if no sequences were collected
  if (length(sequences) == 0) stop("Error: No sequences found for the provided gene IDs.")
  
  # Write sequences to output FASTA file
  con <- file(output_fasta_file, "w")
  for (seq in sequences) writeLines(c(seq$header, seq$sequence), con)
  close(con)
  file.remove(temp_fasta_file)
  
  # Report
  message("Report:")
  message("- Number of gene IDs in node: ", length(collected_from_node_IDs))
  message("- Number of sequences collected: ", length(sequences))
  message("- Output file: ", output_fasta_file)
}