
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
                align = FALSE ,linetype = "dotted", linesize = 0.2, size = label_size, offset = 0.1) +  # Increase label size and adjust horizontal justification
    #theme_tree2() +   # Apply a tree theme
    geom_text2(aes(label = node), hjust = -01.6, vjust = 0.5, size = label_size * 0.7, color = "blue") + # Apply node number labels
    geom_text2(aes(label = ifelse(isTip, NA, label)), hjust = -0.4, vjust = 0.5, size = label_size * 0.7, color = "red")# + # Apply node bootstrap support
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
  num_genes_collected <- length(collected_from_node_IDs)
  message("Number of genes collected: ", num_genes_collected)
  
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




# Function to extract full sequences for gene IDs from a node and save them to a FASTA file
extract_sequences_from_node <- function(tree_file, node_to_collect, fasta_file, output_fasta_file) {
  
  # Validate input files
  if (!file.exists(tree_file)) {
    stop("Error: Tree file does not exist.")
  }
  if (!file.exists(fasta_file)) {
    stop("Error: FASTA file does not exist.")
  }
  
  # Load the tree from the specified file
  tree <- read.tree(tree_file)
  
  # Extract the tip labels (gene IDs) from the specified node
  collected_from_node_raw <- extract.clade(tree, node_to_collect)$tip.label
  
  # Process the labels to extract only the gene IDs
  collected_from_node_IDs <- gsub("([A-Z0-9]+\\.[0-9]+)_.*", "\\1", collected_from_node_raw)
  
  # Get the number of gene IDs collected
  num_genes_collected <- length(collected_from_node_IDs)
  message("Number of gene IDs collected: ", num_genes_collected)
  
  # Initialize a list to store the sequences
  sequences <- list()
  current_seq <- NULL
  capture_sequence <- FALSE
  
  # Read the FASTA file line by line
  fasta_lines <- readLines(fasta_file)
  
  # Loop through the FASTA file to find the sequences corresponding to the collected gene IDs
  for (line in fasta_lines) {
    
    # If the line is a header (starts with ">")
    if (startsWith(line, ">")) {
      
      # If a sequence was being captured, store the current sequence
      if (!is.null(current_seq)) {
        sequences <- c(sequences, list(current_seq))
      }
      
      # Extract the gene ID from the header
      gene_id <- sub("^>(\\S+).*", "\\1", line)
      
      # Check if the gene ID is in the list of collected gene IDs
      if (gene_id %in% collected_from_node_IDs) {
        # Start capturing the sequence for this gene ID
        capture_sequence <- TRUE
        current_seq <- list(header = line, sequence = "")
      } else {
        # Stop capturing sequence if not in the collected IDs
        capture_sequence <- FALSE
        current_seq <- NULL
      }
      
    } else if (capture_sequence) {
      # If capturing sequence, append this line to the current sequence
      current_seq$sequence <- paste0(current_seq$sequence, line)
    }
  }
  
  # Add the last sequence (if any)
  if (!is.null(current_seq)) {
    sequences <- c(sequences, list(current_seq))
  }
  
  # Check if any sequences were captured
  if (length(sequences) == 0) {
    stop("Error: No sequences were found for the provided gene IDs.")
  }
  
  # Write the sequences to the output FASTA file
  con <- file(output_fasta_file, "w")
  for (seq in sequences) {
    writeLines(c(seq$header, seq$sequence), con)
  }
  close(con)
  
  # Final message after saving
  message("Sequences successfully saved to: ", output_fasta_file)
}