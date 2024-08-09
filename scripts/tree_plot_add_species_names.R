
# Function to process a FASTA file and create a table with full names
process_fasta_to_table <- function(fasta_file) {
  fasta_lines <- readLines(fasta_file)
  names <- c()
  sequences <- c()
  current_sequence <- ""
  
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      if (current_sequence != "") {
        sequences <- c(sequences, current_sequence)
        current_sequence <- ""
      }
      
      processed_name <- gsub(" ", "_", substring(line, 2))  # Replace spaces with underscores
      processed_name <- gsub("\\[|\\]", "", processed_name)  # Remove square brackets
      names <- c(names, processed_name)
    } else {
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  if (current_sequence != "") {
    sequences <- c(sequences, current_sequence)
  }
  
  fasta_table <- data.frame(Name = names, Sequence = sequences, stringsAsFactors = FALSE)
  return(fasta_table)
}

# Function to replace gene IDs in a Newick tree with full names from the fasta_table
replace_gene_ids_in_tree <- function(tree_file, fasta_table, output_file) {
  # Read the Newick tree
  tree <- read.tree(tree_file)
  
  # Initialize a vector to store the updated tip labels
  updated_labels <- character(length(tree$tip.label))
  
  # Loop through each gene ID in the tree's tip labels
  for (i in seq_along(tree$tip.label)) {
    gene_id <- tree$tip.label[i]
    
    # Attempt to match the gene ID exactly to the beginning of fasta_table names
    matched_name <- fasta_table$Name[grep(paste0("^", gene_id), fasta_table$Name)]
    
    if (length(matched_name) == 1) {
      # If a match is found, use the full name from the fasta_table
      updated_labels[i] <- matched_name
    } else {
      # If no match is found, keep the original gene ID
      updated_labels[i] <- gene_id
    }
  }
  
  # Update the tree's tip labels with the full names
  tree$tip.label <- updated_labels
  
  # Write the updated tree to a new file
  write.tree(tree, file = output_file)
  
  cat("Gene IDs replaced. Updated tree saved to:", output_file, "\n")
}

# Usage example
fasta_file <- "~/jaman/Myo2_phylogeny/results_Amoebozoan_MyoII_motors/tree_results/Amoebozoan_MyoII_motors_hits_aligned.fasta"
tree_file <- "~/jaman/Myo2_phylogeny/results_Amoebozoan_MyoII_motors/tree_results/Amoebozoan_MyoII_motors_hits_aligned.fasta.treefile"
output_file <- "~/jaman/Myo2_phylogeny/results_Amoebozoan_MyoII_motors/tree_results/Amoebozoan_MyoII_motors_hits_aligned.fasta.cleaned.treefile"

# Process the FASTA file into a table
fasta_table <- process_fasta_to_table(fasta_file)

# Replace gene IDs in the tree with the full names using the updated matching logic
replace_gene_ids_in_tree(tree_file, fasta_table, output_file)


# Function to plot the tree with long labels and extra space
plot_tree_with_long_labels <- function(tree_file, width_scale = 15) {
  # Read the Newick tree
  tree <- read.tree(tree_file)
  
  # Plot the tree using ggtree with extra space for labels
  p <- ggtree(tree) + 
    geom_tiplab(align = TRUE, linesize = 0.5) +  # Align labels and add connecting lines
    theme_tree2() +  # Apply a tree theme
    xlim_tree(width_scale)  # Adjust the width of the tree for better label visibility
  
  # Display the plot
  print(p)
}

# Usage example
output_file <- "~/jaman/Myo2_phylogeny/results_Amoebozoan_MyoII_motors/tree_results/Amoebozoan_MyoII_motors_hits_aligned.fasta.cleaned.treefile"

# Plot the updated tree with extra space for labels
plot_tree_with_long_labels(output_file, width_scale = 15)