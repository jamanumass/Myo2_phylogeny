treePlotterWithExonMaps <- function(tree_file, exon_data = NULL, outgroupID = NULL, label_size = 1, bar_height = 3, y_shift = 0) {
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
  
  
  # Extract tip heights for y-values
  tip_heights <- p_tree$data[p_tree$data$isTip, c("label", "y")]
  colnames(tip_heights) <- c("gene_id", "y_max")
  tip_heights$y_min <- tip_heights$y_max - bar_height
  
  # Match short gene IDs to long gene IDs in the tree
  long_gene_ids <- tree$tip.label
  
  # Create a map of short to long gene IDs
  gene_id_map <- sapply(names(exon_data), function(short_id) {
    matched_long_id <- grep(paste0("^", short_id), long_gene_ids, value = TRUE)
    if (length(matched_long_id) > 0) {
      return(matched_long_id[1])  # Use the first match if multiple
    } else {
      return(NA)
    }
  })
  
  # Update exon_data to use long gene IDs
  names(exon_data) <- gene_id_map
  
  # Filter out any NA matches to avoid issues later
  exon_data <- exon_data[!is.na(names(exon_data))]
  
  
  
    # Modify exon_data to include y_max and y_min from tip_heights
    for (gene_id in names(exon_data)) {
      # Find the corresponding y-value for the gene_id
      y_values <- tip_heights$y_max[match(gene_id, tip_heights$gene_id)]
      
      # Check if y_values is found and valid
      if (!is.na(y_values)) {
        exon_df <- exon_data[[gene_id]]
        
        # Add y_max and y_min to exon_df
        exon_df$y_max <- y_values + y_shift
        exon_df$y_min <- y_values + y_shift - bar_height
        
        exon_data[[gene_id]] <- exon_df
      }
    }
    
    
    for (gene_id in names(exon_data)) {
      print(gene_id)
    }
    
    
    
    
    
    
    
    # Create the exon plot facet
    p_exons <- ggplot() +
      geom_rect(data = do.call(rbind, exon_data),
                aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
                fill = "green", color = "black", size = 0.2) +
      theme_void()
    
    # Extract the y limits from the tree plot
    y_limits <- ggplot_build(p_tree)$layout$panel_params[[1]]$y.range
    
    # Apply these limits to the exon plot
    p_exons <- p_exons +
      coord_cartesian(ylim = y_limits)
    
    # Combine tree and exon plots
    combined_plot <- plot_grid(p_tree, p_exons, ncol = 2, rel_widths = c(2, 1)) # Adjust `rel_widths` as needed
    print(combined_plot)

}