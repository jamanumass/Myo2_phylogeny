# Improved Myosin Phylogeny Visualization Script
# =============================================

# Load required libraries with error handling
required_packages <- c("ape", "ggtree", "dplyr", "ggplot2", "treeio", "gridExtra", "phytools")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed. Please install it using install.packages('", pkg, "')"))
  }
}

#------------------------------------------------------------
# Define parameters and settings
#------------------------------------------------------------
BOOTSTRAP_THRESHOLD <- 75  # Threshold for collapsing nodes
WORK_DIR <- "/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results/sup_figure_run_16/"
OUTPUT_FILE <- "myosin_phylogeny.pdf"
WIDTH <- 12  # PDF width in inches
HEIGHT <- 16  # PDF height in inches

#------------------------------------------------------------
# Set working directory and load tree data
#------------------------------------------------------------
setwd(WORK_DIR)

# Load the first .contree file found
tree_path <- list.files(getwd(), pattern = "\\.contree$", full.names = TRUE)[1]

if (is.na(tree_path)) stop("No .contree file found.")

message("Using tree file: ", tree_path)

#------------------------------------------------------------
# Load and prepare tree
#------------------------------------------------------------
# Read the tree
tree <- tryCatch({
  read.tree(tree_path)
}, error = function(e) {
  stop("Error reading tree file: ", e$message)
})

# Define outgroup for rooting
outgroup_clade_name <- "MyoI" 
outgroup_clade_tips <- c("Acanthamoeba_castellanii_K17_Myo1A", "Naegleria_gruberi_K17_Myo1A")

# Check if outgroup tips exist in the tree
if (!all(outgroup_clade_tips %in% tree$tip.label)) {
  missing_tips <- outgroup_clade_tips[!outgroup_clade_tips %in% tree$tip.label]
  stop("Outgroup tips not found in tree: ", paste(missing_tips, collapse = ", "))
}

# Root the tree
root_node <- getMRCA(tree, outgroup_clade_tips)
if (is.null(root_node)) {
  stop("Could not find MRCA for outgroup tips.")
}
rooted_tree <- root(tree, node = root_node)

