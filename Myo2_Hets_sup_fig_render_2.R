library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)
library(treeio)
library(gridExtra)
library(phytools)

#------------------------------------------------------------
# Define the directory and load the tree data
#------------------------------------------------------------
setwd("/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results/sup_figure_run_16/")
tree_dir <- getwd()
contree_file_path <- list.files(tree_dir, pattern = "\\.contree$", full.names = TRUE)
treefile_path <- list.files(tree_dir, pattern = "\\.treefile$", full.names = TRUE)

tree <- read.tree(contree_file_path)  # Read the tree, decide if coming from .treefile or .contree
# view an unrooted tree to pick the outgroup
#ggtree(tree, layout="equal_angle") + geom_tiplab()


#select root node and produce rooted tree
outgroup_clade_name <- "MyoI" 
outgroup_clade_tips <- c("Acanthamoeba_castellanii_K17_Myo1A", "Naegleria_gruberi_K17_Myo1E")
root_node <- getMRCA(tree, outgroup_clade_tips)
rooted_tree <- root(tree, node = root_node)

# Identify nodes below threshold and collapse them
bootstrap_threshold <- 75


# Use `as.polytomy()` to collapse weak nodes
polytomy_tree <- as.polytomy(tree = rooted_tree, feature = "node.label", fun = function(x) as.numeric(x) < bootstrap_threshold)


# Plot results
#ggtree(rooted_tree) + geom_tiplab() + ggtitle("Original Tree")
#ggtree(polytomy_tree) + ggtitle("Collapsed Tree with Polytomies")

#------------------------------------------------------------
# Define Species Groups and Colors
#------------------------------------------------------------
species_groups <- list(
  Opisthokonts = c("Paraphelidium", "Allomyces", "Rhizophagus", "Homo", "Drosophila", "Capsaspora", "Salpingoeca", "Codosiga", "Corallochytrium", "Diaphanoeca", "Monosiga", "Ichthyophonus", "Amoebidium", "Creolimax", "Batrachochytrium", "Cryptococcus", "Ustilago", "Saccharomyces", "Aspergillus", "Neurospora", "Yarrowia", "Rozella", "Fonticula", "Sphaeroforma"),
  Other_Amorphea = c("Thecamonas", "Rigifila", "Pygsuia"),
  Amoebozoa = c("Acanthamoeba", "Amoeba", "Arcella", "Dictyostelium", "Entamoeba", "Filamoeba", "Mastigamoeba", "Vermamoeba", "Dracoamoeba"),
  Archaeplastida = c("Arabidopsis", "Chlamydomonas", "Selaginella"),
  SAR = c("Toxoplasma", "Trypanosoma", "Bigelowiella", "Aurantiochytrium", "Plasmodiophora", "Nitzschia", "Phytophthora", "Colponemidia"),
  Discoba = c("Soginia", "Naegleria", "Leishmania", "Willaertia", "Neovahlkampfia", "Acrasis", "Pharyngomonas", "Tetramitus", "Vahlkampfia"),
  Other = c()  # This can stay empty or contain any default species name
)

# Group colors
group_colors <- c(
  Opisthokonts = "#007782",
  Other_Amorphea = "#615D44",
  Amoebozoa = "#00A79C",
  Archaeplastida = "#8CC540",
  SAR = "#8266AD",
  Discoba = "#9A2268",
  Other = "black"
)

# Function to assign group and color
assign_group <- function(label) {
  for (grp in names(species_groups)) {
    for (pattern in species_groups[[grp]]) {
      if (grepl(pattern, label, ignore.case = TRUE)) {
        return(c(group = grp, color = group_colors[[grp]]))
      }
    }
  }
  c(group = "Other", color = group_colors[["Other"]])
}

# Create tip data
tip_data <- data.frame(label = polytomy_tree$tip.label, stringsAsFactors = FALSE)

# Assign groups and colors to the tips
group_assignments <- do.call(rbind, lapply(tip_data$label, assign_group))
tip_data <- cbind(tip_data, group_assignments)

# Convert the group to a factor
tip_data$group <- factor(tip_data$group, levels = names(group_colors))




#------------------------------------------------------------
# Create the Base Plot with Colored Tip Labels
#------------------------------------------------------------
base_plot <- ggtree(polytomy_tree, size = 0.1) %<+% tip_data +
  geom_tiplab(aes(label = label, color = color), size = 1.5, align = FALSE, show.legend = FALSE) +
  scale_color_identity() +
  theme_tree2() 
#print(base_plot + ggtitle("base plot"))


#------------------------------------------------------------
# Add Clade Labels
#------------------------------------------------------------
clade_label <- function(label, tree, tips, offset = .4) {
  node <- getMRCA(tree, tips)
  geom_cladelab(node = node, label = label, fontsize = 3, align = F, extend = 0.5, barsize = 0.3, offset = offset) 
}

plot_with_clades <- base_plot +
  clade_label("Myo2", polytomy_tree, c("Aspergillus_nidulans_eukprot_EP00135_P004703", "Pharyngomonas_kirbyi_eukprot_EP00761_P011290")) +
  clade_label("Myo1", polytomy_tree, c("Homo_sapiens_K17_Myo1F", "Dictyostelium_discoideum_K17_Myo1A")) + 
  clade_label("Myo5", polytomy_tree, c("Capsaspora_owczarzaki_K17_Myo5", "Drosophila_melanogaster_K17_Myo5")) + 
  clade_label("Myo8", polytomy_tree, c("Selaginella_moellendorffii_K17_Myo8B", "Chlamydomonas_reinhardtii_K17_Myo8")) + 
  clade_label("Myo11", polytomy_tree, c("Arabidopsis_thaliana_K17_Myo11E", "Chlamydomonas_reinhardtii_K17_Myo11A"))  
  


#print(plot_with_clades + ggtitle("plot with clades"))

#------------------------------------------------------------
# Add Bootstrap values
#------------------------------------------------------------
# Extract bootstrap values and associate with node numbers
bootstrap_values <- as.numeric(polytomy_tree$node.label)
node_numbers <- (length(rooted_tree$tip.label) + 1):(length(polytomy_tree$tip.label) + polytomy_tree$Nnode)
bootstrap_data <- data.frame(node = node_numbers, bootstrap = bootstrap_values, stringsAsFactors = FALSE)

# Merge bootstrap data with the tree data
plot_with_clades$data <- left_join(plot_with_clades$data, bootstrap_data, by = "node")

#find bootstrap for specific nodes of intertest to display as numbers
nodes_of_interest <- c(
  bootstrap_data[bootstrap_data$node == getMRCA(polytomy_tree, c("Thecamonas_trahens_eukprot_EP00031_P006505", "Aspergillus_nidulans_eukprot_EP00135_P004703")),],
  bootstrap_data[bootstrap_data$node == getMRCA(polytomy_tree, c("Filamoeba_sp_ATCC50430_eukprot_EP00027_P007255", "Arcella_intermedia_eukprot_EP01124_P025971")),],
  bootstrap_data[bootstrap_data$node == getMRCA(polytomy_tree, c("Naegleria_fowleri_MhcC_NCBI_XP_044562059_1", "Pharyngomonas_kirbyi_eukprot_EP00761_P011290")),],
  bootstrap_data[bootstrap_data$node == getMRCA(polytomy_tree, c("Aspergillus_nidulans_eukprot_EP00135_P004703", "Pharyngomonas_kirbyi_eukprot_EP00761_P011290")),]
)


# Add circles to internal nodes colored by bootstrap support
plot_with_bootstrap <- plot_with_clades +
  geom_point(
    data = subset(plot_with_clades$data, !isTip & !is.na(bootstrap) & bootstrap > bootstrap_threshold),
    aes(x = x, y = y, fill = bootstrap),
    shape = 21, color = "black", size = 1  # Set size to 1
  ) +
  geom_text( #Bootstrap labels
    data = subset(plot_with_clades$data, node %in% nodes_of_interest),
    aes(x = x, y = y, label = bootstrap),
    nudge_x = 0.04,  
    size = 2,
    color = "red"
  ) +
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap", limits = c(bootstrap_threshold, 100)) +
  theme(legend.position = c(0.8, 0.75))  

#print(plot_with_bootstrap + ggtitle("plot with bootstrap"))


#------------------------------------------------------------
# Add Legend for Species Groups to Upper Left Corner
#------------------------------------------------------------
# Prepare legend data
legend_data <- data.frame(
  group = names(group_colors),
  color = unname(group_colors),
  stringsAsFactors = FALSE
)
legend_data <- legend_data[legend_data$group != "Other", ]

# Add legend text
plot_with_legend <- plot_with_bootstrap +
  xlab("Average substitutions per site") +
  geom_text(
    data = legend_data,
    aes(x = 0, y = Ntip(tree) - (seq_len(nrow(legend_data)) * 3), label = group, color = color),
    hjust = 0,  # Align text to the left
    vjust = 0.5,
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_identity() 

#print(plot_with_legend + ggtitle("plot with legend"))


#------------------------------------------------------------
# Assign to Final Plot and Display
#------------------------------------------------------------
final_plot <- plot_with_legend +
  theme(
    plot.margin = margin(15.5, 15, 5.5, 10)  # Adjust left margin to accommodate legend
  )
#  xlim(0, 8) +
print(final_plot + ggtitle("Eukaryotic Myosins") )
#ggsave(filename = "revised.pdf", plot=last_plot() )
