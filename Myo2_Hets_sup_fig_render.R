# Load Required Libraries
library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)

############################################################
# Load and Process Tree
############################################################
# Define the directory and load the tree data
tree_dir <- "/home/jaman_umass_edu/jaman/Myo2_phylogeny/results/results_Myo2_Pedro_8/tree_inference_result"
contree_file_path <- list.files(tree_dir, pattern = "\\.contree$", full.names = TRUE)
tree <- read.tree(contree_file_path)  # Read the tree
tree <- root(tree, outgroup = tree$tip.label[Descendants(tree, getMRCA(tree, c("Tgon_XP_002368942", "Tgon_XP_002366595")), type = "tips")[[1]]], resolve.root = TRUE) #root tree on myoI

# Cap branch lengths to a maximum value
max_branch_length <- 1  # Set your desired maximum branch length
tree$edge.length <- pmin(tree$edge.length, max_branch_length)

# Convert nodes with low bootstrap support to polytomies
bootstrap_threshold <- 20
tree_poly <- as.polytomy(tree, feature = "node.label",
                         fun = function(x) as.numeric(x) < bootstrap_threshold)

# Prepare bootstrap data for annotation
n_tips <- length(tree_poly$tip.label)
n_nodes <- tree_poly$Nnode
bootstrap_values <- as.numeric(tree_poly$node.label)
bootstrap_data <- data.frame(node = (n_tips + 1):(n_tips + n_nodes),
                             bootstrap = bootstrap_values)

############################################################
# Define Groups and Assign Colors to Tips
############################################################
# Define species groups
species_groups <- list(
  Metazoa = c("Hsap", "Dpul", "Dmel", "Ctel", "Nvec", "Mlei", "Ocar", "Aque", "Homo", "Trichoplax",
              "Caenorhabditis", "Gallus", "Danio", "Drosophila", "Ciona_intestinalis",
              "Caligus_rogercressey", "Amphimedon", "Branchiostoma", "Calanus", "Strongylocentrotus", "Oscarella", "Mnemiopsis"),
  `Unicellular Holozoa` = c("Clim_", "Mbre", "Sros", "Cowc", "Mvib", "Sarc", "Cfra",
                            "Awhi", "Pgem", "Apar", "Capsaspora", "Salpingoeca",
                            "Fonticula", "Sphaeroforma", "Parvularia_atlantis", "Codosiga", "Corallochytrium", "Diaphanoeca", "Monosiga", "Ichthyophonus", "Amoebidium", "Creolimax"),
  Apusozoa = c("Ttra", "Thecamonas"),
  Amoebozoa = c("Ddis", "Ppal", "Ehis", "Acas", "Tieghemostelium", "Polysphondylium",
                "Dictyostelium", "Cavenderia", "Heterostelium", "Pelomyxa", "Entamoeba",
                "Planoprotostelium", "Acanthamoeba", "Stygamoeba", "Vermistella",
                "Vexilliferidae", "Paramoeba", "Vannella", "Cochliopodium",
                "Amoeba_proteus", "Protostelium", "Filamoeba", "Mastigamoeba",
                "Arcella_intermedia", "Vermamoeba_vermiformis", "Dracoamoeba_jomungandri"),
  Archaeplastida = c("Atha", "Acoe", "Sbic", "Smoe", "Ppat", "Crei", "Vcar", "Esil", "Ngad",
              "Ptri", "Tpse", "Madagascaria"),
  SAR = c("Tgon", "Pfal", "Ptet", "Tthe", "Pmar", "Bnat", "Pinf", "Colponemidia"),
  Metamonada = c("Tvag", "Glam", "Anaeramoeba"),
  Discoba = c("Ngru", "Tcru", "Lmaj", "Naegleria", "Willaertia", "Neovahlkampfia",
              "Acrasis", "Pharyngomonas", "Naegr", "Tetramitus", "fowleri"),
  Fungi = c("Aory", "Spom", "Tmel", "Cneo", "Ccin", "Umay", "Mver", "Pbla", "Rory", "Rhizopus",
            "Amac", "Bden", "Spun", "Cryptococcus_neoformans", "Ustilago", "Saccharomyces", "Aspergillus", "Neurospora", "Yarrowia", "Rozella", "Paraphelidium_tribonemae", "Allomyces", "Rhizophagus")
)

# Assign colors to each group
group_colors <- c(
  Metazoa = "#FBBB86",
  `Unicellular Holozoa` = "orange",
  Fungi = "#D86325",
  Apusozoa = "#615D44",
  Amoebozoa = "#00A79C",
  Archaeplastida = "#8CC540",
  SAR = "#8266AD",
  Metamonada = "#F2798E",
  Discoba = "#9A2268",
  Other = "black"
)

# Function to assign group and color based on tip label
assign_group_and_color <- function(label) {
  for (group in names(species_groups)) {
    patterns <- species_groups[[group]]
    for (pattern in patterns) {
      if (grepl(pattern, label, ignore.case = TRUE)) {
        return(c(group = group, color = group_colors[[group]]))
      }
    }
  }
  return(c(group = "Other", color = group_colors[["Other"]]))
}

# Create a data frame with tip labels and their assigned groups and colors
tip_data <- data.frame(label = tree_poly$tip.label, stringsAsFactors = FALSE)
group_and_color <- do.call(rbind, lapply(tip_data$label, assign_group_and_color))
tip_data <- cbind(tip_data, group_and_color)
tip_data$group <- factor(tip_data$group, levels = names(group_colors))

############################################################
# Define Clades and Retrieve Their MRCA Nodes
############################################################
# Define clade specifications with tips and labels
clade_specs <- list(
  list(tips = c("Hsap_NP_510880", "Cowc_CAOG_02007_12"), label = "MyoXVIII"),
  list(tips = c("Hsap_NP_004990", "Cowc_CAOG_02007_4"), label = "MyoVI"),
  list(tips = c("Tgon_XP_002364981", "Tgon_XP_002367823"), label = "MyoXXIII"),
  list(tips = c("Atha_AT2G20290", "Crei_XP_001693409"), label = "MyoXI"),
  list(tips = c("Acas_1_29", "Ddis_XP_645195"), label = "MyoXXXIII"),
  list(tips = c("Hsap_NP_000250", "Bden_Batde5_13697"), label = "MyoV"),
  list(tips = c("Tpse_XP_002297047", "Tgon_XP_002370590"), label = "MyoXXVII"),
  list(tips = c("Tpse_XP_002297174", "Tpse_XP_002287576"), label = "MyoXXI"),
  list(tips = c("Atha_AT3G19960", "Crei_XP_001696752"), label = "MyoVIII"),
  list(tips = c("Hsap_NP_004989", "Ddis_XP_636359"), label = "MyoI"),
  list(tips = c("Hsap_NP_001077084", "Hsap_NP_059129"), label = "MyoIII"),
  list(tips = c("Hsap_AAI46792", "Cowc_CAOG_02007_8"), label = "MyoXVI"),
  list(tips = c("Hsap_AAD49195", "Cowc_CAOG_02007_7"), label = "MyoIX"),
  list(tips = c("Hsap_NP_036466", "Cowc_CAOG_02007_10"), label = "MyoX"),
  list(tips = c("EP00067_Danio_rerio_P025987", "XP_020436862_1_myosin_II_heavy_chain_Heterostelium_album_PN500"), label = "MyoII")
  
)

# Retrieve MRCA nodes for each clade
clade_nodes <- lapply(clade_specs, function(cspec) {
  getMRCA(phy = tree_poly, tip = cspec$tips)
})

# Combine MRCA nodes and labels into a data frame for annotation
clade_annotation_data <- do.call(rbind, lapply(seq_along(clade_specs), function(i) {
  data.frame(node = clade_nodes[[i]], label = clade_specs[[i]]$label)
}))

############################################################
# Plot and Annotate the Tree
############################################################
# Create the base tree plot
base_plot <- ggtree(tree_poly, layout = "rectangular", size = 0.1)

# Merge bootstrap data with tree data for annotations
tree_ggdata <- base_plot$data
tree_ggdata <- merge(tree_ggdata, bootstrap_data, by = "node", all.x = TRUE)

# Define label size for tip labels and node labels
label_size <- 1.5

# Build the annotated tree plot
annotated_plot <- base_plot %<+% tip_data +
  # Add tip labels with colors
  geom_tiplab(aes(label = label, color = color),
              linesize = 0.2, size = label_size,
              offset = -0.5, align = TRUE, show.legend = FALSE) +
  
  # Adjust the x-axis limits to accommodate labels and clades
  xlim(0, max(tree_ggdata$x, na.rm = TRUE) * 1.5) +
  
  # Apply theme and add title
  theme_tree2(base_size = 14) +
  ggtitle("Myo2s and other Myosins from Sebé-Pedrós 2014") +
  
  # Add bootstrap support circles at internal nodes with sufficient support
  geom_point(data = subset(tree_ggdata, !isTip & !is.na(bootstrap) & bootstrap >= bootstrap_threshold),
             aes(x = x, y = y, fill = bootstrap),
             shape = 21, color = "black", size = 1.5) +
  scale_fill_gradientn(colors = c("white", "gray", "black"),
                       limits = c(bootstrap_threshold, 100),
                       name = "Bootstrap") +
  
  # Add internal node numbers
  geom_text2(data = subset(tree_ggdata, !isTip),
             aes(x = x + 0.1, y = y, label = node),
             size = label_size * 0.5, color = "blue") +
  
  # Add clade labels based on computed MRCA nodes
  geom_cladelab(data = clade_annotation_data,
                mapping = aes(node = node, label = label),
                offset = -0.25, fontsize = 2, align = FALSE, barsize = 0.3) +
  
  # Use identity scale for colors to apply the specified group colors
  scale_color_identity()

############################################################
# Add a Legend for Groups
############################################################
# Prepare legend data excluding the 'Other' group
legend_data <- data.frame(
  group = names(group_colors),
  color = unname(group_colors),
  stringsAsFactors = FALSE
)
legend_data <- legend_data[legend_data$group != "Other", ]

# Define legend box position and size
legend_x <- 0.01
legend_y <- 250
legend_width <- 1.5
legend_height <- 50
n_entries <- nrow(legend_data)
entry_height <- legend_height / n_entries

# Calculate positions for each legend entry
legend_data <- legend_data %>%
  mutate(x = legend_x,
         y = legend_y - (seq_len(n_entries) - 0.5) * entry_height)

# Add the legend to the plot
final_plot <- annotated_plot +
  theme(legend.position = c(0.05, 0.5),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = NA, size = 0 )) +

  # Add legend text with corresponding colors
  geom_text(data = legend_data,
            aes(x = x, y = y, label = group, color = color),
            hjust = 0, size = 3, fontface = "bold", show.legend = FALSE) +
  
  # Ensure colors are correctly mapped
  scale_color_identity()

# Display the final annotated tree plot
print(final_plot)
