library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)
library(treeio)
library(gridExtra)


#------------------------------------------------------------
# Define the directory and load the tree data
#------------------------------------------------------------
run <- c(1:22)
setwd("/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results/sup_figure_run_14/")
tree_dir <- getwd()
contree_file_path <- list.files(tree_dir, pattern = "\\.contree$", full.names = TRUE)
treefile_path <- list.files(tree_dir, pattern = "\\.treefile$", full.names = TRUE)

tree <- read.tree(contree_file_path)  # Read the tree

#select root node
outgroup_clade_name <- "MyoI" 
outgroup_clade_tips <- c("Homo_sapiens_SP14_NP_005370", "Homo_sapiens_SP14_NP_004989")

#.contree
contree_file_path <- list.files(tree_dir, pattern = "\\.contree$", full.names = TRUE)
tree <- read.tree(contree_file_path)  # Read the tree
root_node <- getMRCA(tree, outgroup_clade_tips)
rooted_tree <- root(tree, node = root_node)
p1 <- ggtree(rooted_tree) + geom_tiplab(aes(label = label), size = 1.5)

#.treefile
treefile_path <- list.files(tree_dir, pattern = "\\.treefile$", full.names = TRUE)
treefile <- read.tree(treefile_path)  # Read the tree
treefile_root_node <- getMRCA(treefile, outgroup_clade_tips)
rooted_treefile_tree <- root(treefile, node = treefile_root_node)
p2 <- ggtree(rooted_treefile_tree) + geom_tiplab(aes(label = label), size = 1.5)

grid.arrange(p1, p2, ncol = 2)  # Arrange in 2 columns

#------------------------------------------------------------
# Define Species Groups and Colors
#------------------------------------------------------------
species_groups <- list(
  Opisthokonts = c("Homo_sapiens", "Trichoplax", "Caenorhabditis", "Gallus", "Danio", "Drosophila", "Ciona_intestinalis", "Caligus_rogercressey", "Amphimedon", "Branchiostoma", "Calanus", "Strongylocentrotus", "Oscarella", "Mnemiopsis", "Capsaspora", "Salpingoeca", "Fonticula", "Sphaeroforma", "Parvularia_atlantis", "Codosiga", "Corallochytrium", "Diaphanoeca", "Monosiga", "Ichthyophonus", "Amoebidium", "Creolimax", "Rhizopus", "Batrachochytrium", "Cryptococcus_neoformans", "Ustilago", "Saccharomyces", "Aspergillus", "Neurospora", "Yarrowia", "Rozella", "Paraphelidium_tribonemae", "Allomyces", "Rhizophagus"),
  Apusozoa = c("Thecamonas"),
  Amoebozoa = c("Tieghemostelium", "Polysphondylium", "Dictyostelium", "Cavenderia", "Heterostelium", "Pelomyxa", "Entamoeba", "Planoprotostelium", "Acanthamoeba", "Stygamoeba", "Vermistella", "Vexilliferidae", "Paramoeba", "Vannella", "Cochliopodium", "Amoeba_proteus", "Protostelium", "Filamoeba", "Mastigamoeba", "Arcella_intermedia", "Vermamoeba_vermiformis", "Dracoamoeba_jomungandri"),
  Archaeplastida = c("Arabidopsis_thaliana", "Chlamydomonas", "Volvox_carteri", "Ectocarpus_siliculosus", "Nitella_gaditana", "Physcomitrella_patens", "Porphyra_umbilicalis", "Selaginella_moellendorffii", "Madagascaria"),
  SAR = c("Toxoplasma_gondii", "Plasmodium_falciparum", "Phytophthora_infestans", "Colponemidia"),
  Metamonada = c("Trichomonas_vaginalis", "Giardia_lamblia", "Anaeramoeba"),
  Discoba = c("Naegleria", "Willaertia", "Neovahlkampfia", "Acrasis", "Pharyngomonas", "Tetramitus", "Naegleria_fowleri", "Vahlkampfia", "Allovahlkampfia", "Heteramoeba", "Soginia")
)

group_colors <- c(
  Opisthokonts = "#000080",  # Combined Metazoa, Unicellular Holozoa, and Fungi
  Apusozoa = "#615D44",
  Amoebozoa = "#00A79C",
  Archaeplastida = "#8CC540",
  SAR = "#8266AD",
  Metamonada = "#F2798E",
  Discoba = "#9A2268",
  Other = "black"
)

#------------------------------------------------------------
# Assign Groups to Tips
#------------------------------------------------------------
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

tip_data <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
group_assignments <- do.call(rbind, lapply(tip_data$label, assign_group))
tip_data <- cbind(tip_data, group_assignments)
tip_data$group <- factor(tip_data$group, levels = names(group_colors))

#------------------------------------------------------------
# Create the Base Plot with Colored Tip Labels
#------------------------------------------------------------
base_plot <- ggtree(rooted_tree, layout = "roundrect", size = 0.1) %<+% tip_data +
  geom_tiplab(aes(label = label, color = color), size = 1.5, align = FALSE, show.legend = FALSE) +
  scale_color_identity() +
  theme_tree2() 
#print(base_plot + ggtitle("base plot"))


#------------------------------------------------------------
# Add Clade Labels
#------------------------------------------------------------
clade_label <- function(label, tree, tips, offset = .6) {
  node <- getMRCA(tree, tips)
  geom_cladelab(node = node, label = label, fontsize = 3, align = F, barsize = 0.3, offset = offset)
}

plot_with_clades <- base_plot +
  clade_label("MyoII", rooted_tree, c("Colponemidia_sp_eukprot_EP00742_Colp-10_P003461", "Heteramoeba_clara_P25_DN11121_c0_g1_i1_p1"), offset = .2) +
  clade_label("MyoIA-E", rooted_tree, c("Dictyostelium_discoideum_SP14_XP_636359_2", "Homo_sapiens_SP14_NP_005370")) + 
  clade_label("MyoIF", rooted_tree, c("Homo_sapiens_SP14_NP_004989", "Naegleria_gruberi_SP14_XP_002680436")) + 
  clade_label("MyoV", rooted_tree, c("Homo_sapiens_SP14_NP_001073936", "Batrachochytrium_dendrobatidis_SP14_13697")) + 
  clade_label("MyoXI/XXIII", rooted_tree, c("Chlamydomonas_reinhardtii_SP14_XP_001693409", "Dictyostelium_discoideum_SP14_XP_001733030")) + 
  clade_label("MyoVI", rooted_tree, c("Homo_sapiens_SP14_NP_004990", "Capsaspora_owczarzaki_SP14_02007_2")) + 
  clade_label("MyoXIV", rooted_tree, c("Toxoplasma_gondii_SP14_XP_002368942", "Toxoplasma_gondii_SP14_XP_002366825")) + 
  clade_label("MyoXVII", rooted_tree, c("Batrachochytrium_dendrobatidis_SP14_8502", "Batrachochytrium_dendrobatidis_SP14_35424"))  +
  clade_label("MyoXXII", rooted_tree, c("Toxoplasma_gondii_SP14_XP_002367823", "Toxoplasma_gondii_SP14_XP_002364981")) + 
  clade_label("MyoIII-like", rooted_tree, c("Homo_sapiens_SP14_NP_059129", "Capsaspora_owczarzaki_SP14_02007_8")) + 
  clade_label("MyoIX", rooted_tree, c("Homo_sapiens_SP14_NP_004136", "Capsaspora_owczarzaki_SP14_02007_7")) + 
  clade_label("MyoXV", rooted_tree, c("Homo_sapiens_SP14_NP_057323", "Capsaspora_owczarzaki_SP14_02007_13")) + 
  clade_label("MyoXXV", rooted_tree, c("Dictyostelium_discoideum_SP14_XP_644171", "Capsaspora_owczarzaki_SP14_02007_15")) + 
  clade_label("MyoXIII", rooted_tree, c("Naegleria_gruberi_SP14_XP_002681567", "Naegleria_gruberi_SP14_XP_002680898")) + 
  clade_label("MyoI-like", rooted_tree, c("Naegleria_gruberi_SP14_XP_002676946", "Naegleria_gruberi_SP14_XP_002683472"))  
#  clade_label("", rooted_tree, c("", "")) + 

#print(plot_with_clades + ggtitle("plot with clades"))

#------------------------------------------------------------
# Add Bootstrap values
#------------------------------------------------------------
# Extract bootstrap values and associate with node numbers
bootstrap_values <- as.numeric(tree$node.label)
node_numbers <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
bootstrap_data <- data.frame(node = node_numbers, bootstrap = bootstrap_values, stringsAsFactors = FALSE)

# Merge bootstrap data with the tree data
plot_with_clades$data <- left_join(plot_with_clades$data, bootstrap_data, by = "node")

#find bootstrap for specific nodes of intertest to display as numbers
nodes_of_interest <- c(
  bootstrap_data[bootstrap_data$node == getMRCA(rooted_tree, c("Thecamonas_trahens_eukprot_EP00031_P006505", "Aspergillus_nidulans_eukprot_EP00135_P004703")),],
  bootstrap_data[bootstrap_data$node == getMRCA(rooted_tree, c("Filamoeba_sp_ATCC50430_eukprot_EP00027_P007255", "Arcella_intermedia_eukprot_EP01124_P025971")),],
  bootstrap_data[bootstrap_data$node == getMRCA(rooted_tree, c("Naegleria_fowleri_MhcC_NCBI_XP_044562059_1", "Pharyngomonas_kirbyi_eukprot_EP00761_P011290")),],
  bootstrap_data[bootstrap_data$node == getMRCA(rooted_tree, c("Rigifila_ramosa_eukprot_EP00004_P018167", "Pharyngomonas_kirbyi_eukprot_EP00761_P011290")),]
)




# Add circles to internal nodes colored by bootstrap support
bootstrap_threshold <- 75
plot_with_bootstrap <- plot_with_clades +
  geom_point(
    data = subset(plot_with_clades$data, !isTip & !is.na(bootstrap) & bootstrap > bootstrap_threshold),
    aes(x = x, y = y, fill = bootstrap),
    shape = 21, color = "black", size = 1  # Set size to 1
  ) +
  geom_text(
    data = subset(plot_with_clades$data, node %in% nodes_of_interest),
    aes(x = x, y = y, label = bootstrap),
    nudge_x = 0.08,  # Adjust as needed to position text above nodes
    size = 2,
    color = "red"
  ) +
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap", limits = c(bootstrap_threshold, 100)) +
  theme(legend.position = c(0.08, 0.6))  

print(plot_with_bootstrap + ggtitle("plot with bootstrap"))


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
print(final_plot + ggtitle("Eukaryotic Myosins"))

