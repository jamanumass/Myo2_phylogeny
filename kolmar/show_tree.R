getwd()
installed.packages()
source("https://bioconductor.org/biocLite.R")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## see below BiocManager::install("ggtree")
BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)
library(treeio)
library(tidytree)
library(ape)



tr <- read.iqtree(file = "kol_select_sp_align.fas.treefile")
 ggtree(tr)
 
getwd()
file <- system.file("","kol_select_sp_align.fas.treefile", package = "treeio")
iqtree <- read.iqtree(file)


tr <- read.iqtree("kol_select_sp_align.fas.iqtree")

file <- system.file("kol_select_sp_align.fas.iqtree"  ,package="treeio")

#having all sorts of difficulty using the system.file() to get the tree, and then read.iqrtee() to make the object. Moving on

tr<-read.tree(file = "kol_select_sp_align.fas.treefile")
ggtree(tr, layout="daylight", branch.length = 'none') + geom_tiplab(size = 3)

#vignette

library(ape)
set.seed(2017)
tree <- rtree(4)
tree
x <- as_tibble(tree)
x
as.phylo(x)
d <- tibble(label = paste0('t', 1:4),
            trait = rnorm(4))

y <- full_join(x, d, by = 'label')
y
as.treedata(y)
y <- y %>% as.treedata %>% as_tibble

child(y, 5)
# something broken above, moving on


nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

ggtree(tree, color = "red", size = 2, linetype = "dotted")
ggtree(tree) + geom_point(aes(shape = isTip, color = isTip), size = 3)
ggtree(tree)  + 
  geom_nodelab(aes(label = branch.length), offset(1), color = "red") + 
  geom_tippoint(size = 30, color = "yellow", alpha = .2) +
  geom_tiplab(aes(angle = 0),color = "blue",size = 3) + 
  geom_rootedge()
