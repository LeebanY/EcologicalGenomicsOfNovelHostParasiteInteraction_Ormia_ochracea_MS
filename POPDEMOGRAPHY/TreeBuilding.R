
## ------ making a tree
library(ggtree)
library(treeio)
library(tidytree)
library(dplyr)
library(ape)


tree <- read.tree("Ormia_fasttree_phylogeny.tree")

plot(tree)
# Example: Replace tip labels
# Convert the tree to a tibble for easier manipulation
tree_tibble <- as_tibble(tree)

tree$tip.label

labels_tibble <- tree_tibble %>%
  left_join(pca_all, by = c("label" = "individuals")) %>%
  mutate(label = ifelse(!is.na(Sample_ID), Sample_ID, Sample_ID)) %>%
  select(-Sample_ID) %>%
  select(parent, node, branch.length, label)# Remove extra column if desired

labels_tree<-as.phylo(labels_tibble)

# Identify tips that do NOT contain "KAUAI"
tips_to_remove <- labels_tree$tip.label[!grepl("KAUAI", labels_tree$tip.label)]

# Create a new tree with only the desired tips
filtered_tree <- drop.tip(labels_tree, tips_to_remove)


# Plot the tree with branch colors
ggtree(filtered_tree, branch.length='none', layout="rectangular") +
  geom_tiplab(size=1) +
  theme(legend.position = "none")

# Join the tree data with the PCA dataset based on the old tip labels
updated_tree <- tree_tibble %>%
  left_join(pca, by = c("label" = "individuals")) %>%
  mutate(label = ifelse(!is.na(Location), Location, Location)) %>%
  select(-Location) %>%
  select(parent, node, branch.length, label)# Remove extra column if desired

# Replace specified values with "Mainland"
updated_tree$label <- gsub("SANTABARBARA|VENTURA|LOSANGELES", "CALIFORNIA", updated_tree$label)
updated_tree$label <- gsub("YAVAPAI|COCHISE|SANTACRUZ", "ARIZONA", updated_tree$label)

# Convert the updated tibble back to a phylogenetic tree
updated_tree <- as.phylo(updated_tree)

# Find any YAVAPAI individual to use as the outgroup
yavapai_tips <- updated_tree$tip.label[grepl("ARIZONA", updated_tree$tip.label)]

# Check if we found any YAVAPAI individual
if (length(yavapai_tips) > 0) {
  new_root <- yavapai_tips[1]  # Select the first YAVAPAI individual
  updated_tree <- root(updated_tree, outgroup = new_root, resolve.root = TRUE)
} else {
  warning("No YAVAPAI individual found in the tree!")
}

# Plot the rerooted tree to visualize
plot(updated_tree, show.tip.label = TRUE)


# Drop the identified tips
color_map <- c(
  "KAUAI" = "firebrick1",
  "OAHU" = "deepskyblue1",
  "HILO" = "darkgreen",
  "MOLOKAI" = "darkorchid1",
  "ARIZONA" = "orange",
  "CALIFORNIA" = "pink")


# Ensure the labels in the tree match your color map
# If some labels don't have a color, assign "black" as the default
tree_tibble <- as_tibble(updated_tree) %>%
  mutate(color = ifelse(label %in% names(color_map), color_map[label], "black"))

# Convert the tibble back to a tree

tree_colored <- as.phylo(tree_tibble)


# Plot the tree with branch colors
PrettyFastTreeNoBranchLengths<-ggtree(tree_colored, aes(color = label), branch.length='none', layout="rectangular") +
  geom_tiplab(aes(color = label), size=0) +
  scale_color_manual(values = c(color_map, "black")) +
  theme(legend.position = "none")


# Plot the tree with branch lengths
PrettyFastTreeWITHBranchLengths<-ggtree(tree_colored, aes(color = label), layout="rectangular") +
  geom_tiplab(aes(color = label), size=0) +
  scale_color_manual(values = c(color_map, "black")) +
  theme(legend.position = "none")+
  geom_treescale(offset=2)



png("~/Desktop/Ochracea_POSTDOC/FastTree_GTR_ColouredBranchLlengths.png", units="in", width=7.0, height=10.0, res=300)
PrettyFastTreeWITHBranchLengths
dev.off()


#---- pbs tree -- all hits together. 

library(ggtree)
library(treeio)
library(tidytree)
library(dplyr)
library(ape)


tree <- read.tree("PBS_tree.min200.tree")

# Example: Replace tip labels
# Convert the tree to a tibble for easier manipulation
tree_tibble <- as_tibble(tree)

tree$tip.label
# Join the tree data with the PCA dataset based on the old tip labels
updated_tree <- tree_tibble %>%
  left_join(pca, by = c("label" = "individuals")) %>%
  mutate(label = ifelse(!is.na(Location), Location, Location)) %>%
  select(-Location) %>%
  select(parent, node, branch.length, label)# Remove extra column if desired

# Replace specified values with "Mainland"
updated_tree$label <- gsub("YAVAPAI|COCHISE|SANTABARBARA|SANTACRUZ|VENTURA|LOSANGELES", "Mainland", updated_tree$label)


# Convert the updated tibble back to a phylogenetic tree
updated_tree <- as.phylo(updated_tree)


# Drop the identified tips
color_map <- c(
  "KAUAI" = "firebrick1",
  "OAHU" = "deepskyblue1",
  "HILO" = "darkgreen",
  "MOLOKAI" = "darkorchid1",
  "Mainland" = "orange"
  
)


# Ensure the labels in the tree match your color map
# If some labels don't have a color, assign "black" as the default
tree_tibble <- as_tibble(updated_tree) %>%
  mutate(color = ifelse(label %in% names(color_map), color_map[label], "black"))

# Convert the tibble back to a tree

tree_colored <- as.phylo(tree_tibble)


# Plot the tree with branch colors
ggtree(tree_colored, aes(color = label), branch.length='none', layout="rectangular") +
  geom_tiplab(aes(color = label), size=0) +
  scale_color_manual(values = c(color_map, "black")) +
  theme(legend.position = "right")



# function
library(ggtree)
install.packages("ggtree")
library(treeio)
install.packages("treeio")
library(tidytree)

library(dplyr)
library(ape)
install.packages("ape")
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library(ggtree)
library(treeio)
library(tidytree)
library(dplyr)
library(ape)
library(ggtree)
library(treeio)
library(tidytree)
library(dplyr)
library(ape)
library(ggpubr)

plot_tree <- function(tree, pca_data) {
  # Convert the tree to a tibble for easier manipulation
  tree_tibble <- as_tibble(tree)
  
  # Join with PCA dataset and update labels
  updated_tree <- tree_tibble %>%
    left_join(pca_data, by = c("label" = "individuals")) %>%
    mutate(label = ifelse(!is.na(Location), Location, label)) %>%
    select(-Location) %>%
    select(parent, node, branch.length, label)
  
  # Replace specified values with "Mainland"
  updated_tree$label <- gsub("YAVAPAI|COCHISE|SANTABARBARA|SANTACRUZ|VENTURA|LOSANGELES", "Mainland", updated_tree$label)
  
  # Assign "Outgroup" to tips without labels
  tip_nodes <- tree$tip.label
  updated_tree <- updated_tree %>%
    mutate(label = ifelse(is.na(label) & node %in% (1:length(tip_nodes)), "Outgroup", label))
  
  # Convert the updated tibble back to a phylogenetic tree
  updated_tree <- as.phylo(updated_tree)
  
  # Define color map
  color_map <- c(
    "KAUAI" = "firebrick1",
    "OAHU" = "deepskyblue1",
    "HILO" = "darkgreen",
    "MOLOKAI" = "darkorchid1",
    "Mainland" = "orange",
    "Outgroup" = "black"
  )
  
  # Ensure the labels in the tree match your color map
  tree_tibble <- as_tibble(updated_tree) %>%
    mutate(color = ifelse(label %in% names(color_map), color_map[label], "black"))
  
  # Convert back to tree
  tree_colored <- as.phylo(tree_tibble)
  
  # Generate and return the tree plot
  plot <- ggtree(tree_colored, aes(color = label), layout="daylight") +
    geom_tiplab(aes(color = label), size = 0) +
    scale_color_manual(values = c(color_map, "black")) +
    theme(legend.position = "right")
  
  return(plot)
}


#### chr1 tree
chr1tree <- read.tree("separate_pbs_trees/PBS_hits_trees/chr1_hits.min200.fasta.tree")
chr1tree_plot <- plot_tree(chr1tree, pca_all)+labs(title="Chr 1 variants")

#### chr2 tree
chr2tree <- read.tree("separate_pbs_trees/PBS_hits_trees/chr2_hits.min200.fasta.tree")
chr2tree_plot <- plot_tree(chr2tree, pca_all)+labs(title="Chr 2 variants")


#### chr3 tree
chr3tree <- read.tree("separate_pbs_trees/PBS_hits_trees/chr3_hits.min200.fasta.tree")
chr3tree_plot <- plot_tree(chr3tree, pca_all)+labs(title="Chr 3 variants")


#### chr4 tree
chr4tree <- read.tree("separate_pbs_trees/PBS_hits_trees/chr4_hits.min200.fasta.tree")
chr4tree_plot <- plot_tree(chr4tree, pca_all)+labs(title="Chr 4 variants")


TreesPBSOutliersAltogether<-ggarrange(chr1tree_plot, chr2tree_plot, chr3tree_plot, chr4tree_plot, ommon.legend = TRUE, legend="bottom", nrow=2, ncol=2)
