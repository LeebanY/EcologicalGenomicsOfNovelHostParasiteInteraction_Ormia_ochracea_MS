# Mainland Flies -- what kind of population structure do they have? Do they group by what they parasitize owing to local adaptation/cryptic divergence or do they group by location.library(tidyverse)
library(viridis)
library(ggpubr)
library(data.table)
library(ggExtra)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Ok import relevant files into the R script.
###### PCA of Ormia
pca_mainland <- read.table("MainlandStructure/Mainland_Ormia_pruned_input.eigenvec",header = F)
eigenval_mainland <- scan('MainlandStructure/Mainland_Ormia_pruned_input.eigenval')
# import pop file.
pops<-read_tsv("populations_file.txt", col_names = c("Samples", "Population"))

#Get rid of duplicated column 
pca_mainland <- pca_mainland[,-1]
names(pca_mainland)[1] <- "individuals"
names(pca_mainland)[2:ncol(pca_mainland)] <- paste0("PC", 1:(ncol(pca_mainland)-1))

glimpse(pca_mainland$individuals)

#importantly we want to create a variable with real individual names.
pca_mainland$Sample_ID<-pops$Population[match(pca_mainland$individuals, pops$Samples)]
#OK now we need to create an additional column with location names to group them by.
pca_mainland$Location <- str_extract(pca_mainland$Sample_ID, "^[^_]+") 
pca_mainland$Site <- str_extract(pca_mainland$Sample_ID, "[^_]+$")
pca_mainland <- pca_mainland %>%
  mutate(
    Groups = case_when(
      Location == "YAVAPAI" | Location == "COCHISE" | Location == "SANTACRUZ" ~ "Arizona",
      Location == "SANTABARBARA" | Location == "VENTURA" | Location == "LOSANGELES" ~ "California"))


# Now I want to check the PC axes to see how much variation is captured and explained.
pve_main <- data.frame(PC = 1:10, pve_main = eigenval_mainland/sum(eigenval_mainland)*100)
a <- ggplot(pve_main, aes(PC, pve_main)) + geom_bar(stat = "identity")
a <- a + ylab("Percentage variance explained") + theme_light()

# -- make the PCA. 
main_palette <- c(
  "Arizona" = "orange",   # Dark Blue
  "California" = "lightpink1" # Red
)

pca_mainland$Groups <- factor(pca_mainland$Groups , levels = names(main_palette))

ggplot(pca_mainland, aes(PC1, PC2, col = Groups, shapes=Site))+
  geom_vline(xintercept = 0, linetype='dashed', colour='darkgrey', alpha=1)+
  geom_hline(yintercept = 0, linetype='dashed', colour='darkgrey', alpha=1)+
  geom_point(size = 3, alpha=1)

b <- b + coord_equal() + theme_bw() + theme(legend.position="none")
b <- b + scale_colour_manual(values = my_palette)

b <-b + xlab(paste0("PC1 (", signif(pve_all$pve_all[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_all$pve_all[2], 3), "%)"))



