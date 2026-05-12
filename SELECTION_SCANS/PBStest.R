# ==============================================================================
# PBS (POPULATION BRANCH STATISTIC) AND TAJIMA'S D ANALYSIS
# Detecting selective sweeps in Hawaiian Ormia ochracea populations
# Comparing Hawaii vs California vs Arizona populations
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(tidyverse)      # Data manipulation and plotting
library(viridis)        # Colour palettes
library(ggpubr)         # Publication-ready ggplot themes
library(data.table)     # Fast data manipulation (used for interval joins)
library(ggExtra)        # Marginal plots
library(RColorBrewer)   # Colour palettes
library(cowplot)        # Plot composition
library(gtools)         # mixedsort() for natural chromosome ordering
library(rstatix)        # Tidy statistical tests (Wilcoxon, effect sizes)

# ==============================================================================
# SECTION 1: LOAD FST DATA AND COMPUTE PBS
# ==============================================================================
# FST values are computed in windows across the genome between three populations:
#   H = Hawaii, A = Arizona (outgroup), C = California (outgroup)
# PBS isolates the branch leading to Hawaii using the three-population design

# Load windowed FST data (three pairwise comparisons)
FST_threeway <- read_table("Hawaii_California_Arizona.allchrs.filtered.repeatsremoved.fst")

# Clean chromosome names
FST_threeway$CHROM_CvsA <- gsub("chr", "", FST_threeway$CHROM_CvsA)

# Report mean FST across comparisons
mean(na.omit(FST_threeway$WEIGHTED_FST_HvsA))
mean(na.omit(FST_threeway$WEIGHTED_FST_HvsC))

# Compute PBS for each window:
# PBS = [(T_HvsA + T_HvsC) - T_CvsA] / 2
# where T = -log(1 - FST), transforming FST to branch length units
PBS_threeway <- FST_threeway %>%
  select(CHROM_CvsA, BIN_START_CvsA, BIN_END_CvsA,
         WEIGHTED_FST_CvsA, WEIGHTED_FST_HvsA, WEIGHTED_FST_HvsC,
         N_VARIANTS_CvsA, N_VARIANTS_HvsA, N_VARIANTS_HvsC) %>%
  dplyr::rename(
    chromosome = CHROM_CvsA,
    start      = BIN_START_CvsA,
    end        = BIN_END_CvsA
  ) %>%
  # Require at least 100 variants per window in all three comparisons
  filter(N_VARIANTS_HvsA >= 100 & N_VARIANTS_CvsA >= 100 & N_VARIANTS_HvsC >= 100) %>%
  # Convert FST to branch length units
  mutate(
    T_A = -log(1 - WEIGHTED_FST_HvsA),   # Hawaii vs Arizona branch
    T_B = -log(1 - WEIGHTED_FST_HvsC),   # Hawaii vs California branch
    T_C = -log(1 - WEIGHTED_FST_CvsA)    # California vs Arizona branch (outgroup comparison)
  ) %>%
  # Compute PBS for each population branch
  mutate(
    pbs  = ((T_A + T_B) - T_C) / 2,     # Standard PBS for Hawaii
    PBSa = (T_A + T_B - T_C) / 2,       # PBS for Hawaii (branch A)
    PBSb = (T_B + T_C - T_A) / 2,       # PBS for California (branch B)
    PBSc = (T_A + T_C - T_B) / 2        # PBS for Arizona (branch C)
  ) %>%
  # Compute PBSn1A: normalised PBS to control for variable mutation rates
  # PBSn1A = PBSa / (1 + PBSa + PBSb + PBSc)
  mutate(PBSn1A = PBSa / (1 + PBSa + PBSb + PBSc)) %>%
  arrange(-PBSn1A)

# Remove X chromosome from analysis
PBS_threeway <- PBS_threeway %>% filter(chromosome != "X")

# ==============================================================================
# SECTION 2: MANHATTAN PLOT — PBSn1A ACROSS THE GENOME
# ==============================================================================

# Define outlier threshold (95th percentile)
threshold <- 0.2328979

# Manhattan plot: PBSn1A coloured by significance and alternating chromosome shading
PBS_threeway_plot <- PBS_threeway %>%
  na.omit() %>%
  mutate(
    chr_position = (start + end) / 2,
    chromosome   = factor(chromosome, levels = gtools::mixedsort(unique(chromosome))),
    point_color  = case_when(
      PBSn1A > threshold                         ~ "red3",    # Outlier windows
      as.numeric(factor(chromosome)) %% 2 == 0  ~ "grey80",  # Even chromosomes
      TRUE                                      ~ "grey60"   # Odd chromosomes
    )
  ) %>%
  ggplot(aes(x = chr_position, y = PBSn1A, colour = point_color)) +
  geom_jitter(size = 2.5, alpha = 1) +
  facet_grid(cols = vars(chromosome), space = "free_x", scales = "free_x", switch = "x") +
  geom_hline(yintercept = threshold, linetype = 'dashed') +
  scale_colour_identity() +
  labs(x = "", y = expression(PBSn1[Hawaii])) +
  theme_classic() +
  theme(
    axis.text        = element_text(size = 17),
    axis.title       = element_text(size = 17, face = "bold"),
    axis.text.x      = element_text(size = 0),
    axis.ticks.x     = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(size = 17),
    legend.position  = 'none'
  ) +
  ylim(0, 1)

# ==============================================================================
# SECTION 3: IDENTIFY AND EXPORT PBS OUTLIER WINDOWS
# ==============================================================================

# 95th percentile threshold
my_threshold <- quantile(PBS_threeway$PBSn1A, 0.95, na.rm = T)

# Label windows as outlier or background
PBS_threeway <- PBS_threeway %>%
  mutate(outlier = ifelse(PBSn1A > my_threshold, "outlier", "background"))

# Count outlier windows per chromosome
PBS_threeway %>%
  filter(outlier == "outlier") %>%
  group_by(chromosome) %>%
  tally()

# Export outlier windows
OutliersPbsn1a <- PBS_threeway %>% filter(outlier == "outlier")
write.csv(OutliersPbsn1a, "OutliersPbsn1a.csv")

# ==============================================================================
# SECTION 4: PHYLOGENETIC TREE — PBS TOPOLOGY
# ==============================================================================
# Visualise the population tree to show the branching structure underlying PBS

library(ggtree)    # Tree visualisation
library(treeio)    # Tree I/O
library(tidytree)  # Tidy tree manipulation
library(ape)       # Phylogenetic tree utilities

# Load the population tree
tree <- read.tree("PBS_tree.min200.tree")

# Update tip labels using PCA metadata to replace sample IDs with location names
tree_tibble <- as_tibble(tree)

updated_tree <- tree_tibble %>%
  left_join(pca, by = c("label" = "individuals")) %>%
  mutate(label = ifelse(!is.na(Location), Location, Location)) %>%
  select(parent, node, branch.length, label)

# Consolidate all mainland US populations into a single "Mainland" label
updated_tree$label <- gsub(
  "YAVAPAI|COCHISE|SANTABARBARA|SANTACRUZ|VENTURA|LOSANGELES",
  "Mainland",
  updated_tree$label
)

# Convert updated tibble back to phylo object
updated_tree <- as.phylo(updated_tree)

# Colour map for populations
color_map <- c(
  "KAUAI"    = "firebrick1",
  "OAHU"     = "deepskyblue1",
  "HILO"     = "darkgreen",
  "MOLOKAI"  = "darkorchid1",
  "Mainland" = "orange"
)

# Add colour assignments to tree tibble
tree_tibble <- as_tibble(updated_tree) %>%
  mutate(color = ifelse(label %in% names(color_map), color_map[label], "black"))

tree_colored <- as.phylo(tree_tibble)

# Plot coloured tree
ggtree(tree_colored, aes(color = label)) +
  geom_tiplab(aes(color = label), size = 0) +
  scale_color_manual(values = c(color_map, "black")) +
  theme(legend.position = "right")

# ==============================================================================
# SECTION 5: OVERLAY GENE MODELS ON PBS PLOT
# ==============================================================================
# Import GFF file of genes in shared sweep regions and overlay on PBS Manhattan plot
# (Section also contains a commented-out alternative approach using GenomicRanges)

library(rtracklayer)   # GFF/GTF import

# Import GFF of genes overlapping shared selective sweep regions
gff    <- import("EvidenceOfSelectiveSweeps.SharedAcrossTests.genes.gff")
gff_df <- as.data.frame(gff)

# Extract gene-level features and match chromosome factor levels to PBS plot
genes <- gff_df %>%
  filter(type == "gene") %>%
  select(chromosome = seqnames, start, end, strand, Name) %>%
  mutate(
    chromosome = factor(
      chromosome,
      levels = levels(PBS_threeway_plot$data$chromosome)  # Ensures exact match with plot
    )
  )

# Overlay gene positions as horizontal segments below the PBS Manhattan plot
PBS_threeway_plot +
  geom_segment(
    data        = genes,
    aes(x = start, xend = end, y = -0.05, yend = -0.05, colour = "blue"),
    inherit.aes = FALSE,
    size        = 1
  )

# ==============================================================================
# SECTION 6: OVERLAY DIFFERENTIALLY EXPRESSED GENES (DEGs) ON PBS PLOT
# ==============================================================================
# Visualise whether DEG positions coincide with PBS sweep outlier regions

# Load DEG genomic coordinates
deg_gff <- read.table("DEG_ORMIA_HEADS_NOLOGFCCUTOFF_padj0.05.gff",
                      header = FALSE, comment.char = "#", sep = "\t")

deg_gff <- deg_gff %>%
  select(V1, V4, V5) %>%
  mutate(chr_position = (V4 + V5) / 2)

colnames(deg_gff) <- c("chromosome", "start", "end", "chr_position")

# PBS plot with DEG positions overlaid as vertical lines
PBS_withDEG <- PBS_threeway %>%
  na.omit() %>%
  filter(chromosome != "X") %>%
  mutate(chr_position = (start + end) / 2) %>%
  ggplot(aes(x = chr_position, y = PBSn1A)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'lightgrey') +
  labs(x = "Scaffold (MB)", y = expression(PBSn1[A])) +
  facet_grid(cols = vars(as.numeric(chromosome)), space = "free_x", scales = "free_x", switch = "x") +
  geom_hline(yintercept = 0.2328979, linetype = 'dashed') +
  theme_classic() +
  theme(
    axis.text        = element_text(size = 14),
    axis.title       = element_text(size = 13, face = "bold"),
    axis.text.x      = element_text(size = 0),
    axis.ticks.x     = element_blank(),
    legend.position  = 'none'
  ) +
  # Vertical lines mark positions of DEGs
  geom_vline(
    data        = deg_gff,
    aes(xintercept = chr_position),
    inherit.aes = FALSE,
    linetype    = "solid",
    color       = "darkorchid4"
  )

png("Pbs_with_degenes.png", units = "in", width = 14.0, height = 3.0, res = 300)
PBS_withDEG
dev.off()

# PBS plot with GEA significant SNPs overlaid as vertical dashed lines
PBS_threeway %>%
  na.omit() %>%
  filter(chromosome != "X") %>%
  mutate(chr_position = (start + end) / 2) %>%
  ggplot(aes(x = chr_position, y = PBSn1A)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'lightgrey') +
  labs(x = "Scaffold (MB)", y = expression(PBSn1[A])) +
  facet_grid(cols = vars(as.numeric(chromosome)), space = "free_x", scales = "free_x", switch = "x") +
  geom_hline(yintercept = 0.2328979, linetype = 'dashed') +
  theme_classic() +
  theme(
    axis.text       = element_text(size = 14),
    axis.title      = element_text(size = 13, face = "bold"),
    axis.text.x     = element_text(size = 0),
    axis.ticks.x    = element_blank(),
    legend.position = 'none'
  ) +
  # Vertical lines mark positions of GEA silence outlier SNPs
  geom_vline(
    data        = ridge_outliers_silence,
    aes(xintercept = chr_position),
    inherit.aes = FALSE,
    color       = "red",
    linetype    = "dashed",
    alpha       = 0.7
  )

# ==============================================================================
# SECTION 7: GENE ONTOLOGY ENRICHMENT — PBS SWEEP GENES
# ==============================================================================
# Test whether genes in sweep regions are enriched for particular GO functions
# relative to the full annotated genome, using Fisher's exact test

# Load full genome gene function counts and sweep gene function counts
full_gene_set <- read.table("OrmiaCountsofAllGeneFunctions.tsv",
                            header = T, na.strings = "NA", fill = TRUE)
full_gene_set$Function[full_gene_set$Function == ""] <- NA
full_gene_set <- na.omit(full_gene_set)

genes_of_interest <- read.table("EvidenceOfSelectiveSweeps.SharedAcrossTests.genes.enrichment",
                                header = T, na.strings = "NA", fill = TRUE)
genes_of_interest$Function[genes_of_interest$Function == ""] <- NA
genes_of_interest <- na.omit(genes_of_interest)

# Merge and compute expected counts under null (proportional representation)
merged_data <- full_gene_set %>%
  left_join(genes_of_interest, by = "Function", suffix = c("_total", "_interest")) %>%
  mutate(
    Counts_interest = replace_na(Counts_interest, 0),
    Expected        = Counts_total * sum(Counts_interest) / sum(Counts_total)
  )

# Run Fisher's exact test for each GO function
results <- lapply(1:nrow(merged_data), function(i) {
  matrix(
    c(merged_data$Counts_interest[i], merged_data$Expected[i],
      sum(merged_data$Counts_interest) - merged_data$Counts_interest[i],
      sum(merged_data$Counts_total)    - merged_data$Expected[i]),
    nrow = 2, byrow = TRUE
  ) %>% fisher.test()
})

# Collect p-values and apply BH correction
merged_data$p_value     <- sapply(results, function(x) x$p.value)
merged_data$adj_p_value <- p.adjust(merged_data$p_value, method = "BH")

# Extract significantly enriched GO terms and compute enrichment ratio
Significantly_enriched_GO_terms_PosSelection <- merged_data %>%
  filter(adj_p_value < 0.05) %>%
  mutate(enrichment = Counts_interest / Expected)

# Clean up function labels for display
Significantly_enriched_GO_terms_PosSelection$Function <- gsub(
  "_", " ", Significantly_enriched_GO_terms_PosSelection$Function)

# Top 20 most enriched GO terms (by lowest adjusted p-value)
Top30forPlottingPBS <- Significantly_enriched_GO_terms_PosSelection %>%
  top_n(-20, adj_p_value) %>%
  arrange(desc(enrichment))

# Bar chart of enrichment scores for top GO terms
GeneOntologyPBS <- ggplot(Top30forPlottingPBS, aes(x = reorder(Function, enrichment), y = enrichment)) +
  geom_col(fill = 'grey', colour = 'black') +
  coord_flip() +
  theme_pubr() +
  labs(y = "Enrichment score", x = "", title = "") +
  theme(
    axis.text  = element_text(size = 15),
    axis.title = element_text(size = 17),
    strip.text = element_text(size = 17)
  )

# Scatter plot: enrichment vs significance for all enriched GO terms
ggplot(Significantly_enriched_GO_terms_PosSelection,
       aes(x = -log10(adj_p_value), y = enrichment)) +
  geom_point(size = 2) +
  geom_text(aes(label = Function), max.overlaps = 100) +
  coord_flip() +
  theme_pubr() +
  labs(y = "Enrichment score", x = "-log10 (p-value)", title = "") +
  theme(
    axis.text  = element_text(size = 15),
    axis.title = element_text(size = 17),
    strip.text = element_text(size = 17)
  )

png("GeneOntologyPBS.png", units = "in", width = 7.5, height = 6.0, res = 300)
GeneOntologyPBS
dev.off()

# ==============================================================================
# SECTION 8: TAJIMA'S D — OVERLAP WITH SHARED SWEEP REGIONS
# ==============================================================================
# Mark Tajima's D windows that overlap with shared sweep regions (identified
# across multiple sweep detection methods), then compare Tajima's D between
# sweep regions and background genome

######## Tajimas D. 

# Kauai.

TajimasDKauai<-read_table("TajimaD_Downloand/KAUAI.TajimaD.final")

TajimasDKauai$CHROM<-gsub("chr", "", TajimasDKauai$CHROM)

TajimasDKauai<-TajimasDKauai %>%
  filter(N_SNPS > 100 & !CHROM=="X")


TajimasDKauai$TajimaD<-as.double(TajimasDKauai$TajimaD)

mean(TajimasDKauai$TajimaD)

HistogramTajimasD_Kauai<-TajimasDKauai %>%
  ggplot(aes(x=TajimaD))+
  geom_histogram(fill='grey', colour='black')+
  theme_pubr()+
  labs(x=expression(paste("Tajima's D - Kauai", italic(" D"))),
       y="count (10kb windows)")+
  geom_vline(xintercept = 0, linetype='dashed', colour='black')+
  geom_vline(xintercept = 2.979988, size = 1, linetype='dashed', colour='tomato1')



png("TajimasDHistogram.png", units="in", width=5.0, height=3.0, res=300)
HistogramTajimasD_Kauai
dev.off()


GenomewideTajDKauai_plot<-TajimasDKauai %>%
  na.omit() %>%
  ggplot(aes(x=BIN_START, y=TajimaD))+
  geom_jitter(size=2.5,alpha=0.7, colour='lightgrey')+
  geom_smooth(method='loess', alpha=0.6, span=0.1, se=F, colour='tomato1')+
  labs(y="Kaua'i")+
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Scaffold (MB)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=13,face="bold"))+
  theme_classic()+
  theme(legend.position = 'none')+
  theme(axis.text.y=element_text(size=13))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
  #geom_hline(yintercept = 1, linetype='dashed', colour='red', size=1)+
  geom_hline(yintercept = -2, linetype='dashed', colour='red', size=1)


# OAHU 


TajimasDOahu<-read_table("TajimaD_Downloand/OAHU.TajimaD.final")

TajimasDOahu$CHROM<-gsub("chr", "", TajimasDOahu$CHROM)

TajimasDOahu<-TajimasDOahu %>%
  filter(N_SNPS > 100)

TajimasDOahu$TajimaD<-as.double(TajimasDOahu$TajimaD)


mean(TajimasDOahu$TajimaD)

HistogramTajimasD_Oahu<-TajimasDOahu %>%
  ggplot(aes(x=TajimaD))+
  geom_histogram(fill='grey', colour='black')+
  theme_pubr()+
  labs(x=expression(paste("Tajima's D - Oahu", italic(" D"))),
       y="count (10kb windows)")+
  geom_vline(xintercept = 0, linetype='dashed', colour='black')+
  geom_vline(xintercept = 1.539397, size = 1, linetype='dashed', colour='dodgerblue')


png("TajimasDHistogram.png", units="in", width=5.0, height=3.0, res=300)
HistogramTajimasD
dev.off()



TajimasDOahu_plot<-TajimasDOahu %>%
  na.omit() %>%
  ggplot(aes(x=BIN_START, y=TajimaD))+
  geom_jitter(size=2.5,alpha=0.7, colour='lightgrey')+
  geom_smooth(method='loess', alpha=0.6, span=0.1, se=F, colour='dodgerblue')+
  labs(y="O'ahu")+
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Scaffold (MB)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=13,face="bold"))+
  theme_classic()+
  theme(legend.position = 'none')+
  theme(axis.text.y=element_text(size=13))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
  #geom_hline(yintercept = 1, linetype='dashed', colour='red', size=1)+
  geom_hline(yintercept = -2, linetype='dashed', colour='red', size=1)


# OAHU 


TajimasDHilo<-read_table("TajimaD_Downloand/HILO.TajimaD.final")

TajimasDHilo$CHROM<-gsub("chr", "", TajimasDHilo$CHROM)

TajimasDHilo<-TajimasDHilo %>%
  filter(N_SNPS > 100)

TajimasDHilo$TajimaD<-as.double(TajimasDHilo$TajimaD)

mean(TajimasDHilo$TajimaD)

HistogramTajimasD_Hilo<-TajimasDHilo %>%
  ggplot(aes(x=TajimaD))+
  geom_histogram(fill='grey', colour='black')+
  theme_pubr()+
  labs(x=expression(paste("Hilo", italic(" D"))),
       y="count (10kb windows)")+
  geom_vline(xintercept = 0, linetype='dashed', colour='black')+
  geom_vline(xintercept = 1.708403, size = 1, linetype='dashed', colour='forestgreen')


png("TajimasDHistogram.png", units="in", width=5.0, height=3.0, res=300)
HistogramTajimasD_Hilo
dev.off()

Histograms_ofTajimasD<-ggarrange(HistogramTajimasD_Kauai, HistogramTajimasD_Oahu, HistogramTajimasD_Hilo, ncol=3, nrow=1)

png("TajimasDHistogram.png", units="in", width=12.0, height=3.0, res=300)
Histograms_ofTajimasD
dev.off()



TajimasDHilo_plot<-TajimasDHilo %>%
  na.omit() %>%
  ggplot(aes(x=BIN_START, y=TajimaD))+
  geom_jitter(size=2.5,alpha=0.7, colour='lightgrey')+
  geom_smooth(method='loess', alpha=0.6, span=0.1, se=F, colour='forestgreen')+
  labs(y="Hawai'i")+
  facet_wrap(~ CHROM, scales = "free_x", nrow = 1, strip.position = "bottom")+
  labs(x = "Scaffold (MB)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=13,face="bold"))+
  theme_classic()+
  theme(legend.position = 'none')+
  theme(axis.text.y=element_text(size=13))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
  #geom_hline(yintercept = 1, linetype='dashed', colour='red', size=1)+
  geom_hline(yintercept = -2, linetype='dashed', colour='red', size=1)+
  geom_vline(
    data = sharedsweep_midpoints,
    aes(xintercept = chr_position),
    colour = "grey50",
    alpha = 0.6,
    linewidth = 0.5,
    inherit.aes = FALSE
  )

TajimasDKauai_outliers<-TajimasDKauai %>%
  filter(TajimaD < -1) %>%
  mutate(ID = paste(CHROM, "_", BIN_START))

TajimasDOahu_outliers<-TajimasDOahu %>%
  filter(TajimaD < -1) %>%
  mutate(ID = paste(CHROM, "_", BIN_START))

TajimasDHilo_outliers<-TajimasDHilo %>%
  filter(TajimaD < -1) %>%
  mutate(ID = paste(CHROM, "_", BIN_START))

write.csv(TajimasDKauai_outliers, "TajimasDKauai_outliers.csv")
write.csv(TajimasDOahu_outliers, "TajimasDOahu_outliers.csv")
write.csv(TajimasDHilo_outliers, "TajimasDHilo_outliers.csv")

intersect(TajimasDHilo_outliers$ID, TajimasDOahu_outliers$ID, TajimasDKauai_outliers$ID)

patchwork_evidence_selectivesweep<-PBS_threeway_plot/ XPEHH_scan_plot/ GenomewideTajDKauai_plot/ TajimasDOahu_plot/ TajimasDHilo_plot

png("GenomewideScans.png", units="in", width=13.0, height=8, res=300)
patchwork_evidence_selectivesweep
dev.off()



##### Tajimas D for Arizona and California. My expectation if there is any signature is that the peaks of differentiation in flies in Hawaii will show signatures of elevated D in Mainland pops.


### CALIFORNIA

TajimasDCalifornia<-read_table("California.TajimasD")

TajimasDCalifornia$CHROM<-gsub("chr", "", TajimasDCalifornia$CHROM)

TajimasDCalifornia<-TajimasDCalifornia %>%
  filter(N_SNPS > 100)

TajimasDCalifornia$TajimaD<-as.double(TajimasDCalifornia$TajimaD)

mean(TajimasDCalifornia$TajimaD)

HistogramTajimasD_California<-TajimasDCalifornia %>%
  ggplot(aes(x=TajimaD))+
  geom_histogram(fill='grey', colour='black')+
  theme_pubr()+
  labs(x=expression(paste("Tajima's D - Hilo", italic(" D"))),
       y="count (10kb windows)")+
  geom_vline(xintercept = 0, linetype='dashed', colour='black')

TajimasDCalifornia %>%
  na.omit() %>%
  ggplot(aes(x=BIN_START, y=TajimaD))+
  geom_jitter(size=2.5,alpha=0.7, colour='lightgrey')+
  geom_smooth(method='loess', alpha=0.6, span=0.1, se=F, colour='orange')+
  labs(y=expression(paste("California - Tajima's", italic(" D"))))+
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Scaffold (MB)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=13,face="bold"))+
  theme_classic()+
  theme(legend.position = 'none')+
  theme(axis.text.y=element_text(size=13))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
  geom_hline(yintercept = 2, linetype='dashed', colour='red', size=1)+
  geom_hline(yintercept = -2, linetype='dashed', colour='red', size=1)



### ARIZONA samples.
TajimasDArizona<-read_table("Arizona.TajimasD")

TajimasDArizona$CHROM<-gsub("chr", "", TajimasDArizona$CHROM)

TajimasDArizona<-TajimasDArizona %>%
  filter(N_SNPS > 100)

TajimasDArizona$TajimaD<-as.double(TajimasDArizona$TajimaD)

mean(TajimasDArizona$TajimaD)

HistogramTajimasDArizona<-TajimasDArizona %>%
  ggplot(aes(x=TajimaD))+
  geom_histogram(fill='grey', colour='black')+
  theme_pubr()+
  labs(x=expression(paste("Tajima's D - Hilo", italic(" D"))),
       y="count (10kb windows)")+
  geom_vline(xintercept = 0, linetype='dashed', colour='black')

TajimasDArizona %>%
  na.omit() %>%
  ggplot(aes(x=BIN_START, y=TajimaD))+
  geom_jitter(size=2.5,alpha=0.7, colour='lightgrey')+
  geom_smooth(method='loess', alpha=0.6, span=0.1, se=F, colour='orange')+
  labs(y=expression(paste("California - Tajima's", italic(" D"))))+
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Scaffold (MB)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=13,face="bold"))+
  theme_classic()+
  theme(legend.position = 'none')+
  theme(axis.text.y=element_text(size=13))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
  geom_hline(yintercept = 2, linetype='dashed', colour='red', size=1)+
  geom_hline(yintercept = -2, linetype='dashed', colour='red', size=1)



# Load shared sweep BED file (regions significant across multiple tests)
sharedsweeps <- read.table("EvidenceOfSelectiveSweeps.SharedAcrossTests.bed")
colnames(sharedsweeps) <- c("chr", "start", "end", "shared_across", "datasets", "X1", "X2", "X3")

# Convert to data.tables for fast interval joining
taj          <- as.data.table(TajimasDAllCombined)
sharedsweeps <- as.data.table(sharedsweeps)
sharedsweeps$chr <- as.character(sharedsweeps$chr)

# Add window end position (10kb windows) for interval overlap
taj[, BIN_END := BIN_START + 10000 - 1]

# Set keys for interval join
setkey(taj, CHROM, BIN_START, BIN_END)
setkey(sharedsweeps, chr, start, end)

# Perform interval overlap join: find Tajima's D windows overlapping sweep regions
overlaps <- foverlaps(
  taj, sharedsweeps,
  by.x    = c("CHROM", "BIN_START", "BIN_END"),
  by.y    = c("chr", "start", "end"),
  type    = "any",
  nomatch = 0L
)

# Label each Tajima's D window as sweep outlier or background
taj[, outlier := "no"]
taj[unique(overlaps, by = c("CHROM", "BIN_START", "BIN_END")), outlier := "yes"]

TajimasDAllCombined <- taj

# ==============================================================================
# SECTION 9: TAJIMA'S D SUMMARY PLOT — SWEEPS vs BACKGROUND
# ==============================================================================
# Pointrange plot: compare mean Tajima's D in sweep vs background windows
# per population source

# Colour palette for populations
my_palette <- c(
  "Kauai"      = "tomato1",
  "Oahu"       = "deepskyblue",
  "Hilo"       = "forestgreen",
  "Arizona"    = "orange",
  "California" = "lightpink1"
)

# Assign colour group: grey for background, population colour for sweep outliers
TajimasDAllCombined <- TajimasDAllCombined %>%
  mutate(colour_group = ifelse(outlier == "no", "grey", as.character(Source)))

# Pointrange plot: mean ± 1 SD Tajima's D per source, split by outlier status
TajimasD_PBSplot <- ggplot(TajimasDAllCombined, aes(x = Source, y = TajimaD, colour = colour_group)) +
  stat_summary(
    fun.data = mean_sdl, mult = 1,      # mean ± 1 SD
    geom     = "pointrange",
    position = position_dodge(width = 0.6),
    size     = 1
  ) +
  scale_colour_manual(values = c("grey", my_palette)) +
  theme_classic(base_size = 14) +
  labs(colour = "Shared sweep outliers", y = "Tajima's D", x = "") +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey50') +
  theme(
    axis.title      = element_text(size = 17),
    axis.text       = element_text(size = 17),
    legend.title    = element_text(size = 17),
    legend.text     = element_text(size = 17),
    legend.position = 'none'
  ) +
  coord_flip()

png("TajimasDAcrossIslandsForSweepsvsNonSweeps.png", units = "in", width = 5.0, height = 4.5, res = 300)
TajimasD_PBSplot
dev.off()

# ==============================================================================
# SECTION 10: STATISTICAL TEST — TAJIMA'S D SWEEP vs BACKGROUND
# ==============================================================================
# Wilcoxon rank-sum test (non-parametric) comparing Tajima's D between
# sweep and background windows, per population source

wilcox_with_effect <- TajimasDAllCombined %>%
  filter(outlier %in% c("yes", "no")) %>%
  group_by(Source) %>%
  wilcox_test(TajimaD ~ outlier) %>%
  add_significance("p") %>%
  mutate(method = "Wilcoxon Rank-Sum")

# Compute effect sizes (rank-biserial correlation)
effect_sizes <- TajimasDAllCombined %>%
  filter(outlier %in% c("yes", "no")) %>%
  group_by(Source) %>%
  wilcox_effsize(TajimaD ~ outlier)

# Combine test results and effect sizes
final_wilcox_results <- left_join(wilcox_with_effect, effect_sizes, by = "Source")

# ==============================================================================
# SECTION 11: GENOME-WIDE TAJIMA'S D PLOTS PER ISLAND
# ==============================================================================
# Plot Tajima's D across chromosomes for each Hawaiian island separately,
# with LOESS smoothing and shared sweep positions marked as vertical lines

# Select three Hawaiian island sources (exclude mainland outgroups)
sources_to_plot <- TajimasDAllCombined %>%
  pull(Source) %>%
  unique() %>%
  setdiff(c("California", "Arizona")) %>%
  mixedsort() %>%
  head(3)

# Compute shared sweep midpoints for vertical line overlay
sharedsweep_midpoints <- sharedsweeps %>%
  filter(!(chr %in% c("c1", "5"))) %>%
  mutate(
    CHROM        = factor(chr, levels = mixedsort(unique(chr))),
    chr_position = (start + end) / 2
  )

# Shared theme elements for Tajima's D genome plots
tajima_genome_theme <- list(
  facet_wrap(~CHROM, scales = "free_x", nrow = 1, strip.position = "bottom"),
  labs(x = "Scaffold (MB)"),
  theme_classic(),
  theme(
    legend.position = 'none',
    axis.text.y     = element_text(size = 13),
    axis.text.x     = element_text(size = 13),
    axis.ticks.x    = element_blank()
  ),
  geom_hline(yintercept =  0, linetype = 'dashed', colour = 'black', size = 1),   # Neutral expectation
  geom_hline(yintercept = -2, linetype = 'dashed', colour = 'red',   size = 1),   # Strong negative selection threshold
  geom_vline(                                                                       # Shared sweep positions
    data      = sharedsweep_midpoints,
    aes(xintercept = chr_position),
    colour    = "grey50",
    alpha     = 0.6,
    linewidth = 0.5,
    inherit.aes = FALSE
  )
)

# Hilo (Hawai'i island) — green
TajimasDHilo_plot <- TajimasDHilo %>%
  na.omit() %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'lightgrey') +
  geom_smooth(method = 'loess', alpha = 0.6, span = 0.1, se = F, colour = 'forestgreen') +
  labs(y = "Hawai'i") +
  tajima_genome_theme

# Kauai — red
TajimasDKauai_plot <- TajimasDKauai %>%
  na.omit() %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'lightgrey') +
  geom_smooth(method = 'loess', alpha = 0.6, span = 0.1, se = F, colour = 'tomato1') +
  labs(y = "Kauai") +
  tajima_genome_theme

# Oahu — blue
TajimasDOahu_plot <- TajimasDOahu %>%
  na.omit() %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'lightgrey') +
  geom_smooth(method = 'loess', alpha = 0.6, span = 0.1, se = F, colour = 'deepskyblue') +
  labs(y = "Oahu") +
  tajima_genome_theme

# Combine all three island Tajima's D plots
ggarrange(TajimasDKauai_plot, TajimasDOahu_plot, TajimasDHilo_plot, nrow = 3, ncol = 1)
