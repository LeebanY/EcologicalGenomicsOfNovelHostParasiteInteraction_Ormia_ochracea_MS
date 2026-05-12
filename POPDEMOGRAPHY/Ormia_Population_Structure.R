# =============================================================================
# Ormia ochracea Population Structure Analysis
# Investigates whether Ormia populations form a single metapopulation or show
# genetic structure by geography (island / US state) and host.
# Analyses: PCA, FST (windowed), Tajima's D, sequencing depth QC, ADMIXTURE
# =============================================================================

# --- Load required libraries -------------------------------------------------
library(tidyverse)    # Data manipulation and ggplot2
library(viridis)      # Colour scales
library(ggpubr)       # Publication-ready plot helpers (loaded twice in original)
library(data.table)   # Fast data import
library(ggExtra)      # Marginal plots
library(RColorBrewer) # Colour palettes
library(cowplot)      # Plot composition

setwd("~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC")

# =============================================================================
# SECTION 1: PCA — ALL POPULATIONS (all batches combined)
# Input: LD-pruned SNP set; eigenvec/eigenval files produced by PLINK
# =============================================================================

# --- Load PCA output ---------------------------------------------------------
# .eigenvec: sample scores on each PC (one row per individual)
# .eigenval: eigenvalues used to compute % variance explained per PC
pca_all    <- read.table("PCA_pruned_input_allbatches.eigenvec", header = F)
eigenval_all <- scan("PCA_pruned_input_allbatches.eigenval")

# --- Load population metadata ------------------------------------------------
# Contains sample IDs, population labels, and collection year
pops <- read_tsv("populations_file_additional.txt",
                 col_names = c("Samples", "Population", "NA", "Year"))
# Remove the literal "NA" column that arises from the 3-field file format
pops <- pops[, !is.na(names(pops)) & names(pops) != "NA"]

# --- Clean up eigenvec data --------------------------------------------------
# PLINK outputs a redundant first column (family ID = individual ID); remove it
pca_all <- pca_all[, -1]
names(pca_all)[1]              <- "individuals"
names(pca_all)[2:ncol(pca_all)] <- paste0("PC", 1:(ncol(pca_all) - 1))

glimpse(pca_all$individuals)  # Sanity check: confirm individual IDs look correct

# --- Match metadata to PCA scores --------------------------------------------
# Look up the population label for each individual using their sample ID
pca_all$Sample_ID <- pops$Population[match(pca_all$individuals, pops$Samples)]

# Parse the population string (format: LOCATION_..._SITE) to extract fields
pca_all$Location <- str_extract(pca_all$Sample_ID, "^[^_]+")   # Everything before first "_"
pca_all$Site     <- str_extract(pca_all$Sample_ID, "[^_]+$")   # Everything after last "_"
pca_all$Year     <- pops$Year[match(pca_all$individuals, pops$Samples)]

# Assign broad geographic groups for colouring in plots
pca_all <- pca_all %>%
  mutate(
    Groups = case_when(
      Location %in% c("YAVAPAI", "COCHISE", "SANTACRUZ")          ~ "Arizona",
      Location %in% c("SANTABARBARA", "VENTURA", "LOSANGELES")     ~ "California",
      Location == "KAUAI"   ~ "Kaua'i",
      Location == "OAHU"    ~ "O'ahu",
      Location == "MOLOKAI" ~ "Moloka'i",
      Location == "HILO"    ~ "Hawai'i"
    )
  )

# --- Compute % variance explained (PVE) per PC -------------------------------
pve_all <- data.frame(
  PC      = 1:10,
  pve_all = eigenval_all / sum(eigenval_all) * 100
)

# Scree plot (optional diagnostic; stored in 'a' but not saved here)
a <- ggplot(pve_all, aes(PC, pve_all)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()

# --- Define colour palette for geographic groups -----------------------------
my_palette <- c(
  "Kaua'i"     = "tomato1",
  "O'ahu"      = "deepskyblue",
  "Moloka'i"   = "mediumorchid1",
  "Hawai'i"    = "forestgreen",
  "Arizona"    = "orange",
  "California" = "lightpink1"
)

# Set factor levels to control legend and colour order
pca_all$Groups <- factor(pca_all$Groups, levels = names(my_palette))

# --- Build PC1 vs PC2 scatter plot -------------------------------------------
# Colour = geographic group; shape = collection year
b <- ggplot(pca_all, aes(PC1, PC2, col = Groups, shape = as.factor(Year))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "darkgrey", alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey", alpha = 1) +
  geom_point(size = 3, alpha = 1)

b <- b + coord_equal() + theme_bw() + theme(legend.position = "none")
b <- b + scale_colour_manual(values = my_palette)
# Add axis labels showing PVE for PC1 and PC2
b <- b +
  xlab(paste0("PC1 (", signif(pve_all$pve_all[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve_all$pve_all[2], 3), "%)"))

# Save to PNG (ratio = 0.5 gives a wider than tall aspect)
png("PCA_PC1_PC2_Ormia_Filtered_ByIsland.png",
    units = "in", width = 10.0, height = 5.0, res = 300)
b + coord_fixed(ratio = 0.5)
dev.off()


# =============================================================================
# SECTION 2: PCA — HAWAIIAN POPULATIONS ONLY
# Separate PCA run restricted to Hawaiian islands to resolve finer structure
# =============================================================================

pca_hw     <- read.table("PCA_Feb2025/HawaiiOrmia_pruned_input.eigenvec", header = F)
eigenval_hw  <- scan("PCA_Feb2025/HawaiiOrmia_pruned_input.eigenval")

# Reload the simpler populations file (no Year column) for the Hawaii subset
pops <- read_tsv("populations_file.txt", col_names = c("Samples", "Population"))

# Clean eigenvec
pca_hw <- pca_hw[, -1]
names(pca_hw)[1]               <- "individuals"
names(pca_hw)[2:ncol(pca_hw)] <- paste0("PC", 1:(ncol(pca_hw) - 1))

glimpse(pca_hw$individuals)

# Match metadata
pca_hw$Sample_ID <- pops$Population[match(pca_hw$individuals, pops$Samples)]
pca_hw$Location  <- str_extract(pca_hw$Sample_ID, "^[^_]+")
pca_hw$Site      <- str_extract(pca_hw$Sample_ID, "[^_]+$")

pca_hw <- pca_hw %>%
  mutate(
    Groups = case_when(
      Location %in% c("YAVAPAI", "COCHISE", "SANTACRUZ")       ~ "Arizona",
      Location %in% c("SANTABARBARA", "VENTURA", "LOSANGELES")  ~ "California",
      Location == "KAUAI"   ~ "Kaua'i",
      Location == "OAHU"    ~ "O'ahu",
      Location == "MOLOKAI" ~ "Moloka'i",
      Location == "HILO"    ~ "Hawai'i"
    )
  )

# For individuals where Site couldn't be parsed (e.g. single-part ID), use Location
pca_hw$Site[is.na(pca_hw$Site)] <- pca_hw$Location[is.na(pca_hw$Site)]

# PVE for Hawaii PCA
pve_hw <- data.frame(
  PC     = 1:10,
  pve_hw = eigenval_hw / sum(eigenval_hw) * 100
)

# Scree plot (diagnostic only)
a <- ggplot(pve_hw, aes(PC, pve_hw)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()

# NOTE: my_palette is redefined here to only include Hawaiian islands
# (the earlier 6-colour palette is overwritten)
my_palette <- c(
  "Kaua'i" = "tomato1",
  "O'ahu"  = "deepskyblue",
  "Hilo"   = "forestgreen"
)

pca_hw$Groups <- factor(pca_hw$Groups, levels = names(my_palette))
pca_hw$Site   <- factor(pca_hw$Site)   # Ensure Site is a factor for shape mapping

# Build Hawaii-specific PCA plot
# Shape encodes collection site within each island; faceted by island (Location)
b_hw <- ggplot(pca_hw, aes(PC1, PC2, col = Groups, shape = Site)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "lightgrey", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "lightgrey", alpha = 0.5) +
  geom_jitter(size = 3, alpha = 1, width = 0.1) +  # Small jitter to separate overlapping points
  scale_shape_manual(values = 0:17)                 # Up to 18 distinct shapes for sites

b_hw <- b_hw +
  coord_equal() +
  theme(legend.position = "none") +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")
  )

# Facet by island location (4 rows: Kauai, Oahu, Molokai, Hilo)
b_hw <- b_hw +
  scale_colour_manual(values = my_palette) +
  facet_wrap(~Location, nrow = 4)

b_hw <- b_hw +
  xlab(paste0("PC1 (", signif(pve_hw$pve_hw[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve_hw$pve_hw[2], 3), "%)"))

# Remove facet strip labels (islands identified by colour instead)
b_hw <- b_hw +
  theme(strip.text = element_blank()) +
  scale_y_continuous(position = "left")

png("PCA_PC1_PC2_Ormia_Filtered_HAWAIISPECIFIC.png",
    units = "in", width = 6.0, height = 10.0, res = 300)
b_hw
dev.off()

# --- Export Hawaiian PC scores as covariates for downstream modelling --------
# Subset the all-populations PCA to Hawaiian islands only; save for use as
# population structure covariates (e.g. in GEMMA or BayesR)
pca_hw <- pca_all %>%
  filter(Location %in% c("KAUAI", "OAHU", "MOLOKAI", "HILO"))

write.csv(pca_hw, "pca_hawaii_covariates_for_modelling.csv",
          quote = F, row.names = F)


# =============================================================================
# SECTION 3: PCA — ARIZONA POPULATIONS
# Question: do flies cluster by geographic location or by host?
# =============================================================================

pca_Arizona    <- read.table("PCA_Arizona/Arizonan_pca_results.eigenvec", header = FALSE)
eigenval_Arizona <- scan("PCA_Arizona/Arizonan_pca_results.eigenval")

pca_Arizona <- pca_Arizona[, -1]
colnames(pca_Arizona)[1]                  <- "individuals"
colnames(pca_Arizona)[2:ncol(pca_Arizona)] <- paste0("PC", seq_len(ncol(pca_Arizona) - 1))

# Match metadata (pops file last loaded is the 2-column version)
pca_Arizona$Sample_ID <- pops$Population[match(pca_Arizona$individuals, pops$Samples)]
pca_Arizona$Location  <- str_extract(pca_Arizona$Sample_ID, "^[^_]+")
pca_Arizona$Host      <- str_extract(pca_Arizona$Sample_ID, "[^_]+$")  # Host cricket species

pve_Arizona <- data.frame(
  PC  = seq_along(eigenval_Arizona),
  PVE = eigenval_Arizona / sum(eigenval_Arizona) * 100
)

# Optional scree plot
pve_plot_Arizona <- ggplot(pve_Arizona, aes(PC, PVE)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()

arizona_palette <- c(
  "COCHISE"   = "gold3",
  "YAVAPAI"   = "saddlebrown",
  "SANTACRUZ" = "indianred"
)

# Colour = county (location); shape = host cricket species
ArizonanPCAflies <- ggplot(
  pca_Arizona,
  aes(PC1, PC2, colour = Location, shape = Host)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
  geom_point(size = 4, alpha = 0.8) +
  coord_equal() +
  theme_bw() +
  scale_colour_manual(values = arizona_palette) +
  xlab(paste0("PC1 (", signif(pve_Arizona$PVE[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve_Arizona$PVE[2], 3), "%)")) +
  labs(title = "Arizona samples")


# =============================================================================
# SECTION 4: PCA — CALIFORNIA POPULATIONS
# Same structure as Arizona section above
# =============================================================================

pca_California    <- read.table("California_PCA/Californian_pca_results.eigenvec", header = FALSE)
eigenval_California <- scan("California_PCA/Californian_pca_results.eigenval")

pca_California <- pca_California[, -1]
colnames(pca_California)[1]                     <- "individuals"
colnames(pca_California)[2:ncol(pca_California)] <- paste0("PC", seq_len(ncol(pca_California) - 1))

pca_California$Sample_ID <- pops$Population[match(pca_California$individuals, pops$Samples)]
pca_California$Location  <- str_extract(pca_California$Sample_ID, "^[^_]+")
pca_California$Host      <- str_extract(pca_California$Sample_ID, "[^_]+$")

pve_California <- data.frame(
  PC  = seq_along(eigenval_California),
  PVE = eigenval_California / sum(eigenval_California) * 100
)

pve_plot_California <- ggplot(pve_California, aes(PC, PVE)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()

california_palette <- c(
  "VENTURA"      = "deeppink2",
  "SANTABARBARA" = "purple3",
  "LOSANGELES"   = "plum"
)

CaliforniaPCAflies <- ggplot(
  pca_California,
  aes(PC1, PC2, colour = Location, shape = Host)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
  geom_point(size = 4, alpha = 0.8) +
  coord_equal() +
  theme_bw() +
  scale_colour_manual(values = california_palette) +
  xlab(paste0("PC1 (", signif(pve_California$PVE[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve_California$PVE[2], 3), "%)")) +
  labs(title = "California samples")

# Interpretation note (from original):
# Little evidence of host-driven structure in either state.
# Slight geographic signal: LA and Ventura flies are similar; Santa Barbara is more distinct.
# Arizona shows a possible east-west gradient (Cochise → Yavapai → Santa Cruz) but weak.

# --- Save individual and combined PCA plots ----------------------------------
png("PCA_PC1_PC2_CaliforniaPCAflies.png", units = "in", width = 8.0, height = 6.0, res = 300)
CaliforniaPCAflies
dev.off()

png("PCA_PC1_PC2_ArizonaPCAflies.png", units = "in", width = 8.0, height = 6.0, res = 300)
ArizonanPCAflies
dev.off()

# Side-by-side panel figure for California and Arizona
png("PCA_PC1_PC2_ArizonaCaliforniaPCAflies.png", units = "in", width = 12.0, height = 6.0, res = 300)
ggarrange(
  CaliforniaPCAflies, ArizonanPCAflies,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = "right"
)
dev.off()


# =============================================================================
# SECTION 5: WINDOWED FST — PAIRWISE ISLAND COMPARISONS
# Preliminary FST estimated by VCFtools (Weir & Cockerham) in sliding windows.
# Three comparisons: Kauai-Oahu, Hilo-Oahu, Hilo-Kauai
# =============================================================================

# Helper: rename raw scaffold IDs (OY7271xx.1) to readable chromosome labels.
# This block is repeated for each comparison dataset — a function would reduce
# repetition but the inline approach is retained here for transparency.
rename_chroms <- function(df) {
  df$CHROM <- substring(df$CHROM, 14)  # Strip leading assembly prefix characters
  df$CHROM <- gsub("OY727105.1", "1",  df$CHROM)
  df$CHROM <- gsub("OY727106.1", "2",  df$CHROM)
  df$CHROM <- gsub("OY727107.1", "3",  df$CHROM)
  df$CHROM <- gsub("OY727108.1", "4",  df$CHROM)
  df$CHROM <- gsub("OY727109.1", "5",  df$CHROM)
  df$CHROM <- gsub("OY727110.1", "X",  df$CHROM)
  df$CHROM <- gsub("OY727111.1", "MT", df$CHROM)
  df
}

# --- 5a. Kauai vs Oahu -------------------------------------------------------
FST_ormia_kauai_oahu <- read_table(
  "~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/Kauai_oahu.fst.windowed.weir.fst"
)

# Floor negative FST values at 0 (artefacts of small sample sizes in VCFtools)
# and require ≥100 variants per window for reliable estimates
FST_ormia_kauai_oahu <- FST_ormia_kauai_oahu %>%
  mutate(WEIGHTED_FST = ifelse(WEIGHTED_FST < 0, 0, WEIGHTED_FST)) %>%
  filter(N_VARIANTS > 100)

mean(FST_ormia_kauai_oahu$WEIGHTED_FST)  # Quick check of genome-wide mean FST

ggplot(FST_ormia_kauai_oahu, aes(x = WEIGHTED_FST)) + geom_histogram()  # Distribution check

FST_ormia_kauai_oahu <- rename_chroms(FST_ormia_kauai_oahu)

# Manhattan-style FST plot across chromosomes
KAUAI_OAHU <- FST_ormia_kauai_oahu %>%
  dplyr::filter(WEIGHTED_FST > 0) %>%                      # Remove zero-FST windows
  mutate(POS = (BIN_START + BIN_END) / 2) %>%              # Window midpoint as x position
  ggplot(aes(x = POS / 1e6, y = WEIGHTED_FST)) +           # Position in Mb
  geom_jitter(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.1, se = FALSE) +  # Local trend
  labs(y = "FST (Kauai - Oahu)", x = "Scaffold (MB)") +
  facet_grid(cols = vars(CHROM),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 13, face = "bold"),
        axis.text.y  = element_text(size = 12),
        strip.text.x = element_text(size = 14))

# --- 5b. Hilo vs Oahu --------------------------------------------------------
FST_ormia_hilo_oahu <- read_table(
  "~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/Hilo_oahu.fst.windowed.weir.fst"
)

FST_ormia_hilo_oahu <- FST_ormia_hilo_oahu %>%
  mutate(WEIGHTED_FST = ifelse(WEIGHTED_FST < 0, 0, WEIGHTED_FST)) %>%
  filter(N_VARIANTS > 100)

mean(FST_ormia_hilo_oahu$WEIGHTED_FST)
ggplot(FST_ormia_hilo_oahu, aes(x = WEIGHTED_FST)) + geom_histogram()

FST_ormia_hilo_oahu <- rename_chroms(FST_ormia_hilo_oahu)

HILO_OAHU <- FST_ormia_hilo_oahu %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>%
  ggplot(aes(x = POS / 1e6, y = WEIGHTED_FST)) +
  geom_jitter(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.2, se = FALSE) +
  labs(y = "FST (Hilo - Oahu)", x = "Scaffold (MB)") +
  facet_grid(cols = vars(CHROM),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 13, face = "bold"),
        axis.text.y  = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# --- 5c. Hilo vs Kauai -------------------------------------------------------
FST_ormia_hilo_kauai <- read_table(
  "~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/Kauai_Hilo.fst.windowed.weir.fst"
)

FST_ormia_hilo_kauai <- FST_ormia_hilo_kauai %>%
  mutate(WEIGHTED_FST = ifelse(WEIGHTED_FST < 0, 0, WEIGHTED_FST)) %>%
  filter(N_VARIANTS > 100)

mean(FST_ormia_hilo_kauai$WEIGHTED_FST)
ggplot(FST_ormia_hilo_kauai, aes(x = WEIGHTED_FST)) + geom_histogram()

FST_ormia_hilo_kauai <- rename_chroms(FST_ormia_hilo_kauai)

HILO_KAUAI <- FST_ormia_hilo_kauai %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>%
  ggplot(aes(x = POS / 1e6, y = WEIGHTED_FST)) +
  geom_jitter(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.1, se = FALSE) +
  labs(y = "FST (Hilo - Kauai)", x = "Scaffold (MB)") +
  facet_grid(cols = vars(CHROM),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 13, face = "bold"),
        axis.text.y  = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0, 1))

# Combined FST panel
ggarrange(KAUAI_OAHU, HILO_KAUAI, HILO_OAHU)

# --- Export high-FST outlier windows (Kauai–Oahu) ----------------------------
# Windows with FST > 0.5 are strong candidates for locally adaptive loci
FST_outliers_kauai <- FST_ormia_kauai_oahu %>%
  filter(WEIGHTED_FST > 0.5)

write_delim(FST_outliers_kauai, "FST_kauai_bigoutlier.tsv", delim = "\t")


# =============================================================================
# SECTION 6: TAJIMA'S D — KAUAI (and multi-population comparisons)
# Tajima's D distinguishes selection signatures from neutral evolution:
#   D << 0 → excess rare alleles → positive selection / population expansion
#   D >> 0 → excess intermediate alleles → balancing selection / bottleneck
# Computed by VCFtools in 10kb windows
# =============================================================================

TajimasDKauai <- read.table(
  "~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/Kauai_TajimasD.txt.Tajima.D",
  header = TRUE
)

# Rename scaffold IDs to readable chromosome labels (same mapping as FST)
TajimasDKauai <- rename_chroms(TajimasDKauai)

# Retain only windows with ≥100 SNPs for reliable D estimates
TajimasDKauai <- TajimasDKauai %>%
  filter(N_SNPS > 100)

TajimasDKauai$TajimaD <- as.double(TajimasDKauai$TajimaD)

# --- Histogram of Tajima's D distribution ------------------------------------
# Red dashed lines at ±1 flag windows of interest (|D| > 1 is often used
# as a preliminary outlier threshold, though formal tests differ)
HistogramTajimasD <- TajimasDKauai %>%
  ggplot(aes(x = TajimaD)) +
  geom_histogram(fill = "black", colour = "black") +
  theme_pubr() +
  labs(x = expression(paste("Tajima's", italic(" D"))),
       y = "count (10kb windows)") +
  geom_vline(xintercept =  0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept =  1, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = -1, linetype = "dashed", colour = "red")

png("TajimasDHistogram.png", units = "in", width = 5.0, height = 3.0, res = 300)
HistogramTajimasD
dev.off()

# --- Genome-wide Tajima's D plot for Kauai -----------------------------------
# Autosomes only (X and MT excluded to avoid confounds from sex-linkage /
# maternal inheritance)
GenomewideTajD <- TajimasDKauai %>%
  na.omit() %>%
  filter(!CHROM %in% c("NA", "X", "MT")) %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.5, se = FALSE) +
  labs(y = expression(paste("KAUAI - Tajima's", italic(" D"))),
       x = "Scaffold (MB)") +
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 13),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 13, face = "bold")) +
  geom_hline(yintercept =  0, linetype = "dashed", colour = "black", size = 1) +
  geom_hline(yintercept =  1, linetype = "dashed", colour = "red",   size = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", colour = "red",   size = 1)

png("TajimasDGenomewide.png", units = "in", width = 10.0, height = 4.0, res = 300)
GenomewideTajD
dev.off()

# Fine-resolution version (span = 0.01) — very wiggly; used for local inspection
TajimasDKauai %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.01, se = FALSE) +
  labs(y = "Tajima's D (Kauai)", x = "Scaffold (MB)") +
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text   = element_text(size = 14),
        axis.title  = element_text(size = 13, face = "bold")) +
  geom_hline(yintercept =  0, linetype = "dashed", colour = "black") +
  geom_hline(yintercept =  1, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = -1, linetype = "dashed", colour = "red")

# Export windows with strongly negative Tajima's D (putative selection targets)
TajimasDKauai_outliers <- TajimasDKauai %>%
  filter(TajimaD < -1)

write_tsv(TajimasDKauai_outliers, "TajimasDKauai_outliers.tsv")

# Note: TajimasDPlot object referenced below is not yet defined at this point —
# this png() call will fail unless TajimasDPlot is created earlier in the session
png("~/Desktop/Ochracea_POSTDOC/TajimasD_OrmiaKauai_Preliminary.png",
    units = "in", width = 8.0, height = 5.0, res = 300)
TajimasDPlot
dev.off()

# --- Genes in Tajima's D outlier windows (D < -1) ----------------------------
# These were identified by intersecting outlier windows with gene annotations
intersecting_tajd_genes <- c(
  "PLAAT1", "MAP1B", "GPR171", "NTN1", "dkk3a", "C19orf47", "C2orf66",
  "WDR5", "INPP5B", "ZNF800", "HPRT1", "RESF1", "SLC25A21", "TMEM30A",
  "WDR35", "MBOAT2", "POLR3A", "FTSJ3", "SLC6A2", "ZNF236", "ENPP2",
  "RAD51AP2", "CBLIF", "gcna", "MAF1", "AP3S1", "TGFB2", "C21orf62",
  "PSMD3", "USP14", "TMCO1", "DCTN6", "POLR2C", "KNL1", "gucy2g",
  "UMAD1", "GRHPR", "WDR19"
)

# Map human gene symbols to Drosophila melanogaster orthologs
# (orthologs() presumably from a package such as babelgene or homologene)
tajd_genes_orthologs <- orthologs(
  genes   = intersecting_tajd_genes,
  species = "fruit fly",
  human   = TRUE
)

# Check overlap with a pre-loaded set of fly orthologs (FlyOrthologs object
# must be defined elsewhere in the pipeline)
intersect(tajd_genes_orthologs$symbol, FlyOrthologs$symbol)


# =============================================================================
# SECTION 7: MULTI-POPULATION TAJIMA'S D COMPARISON
# Combines Tajima's D from all populations for a raincloud / boxplot comparison
# =============================================================================

# Add Source labels to each population's Tajima's D dataframe before binding
# (these dataframes — TajimasDOahu, TajimasDHilo, TajimasDCalifornia,
# TajimasDArizona — are loaded elsewhere in the pipeline)
TajimasDOahu    <- TajimasDOahu    %>% mutate(Source = "Oahu")
TajimasDKauai   <- TajimasDKauai   %>% mutate(Source = "Kauai")
TajimasDHilo    <- TajimasDHilo    %>% mutate(Source = "Hilo")
TajimasDCali    <- TajimasDCalifornia %>% mutate(Source = "California")
TajimasDArizona <- TajimasDArizona %>% mutate(Source = "Arizona")

TajimasDAllCombined <- rbind(
  TajimasDOahu, TajimasDKauai, TajimasDHilo, TajimasDCali, TajimasDArizona
)

# Install visualisation packages if not already present
install.packages("ggdist")
install.packages("gghalves")
install.packages("ggforce")
library(ggforce)
library(gghalves)
library(ggdist)   # stat_halfeye: half violin + point interval ("raincloud" style)

# Raincloud + boxplot: shows full distribution shape alongside median/IQR
TajimasDPlot <- ggplot(
  TajimasDAllCombined,
  aes(Source, TajimaD, colour = Source, fill = Source)
) +
  # Half-violin (density) shifted to the right of the axis
  ggdist::stat_halfeye(
    adjust = .5, width = .7, .width = 0,
    justification = -.2, point_colour = NA
  ) +
  # Narrow boxplot overlaid for median and IQR; outliers hidden
  geom_boxplot(width = .25, outlier.shape = NA, colour = "black") +
  theme_pubr() +
  coord_flip() +                                             # Horizontal orientation
  labs(x = "", y = "Tajima's D") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  theme(
    axis.text  = element_text(size = 17),
    axis.title = element_text(size = 17),
    strip.text = element_text(size = 17),
    legend.position = "none"
  ) +
  # Colours match the palette used in the PCA plots
  scale_colour_manual(values = c("orange", "lightpink1", "forestgreen", "tomato1", "deepskyblue")) +
  scale_fill_manual(values   = c("orange", "lightpink1", "forestgreen", "tomato1", "deepskyblue"))

png("~/Desktop/Ochracea_POSTDOC/TajimasD_AltogetherIslandsState.png",
    units = "in", width = 7.0, height = 5.0, res = 300)
TajimasDPlot
dev.off()


# =============================================================================
# SECTION 8: TEMPORAL COMPARISON OF TAJIMA'S D (2022 vs 2024 samples)
# Tests whether Tajima's D has shifted between sampling years at Hilo and Oahu,
# which could reflect rapid evolutionary change or sampling differences
# =============================================================================

# Load 2024 samples for Big Island (Hilo) and Oahu
Hilo24 <- read.table("2024_Bigisland_TajimasD.Tajima.D", header = TRUE)
Oahu24 <- read.table("2024_Oahu_TajimasD.Tajima.D",      header = TRUE)

# Filter: remove windows with NAs and require ≥100 SNPs
Hilo24 <- Hilo24 %>% na.omit() %>% filter(N_SNPS > 100)
Oahu24 <- Oahu24 %>% na.omit() %>% filter(N_SNPS > 100)

# Tag with source population and year
Hilo24 <- Hilo24 %>% mutate(Source = "Hilo", Year = "2024")
Oahu24 <- Oahu24 %>% mutate(Source = "Oahu", Year = "2024")

# Combine with earlier-year data (TajimasDHilo and TajimasDOahu must already
# have a Year column from earlier in the pipeline)
Comparison_Between_Years_TajimasD <- rbind(TajimasDOahu, Oahu24, Hilo24, TajimasDHilo)

# Quick mean comparison
mean(Hilo24$TajimaD)
mean(TajimasDHilo$TajimaD)

# --- Non-parametric tests for year-to-year differences ----------------------
# Wilcoxon rank-sum test: does the Tajima's D distribution differ between years?

Comparison_Between_Years_TajimasD_Hilo <- rbind(Hilo24, TajimasDHilo)
wilcox.test(TajimaD ~ Year, data = Comparison_Between_Years_TajimasD_Hilo)

Comparison_Between_Years_TajimasD_Oahu <- rbind(TajimasDOahu, Oahu24)
wilcox.test(TajimaD ~ Year, data = Comparison_Between_Years_TajimasD_Oahu)

# Visualise temporal shift: faceted by population, x = year
ggplot(Comparison_Between_Years_TajimasD, aes(Year, TajimaD)) +
  ggdist::stat_halfeye(
    adjust = .5, width = .7, .width = 0,
    justification = -.2, point_colour = NA
  ) +
  geom_boxplot(width = .25, outlier.shape = NA) +
  facet_wrap(~Source) +
  theme_pubr() +
  coord_flip() +
  labs(x = "", y = "Tajima's D") +
  theme(axis.text.y  = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")


# =============================================================================
# SECTION 9: SEQUENCING DEPTH QC — ALL SAMPLES
# Mean per-site depth computed by VCFtools (--site-mean-depth)
# Checks for coverage anomalies, especially on the X chromosome
# =============================================================================

MeanDepthAcrossOrmia <- read_table(
  "~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/OrmiaSeqDepth.subsetofsamples.ldepth.mean"
)

MeanDepthAcrossOrmia <- rename_chroms(MeanDepthAcrossOrmia)
glimpse(MeanDepthAcrossOrmia)

# Remove mitochondrial chromosome (often extreme depth due to copy number)
MeanDepthAcrossOrmia <- MeanDepthAcrossOrmia %>% filter(CHROM != "MT")

# Genome-wide depth profile (LOESS smoothed)
Mean_depth_Across_Ormia_Genome <- MeanDepthAcrossOrmia %>%
  ggplot(aes(x = POS, y = MEAN_DEPTH)) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.1, se = FALSE) +
  labs(y = "Mean depth across all samples", x = "Scaffold (MB)") +
  facet_grid(cols = vars(CHROM),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 13, face = "bold"),
        axis.text.y  = element_text(size = 12))

# X chromosome depth in isolation — checking for elevated depth (paralogy,
# mis-mapping) or reduced depth (hemizygosity signal in males)
X_depth <- MeanDepthAcrossOrmia %>%
  filter(CHROM == "X") %>%
  ggplot(aes(x = POS, y = MEAN_DEPTH)) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.1, se = FALSE) +
  labs(y = "Mean depth across all samples", x = "Scaffold (MB)") +
  facet_grid(cols = vars(CHROM),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text   = element_text(size = 14),
        axis.title  = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 12))

X_depth  # Display in console for quick review

png("~/Desktop/Ochracea_POSTDOC/OrmiaMeanDepth.png",
    units = "in", width = 8.0, height = 5.0, res = 300)
Mean_depth_Across_Ormia_Genome
dev.off()

png("~/Desktop/Ochracea_POSTDOC/X_chromosomeDepthOrmia.png",
    units = "in", width = 8.0, height = 5.0, res = 300)
X_depth
dev.off()

# Depth variance per chromosome — identifies chromosomes with unusually variable
# coverage that could indicate mapping artefacts
MeanDepthAcrossOrmia %>%
  ggplot(aes(x = VAR_DEPTH, fill = CHROM)) +
  geom_histogram(colour = "black") +
  facet_wrap(~CHROM, scales = "free_x")

# --- X chromosome depth — Oahu samples only ----------------------------------
# Isolates whether anomalous X depth is sample-wide or Oahu-specific
OahuReadDepth <- read_table(
  "~/Desktop/Ochracea_POSTDOC/Oahu_Site_Depth.subsetsamples.ldepth.mean"
)

OahuReadDepth <- rename_chroms(OahuReadDepth)
OahuReadDepth <- OahuReadDepth %>% filter(CHROM != "MT")

OahuReadDepth %>%
  filter(CHROM == "X") %>%
  ggplot(aes(x = POS, y = MEAN_DEPTH)) +
  geom_smooth(method = "loess", alpha = 0.6, span = 0.1, se = FALSE) +
  labs(y = "Mean depth across all samples", x = "Scaffold (MB)") +
  facet_grid(cols = vars(CHROM),
             space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text   = element_text(size = 14),
        axis.title  = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 12))


