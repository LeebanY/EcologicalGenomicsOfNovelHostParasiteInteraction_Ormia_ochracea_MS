library(tidyverse)    # Core data manipulation and plotting (dplyr, ggplot2, tidyr, etc.)
library(jsonlite)     # JSON parsing
library(rjson)        # Alternative JSON parser (used in the commented-out import section)
library(viridis)      # Colour palettes for plots
library(dunn.test)    # Post-hoc Dunn test following Kruskal-Wallis
source("cheridis.R")  # Custom colour palette functions (scale_fill_cherulean etc.)
library(ggridges)     # Ridge/joy plots for distribution visualisation

# ==============================================================================
# SECTION 1: DATA IMPORT (COMMENTED OUT - already processed and saved as RDS)
# ==============================================================================
# This block loops through RELAX JSON output files, extracts "test results" 
# from each, adds a filename identifier column, and combines into one dataframe.
# The final combined dataframe was saved as:
#   RELAX_full_dataset_060324.csv / .rds
# This section does not need to be re-run unless re-processing raw JSON files.

# ==============================================================================
# SECTION 2: DATA LOADING AND CLEANING
# ==============================================================================

# Load the pre-processed RELAX dataset
RELAX_full_dataset <- readRDS("RELAX_full_dataset_060324.rds")

# Convert rownames (which encode branch/species info) into a proper column
RELAX_full_dataset <- tibble::rownames_to_column(RELAX_full_dataset, "Branches")

# Split the "Branches" column into Species and GeneIdentifier
# Format is: Species_GeneIdentifier
RELAX_full_dataset[c('Species', 'GeneIdentifier')] <- str_split_fixed(RELAX_full_dataset$Branches, '_', 2)

# Clean up GeneIdentifier by removing trailing "_1" suffixes (isoform artifacts)
RELAX_full_dataset[c('GeneIdentifierKeep', 'Garbage')] <- str_split_fixed(RELAX_full_dataset$GeneIdentifier, '_1', 2)

# Select and rename relevant columns, then filter to only the 10 focal species:
# NowFer, OrmOch, PhaObe, TacFer, TacGro, TacLur, TheAcu, CisGlo, EpiSuc, GymRot
RELAX_final_dataset <- RELAX_full_dataset %>%
  dplyr::select(filename, Species, GeneIdentifierKeep, Corrected.p.value, MLE, LRT, p.value) %>%
  dplyr::rename(
    Orthogroup        = filename,
    Species           = Species,
    GeneIdentifier    = GeneIdentifierKeep,
    Corrected.p.value = Corrected.p.value,
    K_parameter       = MLE,      # K = selection intensity parameter from RELAX
    LRT               = LRT,      # Likelihood ratio test statistic
    p.value           = p.value
  ) %>%
  filter(Species == "NowFer" | Species == "OrmOch" | Species == "PhaObe" | Species == "TacFer" |
           Species == "TacGro" | Species == "TacLur" | Species == "TheAcu" | Species == "CisGlo" |
           Species == "EpiSuc" | Species == "GymRot")

# ==============================================================================
# SECTION 3: ANNOTATION MATCHING
# ==============================================================================

# Load functional annotations for all species and Ormia-specific gene IDs
Annotations  <- read_table("FINAL_FULL_ANNOTATIONS_ALLSPECIES.txt")
Ormia_genes  <- read_table("Ormia_GENES_ID.txt")

# Match functional annotations to RELAX dataset by gene identifier
RELAX_final_dataset$FunctionalAnnotation <- Annotations$Functional_annotation[
  match(RELAX_final_dataset$GeneIdentifier, Annotations$ID)]

# Match chromosomal position and gene name from Ormia gene list
RELAX_final_dataset$OrmiaChromosome <- Ormia_genes$Chr[
  match(RELAX_final_dataset$GeneIdentifier, Ormia_genes$GeneID)]

RELAX_final_dataset$Genename <- Ormia_genes$GeneName[
  match(RELAX_final_dataset$GeneIdentifier, Ormia_genes$GeneID)]

# ==============================================================================
# SECTION 4: SUMMARY STATISTICS
# ==============================================================================

# Helper function: format numbers with commas (no scientific notation)
commapos <- function(x, ...) {
  format(abs(x), big.mark = ",", trim = TRUE, scientific = FALSE, ...)
}

# Count number of significantly selected genes per species (Corrected p < 0.05)
RELAX_final_dataset %>%
  group_by(Species) %>%
  filter(Corrected.p.value < 0.05) %>%
  count(Species, sort = TRUE)

# ==============================================================================
# SECTION 5: RELAXED vs INTENSIFIED SELECTION BAR PLOT
# ==============================================================================

# For each species, calculate the percentage of significantly selected genes
# that are relaxed (K < 1) vs intensified (K > 1), then plot as stacked bar chart
IntensifiedvsRelaxedPlotBarPlot <- RELAX_final_dataset_PBS_analysis %>%
  group_by(Species) %>%
  filter(Corrected.p.value < 0.05) %>%
  summarise(
    median_k_selection = median(K_parameter),
    mean_k_selection   = mean(K_parameter),
    count_of_perspp    = n(),
    Relaxed            = (sum(K_parameter < 1) / count_of_perspp) * 100,   # % genes with K < 1
    Intensified        = (sum(K_parameter > 1) / count_of_perspp) * 100    # % genes with K > 1
  ) %>%
  # Pivot to long format so Relaxed/Intensified become a grouping variable
  pivot_longer(cols = 5:6, names_to = "selection_regime_count", values_to = "proportion") %>%
  ggplot(aes(x = Species, y = proportion, fill = selection_regime_count)) +
  geom_col(colour = 'black') +
  labs(
    y      = "Percentage of genes (%) under relaxed v.s. intensified selection",
    x      = "",
    legend = 'Selection Regime'
  ) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.spacing.x = unit(3, "mm")) +
  theme_pubclean() +
  theme(legend.position = 'none') +
  theme(axis.text.y  = element_text(size = 16)) +
  theme(axis.text.x  = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 16))

# Save plot
png("IntensifiedvsRELAXED.png", units = "in", width = 6, height = 4.0, res = 300)
IntensifiedvsRelaxedPlotBarPlot
dev.off()

# ==============================================================================
# SECTION 6: EXPLORATORY K PARAMETER DISTRIBUTION PLOTS
# ==============================================================================

# Histogram of K parameter per species (all genes, not filtered)
ggplot(RELAX_final_dataset, aes(x = K_parameter, colour = Species)) +
  geom_histogram(fill = "white") +
  geom_vline(xintercept = 1, linetype = 'dashed') +  # K=1 = no selection change
  scale_colour_viridis(option = "D", discrete = T) +
  facet_wrap(~Species)

# Boxplot of K parameter per species
ggplot(RELAX_final_dataset, aes(x = Species, y = K_parameter, colour = Species)) +
  geom_boxplot()

# Histogram of K for Ormia only, faceted by chromosome
RELAX_final_dataset %>%
  filter(Species == "OrmOch") %>%
  ggplot(aes(x = K_parameter)) +
  geom_histogram() +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_colour_viridis(option = "D", discrete = T) +
  facet_wrap(~OrmiaChromosome)

# ==============================================================================
# SECTION 7: IDENTIFY ORMIA-SPECIFIC RELAXED / INTENSIFIED GENES
# ==============================================================================
# These functions identify orthogroups where all other species have intensified
# selection (K > 1) but Ormia is relaxed (K < 1), or vice versa.
# This tests whether Ormia shows unique selection regimes compared to relatives.

# Function: flags genes relaxed ONLY in Ormia (all others intensified)
RelaxedinOrmiaOnly <- function(K_parameter, Species) {
  if (sum(K_parameter > 1) == sum(Species != "OrmOch") && sum(Species == "OrmOch") == 1) {
    if_else(Species == "OrmOch",
            ifelse(K_parameter < 1, "ORMIA_RELAX", as.character(K_parameter)),
            as.character(K_parameter))
  } else {
    "INTENSIFIED IN ALL"
  }
}

# Function: flags genes intensified ONLY in Ormia (all others relaxed)
IntensifiedinOrmiaOnly <- function(K_parameter, Species) {
  if (sum(K_parameter < 1) == sum(Species != "OrmOch") && sum(Species == "OrmOch") == 1) {
    if_else(Species == "OrmOch",
            ifelse(K_parameter > 1, "ORMIA_INTENSIFIED", as.character(K_parameter)),
            as.character(K_parameter))
  } else {
    "RELAXED_IN_ALL"
  }
}

# Apply functions to identify Ormia-unique relaxed genes
Ormia_Relaxed_Genes <- RELAX_final_dataset %>%
  filter(Corrected.p.value < 0.05) %>%
  group_by(Orthogroup) %>%
  mutate(Relax_Ormia = RelaxedinOrmiaOnly(K_parameter, Species)) %>%
  filter(Relax_Ormia == "ORMIA_RELAX") %>%
  ungroup()

# Apply functions to identify Ormia-unique intensified genes
Ormia_Intensified_GenesTest <- RELAX_final_dataset %>%
  filter(Corrected.p.value < 0.05) %>%
  group_by(Orthogroup) %>%
  mutate(Intensified_Ormia = IntensifiedinOrmiaOnly(K_parameter, Species)) %>%
  filter(Intensified_Ormia == "ORMIA_INTENSIFIED") %>%
  ungroup()

# Check: how many orthogroups have ALL species with K > 1 and p < 0.05?
# (i.e. intensified selection across the whole group, not Ormia-specific)
RELAX_final_dataset %>%
  group_by(Orthogroup) %>%
  summarize(
    all_K_gt_1   = all(K_parameter > 1),
    all_p_lt_005 = all(Corrected.p.value < 0.05)
  ) %>%
  filter(all_K_gt_1 == TRUE & all_p_lt_005 == TRUE)

# Extract all Ormia genes under intensified selection (not requiring other spp to be relaxed)
Ormia_Intensified_Transcripts <- RELAX_final_dataset %>%
  filter(K_parameter > 1 & Corrected.p.value < 0.05 & Species == "OrmOch")

write.csv(Ormia_Intensified_Transcripts, "Ormia_Intensified_Transcripts_ALL.RELAXANALYSIS.csv")

# ==============================================================================
# SECTION 8: GENOMIC COORDINATES AND MANHATTAN-STYLE PLOT
# ==============================================================================

# Join RELAX results with genomic coordinates
RELAX_final_datasetCoordinates <- left_join(
  Coordinates, RELAX_final_dataset,
  by = join_by(Orthogroups == Orthogroup)
)

# Plot K parameter across the genome for a focal species (CisGlo shown here)
# x-axis = midpoint position of each gene in Mb; faceted by chromosome
RELAX_final_datasetCoordinates %>%
  filter(Corrected.p.value < 0.05) %>%
  filter(Species == "CisGlo") %>%
  mutate(POS = (Start + End) / 2) %>%        # Use gene midpoint as position
  ggplot(aes(x = POS / 1000000, y = K_parameter)) +
  geom_jitter(size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'loess', alpha = 0.5, span = 0.5, se = F) +
  labs(y = "dN/dS") +
  facet_grid(cols  = vars(Chr), space  = "free_x", scales = "free_x", switch = "x") +
  labs(x = "Scaffold (MB)") +
  theme(axis.text  = element_text(size = 14), axis.title = element_text(size = 13, face = "bold")) +
  theme_classic() +
  theme(legend.position = 'none') +
  theme(axis.text.y  = element_text(size = 12)) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  theme(strip.text.x = element_text(size = 14))

# ==============================================================================
# SECTION 9: ODORANT BINDING PROTEIN SUBSET ANALYSIS
# ==============================================================================
# Test whether genes involved in odorant detection show unusual selection patterns
# in Ormia, consistent with sensory adaptation to host cricket signals.

# Get Ormia odorant binding/receptor genes from the sensory gene HOG (hierarchical orthogroup) list
Ormia_OdorantBindingReceptors <- Sensory_genes_HOG %>%
  filter(species == "Ormiaochracea") %>%
  separate_rows(pid)   # Expand rows where multiple protein IDs are in one cell

# Clean up protein IDs to match format in RELAX dataset
Ormia_OdorantBindingReceptors$pid <- gsub('\\.1', '', Ormia_OdorantBindingReceptors$pid)
Ormia_OdorantBindingReceptors$pid <- gsub('ENSEGHP', 'ENSEGHT', Ormia_OdorantBindingReceptors$pid)

# Match odorant receptor genes to their orthogroup in the RELAX dataset
Ormia_OdorantBindingReceptors$TranscriptOrthogroup <- RELAX_final_dataset$Orthogroup[
  match(Ormia_OdorantBindingReceptors$pid, RELAX_final_dataset$GeneIdentifier)]

# Subset RELAX data to only odorant receptor orthogroups
Ormia_OdorantBindingReceptors_Test <- RELAX_final_dataset[
  RELAX_final_dataset$Orthogroup %in% Ormia_OdorantBindingReceptors$TranscriptOrthogroup, ]

# Add gene names to subset
Ormia_OdorantBindingReceptors_Test$GeneNames <- Ormia_OdorantBindingReceptors$GeneNames[
  match(Ormia_OdorantBindingReceptors_Test$Orthogroup, Ormia_OdorantBindingReceptors$TranscriptOrthogroup)]

# Summarise mean K parameter for significantly selected odorant receptor genes per species
Ormia_OdorantBindingReceptors_Test %>%
  filter(Corrected.p.value < 0.05) %>%
  group_by(Species) %>%
  summarise(K_param_mean = mean(K_parameter))

# Bar plot: -log10(corrected p-value) for odorant receptor genes, per orthogroup and species
Ormia_OdorantBindingReceptors_Test %>%
  ggplot(aes(x = Species, y = -log10(Corrected.p.value), fill = Species)) +
  geom_col() +
  scale_fill_viridis(option = 'D', discrete = T) +
  facet_wrap(~as.character(Orthogroup), scales = "free") +
  xlab("Species") +
  ylab("K parameter - Odorant receptors/binding genes") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_flip()

# Density plot: log2(K) for significantly selected odorant receptor genes, per species
Ormia_OdorantBindingReceptors_Test %>%
  filter(Corrected.p.value < 0.05) %>%
  ggplot(aes(x = log2(K_parameter), fill = Species)) +
  scale_fill_cherulean(palette = "spinus") +
  geom_density(alpha = 0.8, colour = 'black') +
  ylab("density") +
  xlab("K parameter - Odorant receptors/binding genes") +
  facet_wrap(~Species) +
  geom_vline(xintercept = log2(1)) +   # K=1 reference line
  geom_rug(sides = "b") +
  theme_bw() +
  theme(legend.position = "none")

# Ridge plot version of density (using ggridges)
Ormia_OdorantBindingReceptors_Test %>%
  ggplot(aes(x = log2(K_parameter), y = Species)) +
  geom_density_ridges(scale = 4) +
  scale_y_discrete(expand  = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme_ridges()

# ==============================================================================
# SECTION 10: GENOME-WIDE K PARAMETER BOXPLOT ACROSS SPECIES
# ==============================================================================

# Boxplot of log(K) for significantly selected genes across all species
png("~/Desktop/Ochracea_POSTDOC/RELAX_across_the_species.png", units = "in", width = 8.0, height = 5.0, res = 300)
RELAX_final_dataset %>%
  group_by(Species) %>%
  filter(Corrected.p.value < 0.05) %>%
  ggplot(aes(x = Species, y = log(K_parameter), fill = Species)) +
  geom_boxplot(outliers = F, width = 0.4, size = 1, colour = 'black') +
  coord_flip() +
  theme_classic() +
  labs(
    y = expression('log(selection intensity parameter' ~ italic(K) ~ ')'),
    x = ""
  ) +
  scale_fill_cherulean(palette = "spinus", discrete = T, reverse = F, alpha = 0.8) +
  theme(axis.text  = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.position = 'none')
dev.off()

# ==============================================================================
# SECTION 11: STATISTICAL TESTING OF K ACROSS SPECIES
# ==============================================================================

# Filter to significantly selected genes for statistical testing
StatTest_RELAX_df <- RELAX_final_dataset %>%
  group_by(Species) %>%
  filter(Corrected.p.value < 0.05)

# Kruskal-Wallis test: is K significantly different across species?
kruskal.test(K_parameter ~ Species, StatTest_RELAX_df)

# Post-hoc Dunn test with Benjamini-Hochberg correction to identify which pairs differ
res <- dunn.test::dunn.test(StatTest_RELAX_df$K_parameter, StatTest_RELAX_df$Species, method = 'bh')
res <- data.frame(res)
res <- res %>% filter(P.adjusted < 0.05)   # Keep only significant pairwise comparisons

# ==============================================================================
# SECTION 12: PBS OUTLIER GENE ANALYSIS
# ==============================================================================
# Compare selection intensity (K) between genes in PBS sweep outlier regions
# vs background genes, to test whether sweep loci show different selection regimes.

# Load PBS outlier gene list and label genes as PBS outlier vs background
PBS_genes_dnds <- read.table("PBS_genes_postprocessing.txt", header = T)

RELAX_final_dataset_PBS_analysis <- RELAX_final_dataset %>%
  mutate(PBS = if_else(Orthogroup %in% PBS_genes_dnds$Orthogroups, "PBS outlier", "background")) %>%
  filter(Corrected.p.value < 0.05)

# Set species factor order for consistent plotting
RELAX_final_dataset_PBS_analysis <- RELAX_final_dataset_PBS_analysis %>%
  mutate(Species = factor(Species, levels = c(
    "PhaObe", "CisGlo", "GymRot", "NowFer", "TacLur",
    "TacFer", "TacGro", "OrmOch", "TheAcu", "EpiSuc")))

# Boxplot: log(K) for PBS outlier vs background genes, per species
RelaxPlotForPBSgenes <- RELAX_final_dataset_PBS_analysis %>%
  group_by(Species) %>%
  ggplot(aes(x = Species, y = log(K_parameter), fill = PBS)) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.5,
    size          = 0.8,
    colour        = 'black',
    position      = position_dodge(width = 0.5)
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    y = expression('log(selection intensity parameter' ~ italic(K) ~ ')'),
    x = ""
  ) +
  scale_fill_manual(values = c("background" = "grey50", "PBS outlier" = "tomato1")) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5, colour = "grey50") +
  theme(
    axis.text      = element_text(size = 17),
    axis.title     = element_text(size = 17),
    legend.position = "none"
  ) +
  ylim(-5, 5)

# Summarise median log(K) per species and PBS group
RELAX_final_dataset_PBS_analysis %>%
  filter(is.finite(log(K_parameter))) %>%
  group_by(Species, PBS) %>%
  summarise(
    median_logK = median(log(K_parameter), na.rm = TRUE),
    n           = n(),
    .groups     = "drop"
  )

# Save PBS plot
png("RELAX_across_the_species_PBSgenes.png", units = "in", width = 4.5, height = 6.0, res = 300)
RelaxPlotForPBSgenes
dev.off()

# Wilcoxon test per species: does K differ between PBS outliers and background genes?
results <- RELAX_final_dataset_PBS_analysis %>%
  group_by(Species) %>%
  summarise(
    t_test  = list(wilcox.test(log(K_parameter + 1e-6) ~ PBS)),   # Small offset avoids log(0)
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    p_adj   = p.adjust(p_value, method = "BH")   # BH correction across species
  )

# Mean K per species and PBS group (untransformed)
Relax_PBS_May2026 <- RELAX_final_dataset_PBS_analysis %>%
  group_by(Species, PBS) %>%
  summarise(
    mean_K  = mean(K_parameter, na.rm = TRUE),
    n       = n(),
    .groups = "drop"
  )

# Mean log(K) per species and PBS group
RELAX_final_dataset_PBS_analysis %>%
  group_by(Species, PBS) %>%
  summarise(
    mean_logK = log(mean(K_parameter, na.rm = TRUE)),
    n         = n(),
    .groups   = "drop"
  )

# ==============================================================================
# SECTION 13: COMPARE K IN ORMIA vs OTHER SPECIES FOR PBS OUTLIER GENES
# ==============================================================================

# For each orthogroup, compute mean K in Ormia and mean K in all other species separately
# This allows direct comparison of Ormia-specific vs background selection intensity
Compare_K_params <- RELAX_final_dataset_PBS_analysis %>%
  group_by(Orthogroup) %>%
  mutate(
    K_non_OrmOch = mean(K_parameter[Species != "OrmOch"], na.rm = TRUE),  # Mean K across non-Ormia spp
    K_OrmOch     = mean(K_parameter[Species == "OrmOch"], na.rm = TRUE)   # Mean K in Ormia
  ) %>%
  ungroup() %>%
  select(Orthogroup, K_non_OrmOch, K_OrmOch, FunctionalAnnotation, Genename, PBS)

# Boxplot: Ormia K parameter for PBS outlier vs background genes
ggplot(Compare_K_params, aes(x = PBS, y = log(K_OrmOch), fill = PBS)) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.5,
    size          = 0.8,
    colour        = 'black',
    position      = position_dodge(width = 0.5)
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    y = expression('log(selection intensity parameter' ~ italic(K) ~ ')'),
    x = ""
  ) +
  scale_fill_manual(values = c("background" = "grey50", "PBS outlier" = "tomato1")) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5, colour = "grey50") +
  theme(
    axis.text       = element_text(size = 17),
    axis.title      = element_text(size = 17),
    legend.position = "none"
  ) +
  ylim(-5, 5)