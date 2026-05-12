# ==============================================================================
# SILENCE GEA (GENOTYPE-ENVIRONMENT ASSOCIATION) ANALYSIS
# Testing for genomic associations with host cricket silence phenotype
# in Hawaiian Ormia ochracea populations
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(tidyverse)   # Data manipulation and plotting
library(ggpubr)      # Publication-ready ggplot themes
library(dplyr)       # Data wrangling
library(purrr)       # Functional programming tools
library(stringr)     # String manipulation
library(LEA)         # Landscape and ecological association analysis (SNMF, LFMM)
library(algatr)      # Landscape genomics toolkit (wraps LFMM)
library(patchwork)   # Combining ggplots
library(viridis)     # Colour palettes
library(vcfR)        # Reading and manipulating VCF files
library(adegenet)    # Population genetics tools
library(gtools)      # mixedsort() for natural chromosome ordering
library(geodata)     # WorldClim environmental data download
library(terra)       # Spatial raster operations
library(babelgene)   # Gene ortholog lookup
library(scales)      # Axis formatting (percent_format etc.)

# ==============================================================================
# SECTION 1: LOAD AND PROCESS PHENOTYPE DATA (PROPORTION SINGING MALES)
# ==============================================================================

# Load field data: counts of male morphs per site
SelectionPressDF <- read_csv("../LFMM_Silence_SelectionPressure.csv", col_names = T)

# Calculate average total male density across all sites (used for normalisation)
avg_male_density <- SelectionPressDF %>%
  summarise(avg_males = mean(Nw_male + CwNw_male + SwNw_male + Fw_male + CwFw_male + SwFw_male))

# Compute proportion of singing (normal-wing) males per site, and a
# density-normalised version to account for variable sampling effort
SelectionPressDF <- SelectionPressDF %>%
  group_by(Site) %>%
  mutate(
    total_males             = Nw_male + CwNw_male + SwNw_male + Fw_male + CwFw_male + SwFw_male,
    Prop_singing            = Nw_male / total_males,
    Normalized_Prop_singing = Prop_singing / avg_male_density$avg_males
  ) %>%
  ungroup()

# ==============================================================================
# SECTION 2: STACKED BAR CHART OF MALE MORPH PROPORTIONS PER SITE
# ==============================================================================

# Reshape to long format, compute within-site proportions, and order by island
df_plot <- SelectionPressDF %>%
  pivot_longer(cols = 4:9, names_to = "Male_morphs", values_to = "Count") %>%
  filter(!is.na(Count)) %>%
  group_by(Island, Site) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  arrange(Island, Site) %>%
  mutate(Site = factor(Site, unique(Site)))   # Preserve island-grouped ordering on axis

# Clean up morph labels by removing "_male" suffix
df_plot$Male_morphs <- gsub("_male", "", df_plot$Male_morphs)

# Stacked proportion bar chart
BarChartShowingProportionsCricketSampling <- ggplot(df_plot, aes(x = Site, y = Proportion, fill = Male_morphs)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "plasma", discrete = T, direction = -1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "", y = "Proportion of Male Morphs", fill = "Male morph") +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1)) +
  coord_flip() +
  theme(
    legend.position  = "top",
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_blank(),
    axis.line.y      = element_blank(),
    axis.text.x      = element_text(size = 18),
    axis.title.x     = element_text(size = 18),
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 16)
  )

# Assign custom colours so that "Nw" (normal-wing/singing) morph gets the
# yellow end of the plasma palette for visual clarity
morph_levels  <- levels(df_plot$Male_morphs)
n_levels      <- length(morph_levels)
plasma_colors <- viridis(n_levels, option = "plasma", direction = 1)
colors_named  <- setNames(plasma_colors, morph_levels)
colors_named["Nw"] <- plasma_colors[1]   # Force Nw to yellow

# Alternative plot using custom colours
ggplot(df_plot, aes(x = Site, y = Proportion, fill = Male_morphs)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_named) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = "", y = "Proportion of Male Morphs", fill = "Male morph") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x     = element_text(size = 18),
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 16),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_blank(),
    axis.line.y      = element_blank(),
    legend.position  = "top"
  ) +
  coord_flip()

png("BarChartShowingProportionsCricketSampling.png", units = "in", width = 5.0, height = 3.5, res = 300)
BarChartShowingProportionsCricketSampling
dev.off()

# ==============================================================================
# SECTION 3: ASSIGN CATEGORICAL SINGING VARIABLE PER SITE
# ==============================================================================
# Sites are grouped into three categories based on proportion of singing males:
#   0   = no singers (AgStation)
#   0.5 = intermediate (mixed sites)
#   1   = fully silent cricket populations (high song attenuation)

SelectionPressDF <- SelectionPressDF %>%
  mutate(Categorical_Singing_Variable = case_when(
    Site %in% c("AgStation")                                                   ~ 0,
    Site %in% c("CommonGround", "Hanalei", "KauaiCommunitycollege",
                "Princeville", "CommunityCenter", "HumaneSociety", "Breadfruit") ~ 0.5,
    Site %in% c("BYU", "PonoKai", "UH", "ChurchLawn")                         ~ 1
  ))

# Alternative categorisation (Princeville moved to high-song group)
SelectionPressDF_ALTERNATE <- SelectionPressDF %>%
  mutate(Categorical_Singing_Variable = case_when(
    Site %in% c("AgStation")                                                    ~ 0,
    Site %in% c("CommonGround", "Hanalei", "KauaiCommunitycollege",
                "CommunityCenter", "HumaneSociety", "Breadfruit")              ~ 0.5,
    Site %in% c("BYU", "PonoKai", "UH", "ChurchLawn", "Princeville")          ~ 1
  ))

write.csv(SelectionPressDF, "SelectionPressDataframe_Silence.csv", quote = F, row.names = F)

# Add k-means cluster assignments (km object assumed to be pre-computed)
SelectionPressDF$cluster <- km$cluster

# Plot proportion singing per site, faceted by categorical singing group
ExampleSilenceGWASPlot <- ggplot(SelectionPressDF_ALTERNATE,
                                 aes(x = reorder(Site, Prop_singing), y = Prop_singing, colour = Island)) +
  geom_point(size = 4) +
  geom_smooth(se = FALSE) +
  theme_classic() +
  facet_wrap(~Categorical_Singing_Variable) +
  labs(x = "Sites", y = "Proportion of Singing Males") +
  theme(
    axis.text      = element_text(size = 17),
    axis.title     = element_text(size = 17),
    axis.text.y    = element_text(size = 17),
    strip.text     = element_text(size = 17),
    legend.position = 'none'
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

png("ExampleSilenceGWASPlot.png", units = "in", width = 12.0, height = 7.0, res = 300)
ExampleSilenceGWASPlot
dev.off()

# ==============================================================================
# SECTION 4: LOAD SAMPLE METADATA AND MATCH PHENOTYPE TO INDIVIDUALS
# ==============================================================================

# Load sample metadata: columns = sample ID, population
MetadataSilence <- read.table("TestForSilenceMetadata.txt", header = F)
colnames(MetadataSilence) <- c("Samples", "Population")

# Match normalised proportion singing to each individual by population
MetadataSilence$Normalized_Prop_singing <- SelectionPressDF$Normalized_Prop_singing[
  match(MetadataSilence$Population, SelectionPressDF$Site)]
MetadataSilence$Normalized_Prop_singing[is.na(MetadataSilence$Normalized_Prop_singing)] <- 0

# Create a binary variable: 1 = no singers, 0 = some singers
MetadataSilence <- MetadataSilence %>%
  mutate(Binary_Var = case_when(
    PropOfSilence == 0 ~ 1,
    PropOfSilence > 0  ~ 0
  ))

# Environmental variable for LFMM (normalised proportion singing)
Silence_Env_Variable <- MetadataSilence %>%
  select(Normalized_Prop_singing)

# Assign categorical singing variable at individual level to match population
MetadataSilence <- MetadataSilence %>%
  mutate(Categorical_Singing_Variable = case_when(
    Population %in% c("AgStation")                                                    ~ 0,
    Population %in% c("CommonGround", "Hanalei", "KauaiCommunitycollege",
                      "Princeville", "CommunityCenter", "HumaneSociety", "Breadfruit") ~ 0.5,
    Population %in% c("BYU", "PonoKai", "UH", "ChurchLawn")                          ~ 1
  ))

# Select the categorical variable for use in LFMM
Silence_Categorical_Variable <- MetadataSilence %>%
  select(Categorical_Singing_Variable)

# ==============================================================================
# SECTION 5: LOAD VCF AND IMPUTE GENOTYPES
# ==============================================================================

# Load autosomal SNPs filtered to Hawaii-only individuals with no missing data
vcf_silent <- read.vcfR("Autosomal_Variants.filtered.norepeats.pruned.Hawaiionly.variants.augemented.NOMISSING.vcf.gz")

# Convert VCF to dosage format (0/1/2 allele count matrix) for LFMM
vcf_silent <- vcf_to_dosage(vcf_silent)
summary(vcf_silent)

# ==============================================================================
# SECTION 6: POPULATION STRUCTURE — SNMF
# ==============================================================================
# Estimate latent factors (K) to control for population structure in LFMM

# Convert genotype matrix to .geno format and run SNMF for K = 1 to 5
geno      <- gen_to_geno(vcf_silent)
LEA::write.geno(geno, "mygeno.geno")
snmf_proj <- LEA::snmf("mygeno.geno", K = 1:5, repetitions = 3, entropy = TRUE)

# Select best K by minimum cross-entropy
ce_values <- sapply(1:5, function(k) min(LEA::cross.entropy(snmf_proj, K = k)))
best_k    <- which.min(ce_values)

# ==============================================================================
# SECTION 7: LFMM — GENOTYPE-ENVIRONMENT ASSOCIATION TEST
# ==============================================================================
# Test for associations between SNP genotypes and the categorical singing variable,
# while controlling for population structure using K = 3 latent factors

ridge_results_silence <- lfmm_run(vcf_silent, Silence_Categorical_Variable, K = 3, lfmm_method = "ridge")

# Inspect top hits
lfmm_table(ridge_results_silence$df, order = TRUE)

# Extract results to a dataframe
ridge_silence_df <- ridge_results_silence$df

# Parse chromosome and position from the SNP identifier string
ridge_silence_df$chromosome    <- str_extract(ridge_silence_df$snp, "^[^_]+")     # Everything before first "_"
ridge_silence_df$chr_position  <- str_extract(ridge_silence_df$snp, "(?<=_).*")   # Everything after first "_"
ridge_silence_df$chromosome    <- gsub("chr", "", ridge_silence_df$chromosome)    # Remove "chr" prefix
ridge_silence_df$chromosome    <- as.numeric(ridge_silence_df$chromosome)
ridge_silence_df$chr_position  <- as.numeric(ridge_silence_df$chr_position)

# Check the 95th percentile of -log10(p) to inform threshold decisions
quantile(-log10(ridge_silence_df$calibrated.pvalue), 0.95, na.rm = T)

# Extract significant outliers (adjusted p < 0.05)
ridge_outliers_silence <- ridge_silence_df %>%
  filter(adjusted.pvalue < 0.05)

write.csv(ridge_outliers_silence, "ridge_outliers_silence.csv")

# ==============================================================================
# SECTION 8: MANHATTAN PLOT — GENOME-WIDE GEA RESULTS
# ==============================================================================

# Significance threshold
threshold_SILENCE <- 2.092974

# Genome-wide Manhattan plot coloured by significance and chromosome
SilencePlot <- ridge_silence_df %>%
  na.omit() %>%
  mutate(
    chromosome  = factor(chromosome, levels = mixedsort(unique(chromosome))),
    neg_log_p   = -log10(adjusted.pvalue),
    point_color = case_when(
      neg_log_p > -log10(0.05)                          ~ "red3",    # Significant hits
      as.numeric(factor(chromosome)) %% 2 == 0          ~ "grey60",  # Alternate chromosomes
      TRUE                                              ~ "grey80"
    )
  ) %>%
  ggplot(aes(x = chr_position, y = neg_log_p, colour = point_color)) +
  geom_jitter(size = 3, alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', colour = 'black') +
  facet_grid(cols = vars(chromosome), space = "free_x", scales = "free_x", switch = "x") +
  scale_colour_identity() +
  labs(x = "Scaffold (MB)", y = expression(-log[10] ~ "(P)")) +
  theme_classic() +
  theme(
    axis.text       = element_text(size = 17),
    axis.title      = element_text(size = 17),
    axis.text.x     = element_text(size = 0),
    axis.ticks.x    = element_blank(),
    strip.text      = element_text(size = 17),
    strip.background = element_blank(),
    legend.position = 'none'
  )

# Combine GEA Manhattan plot with PBS sweep plot (assumed pre-computed)
EvidencePropOfSingingMales <- ggarrange(SilencePlot, PBS_threeway_plot, nrow = 2, ncol = 1)

png("EvidenceOfPropOfSingingMales.png", units = "in", width = 12.0, height = 5, res = 300)
EvidencePropOfSingingMales
dev.off()

# ==============================================================================
# SECTION 9: ANNOTATE MANHATTAN PLOT WITH GENE NAMES
# ==============================================================================

library(org.Dm.eg.db)   # Drosophila gene ontology database
library(ggrepel)        # Non-overlapping text labels

# Load gene annotation positions for significant GEA loci
GeneIDsForPlotting <- read.csv("CategoricalSilenceVariableLFMMgeneAnnotationsForPlotting.csv")

# Prepare annotation dataframe with gene midpoints for plotting
GeneIDsForPlotting_df <- GeneIDsForPlotting %>%
  mutate(Gene = sub("-.*", "", Genename)) %>%    # Trim isoform suffixes
  dplyr::select(chr, start, end, Gene) %>%
  mutate(
    midpoint     = (start + end) / 2,
    chr_position = midpoint / 1e6,               # Convert to Mb if plot uses Mb
    chromosome   = factor(chr, levels = mixedsort(unique(chr)))
  ) %>%
  distinct(chromosome, chr_position, .keep_all = TRUE) %>%
  rename(chromosome = chr)

# Overlay gene positions on Manhattan plot as vertical lines + labels
SilencePlot +
  geom_vline(
    data        = GeneIDsForPlotting_df,
    aes(xintercept = midpoint),
    colour      = "black",
    linetype    = "solid",
    size        = 0.6,
    alpha       = 0.2,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data        = GeneIDsForPlotting_df,
    aes(x = midpoint, y = max(ridge_silence_df$neg_log_p, na.rm = TRUE) + 0.5, label = Gene),
    colour      = "black",
    size        = 3.5,
    segment.size  = 0.2,
    segment.color = "black",
    inherit.aes   = FALSE,
    max.overlaps  = 10
  )

# ==============================================================================
# SECTION 10: ZOOM-IN PLOTS — CHR2 AND CHR4 FOCAL REGIONS
# ==============================================================================
# Combine GEA p-values with PBS/Tajima's D to show convergent signals
# at candidate loci on chromosomes 2 and 4

# --- Rescaling functions for dual-axis plots ---

# Get value ranges for rescaling between GEA p-values and PBS/Tajima's D
range_pvals <- range(-log10(ridge_silence_df$adjusted.pvalue), na.rm = TRUE)
range_PBS   <- range(PBS_threeway$PBSn1A, na.rm = TRUE)

scale_PBS_to_p  <- function(x) (x - range_PBS[1])   / diff(range_PBS)   * diff(range_pvals) + range_pvals[1]
scale_p_to_PBS  <- function(x) (x - range_pvals[1]) / diff(range_pvals) * diff(range_PBS)   + range_PBS[1]

# Add genomic midpoint position to PBS data
PBS_threeway <- PBS_threeway %>%
  mutate(chr_position = (start + end) / 2)

# --- Chr2 zoom: GEA p-values overlaid with PBSn1A sweep signal ---
Chr2_SilentSignificantSNPs <- ggplot() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  # PBS sweep signal as a red line (rescaled to p-value axis)
  geom_line(
    data  = PBS_threeway %>% filter(chromosome == 2),
    aes(x = chr_position, y = scale_PBS_to_p(PBSn1A)),
    color = "red2",
    size  = 1
  ) +
  
  # GEA p-values as points
  geom_jitter(
    data  = ridge_silence_df %>% filter(chromosome == 2),
    aes(x = chr_position, y = -log10(adjusted.pvalue)),
    size  = 2.5, alpha = 1, color = "grey20", shape = 17
  ) +
  
  scale_y_continuous(
    name     = "-log10(p-value)",
    sec.axis = sec_axis(~scale_p_to_PBS(.), name = "PBSn1A")
  ) +
  scale_x_continuous(
    name   = "Genomic position (Mb)",
    limits = c(52000000, 59000000),
    labels = function(x) sprintf("%.1f", x / 1e6)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 14),
    axis.title.y.right = element_text(size = 14),
    legend.position    = "none"
  ) +
  geom_vline(xintercept = 58829212, colour = 'violet') +   # GRIK1 position
  geom_vline(xintercept = 58561647, colour = 'blue')       # TLR3/Toll-7 position

# --- Chr4 zoom: GEA p-values overlaid with Tajima's D ---
Chr4_SilentSignificantSNPs <- ggplot() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  geom_jitter(
    data  = ridge_silence_df %>% filter(chromosome == 4),
    aes(x = chr_position, y = -log10(adjusted.pvalue)),
    size  = 2.5, alpha = 1, color = "grey20", shape = 17
  ) +
  
  # Tajima's D for sweep outlier populations
  geom_point(
    data  = TajimasDAllCombined %>% filter(CHROM == 4 & Source %in% names(my_palette)),
    aes(x = BIN_START, y = scale_tajima_to_p(TajimaD), color = Source),
    size  = 2.5, alpha = 1
  ) +
  scale_color_manual(values = my_palette) +
  
  scale_y_continuous(
    name     = "-log10(p-value)",
    sec.axis = sec_axis(~scale_p_to_tajima(.), name = "Tajima's D")
  ) +
  scale_x_continuous(
    name   = "Genomic position (Mb)",
    limits = c(53400000, 53800000),
    labels = function(x) sprintf("%.1f", x / 1e6)
  ) +
  # Shade the sweep region on Chr4
  annotate("rect", xmin = 53601296, xmax = 53703551,
           ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.2) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 14),
    axis.title.y.right = element_text(size = 14),
    legend.position    = "none"
  )

# Combine Chr2 and Chr4 zoom plots
SignificantSNPs_GeneWidPlot <- ggarrange(Chr2_SilentSignificantSNPs, Chr4_SilentSignificantSNPs)

png("SilenceGEA_GeneWideCHR2.png", units = "in", width = 5, height = 3.5, res = 300)
Chr2_SilentSignificantSNPs
dev.off()

# ==============================================================================
# SECTION 11: GENOTYPE PLOTS FOR SIGNIFICANT GEA HITS
# ==============================================================================
# For each significant SNP, plot the relationship between genotype (0/0, 0/1, 1/1)
# and the proportion of singing males at each population

# Load genotype calls for significant SNPs
Significanthits_Genotypes <- read.table("GenotypesSilenceGxP.txt", header = T)

# Match categorical singing variable to each individual
Significanthits_Genotypes$Categorical_Singing_Variable <- MetadataSilence$Categorical_Singing_Variable[
  match(Significanthits_Genotypes$sample, MetadataSilence$Samples)]

glimpse(Significanthits_Genotypes)

# Create a unique ID for each SNP (chromosome + position)
Significanthits_Genotypes$ID <- paste(Significanthits_Genotypes$Chromosome, Significanthits_Genotypes$pos)

# Function: generate and save one boxplot per SNP showing genotype vs phenotype
plot_genotypes <- function(data, output_dir = "GenotypePlots") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  unique_ids <- unique(data$ID)
  
  for (id in unique_ids) {
    subset_data <- subset(data, ID == id)
    if (nrow(subset_data) < 2) next    # Skip SNPs with too few observations
    
    outfile <- file.path(output_dir, paste0("Genotypes_", id, ".png"))
    
    p <- ggplot(subset_data, aes(x = genotype, y = Categorical_Singing_Variable, colour = genotype)) +
      geom_boxplot(size = 1, width = 0.5) +
      geom_jitter(alpha = 0.3, size = 2, width = 0.05, height = 0.05) +
      scale_colour_manual(values = c("#D55E00", "#0072B2", "#009E73")) +
      theme_classic2() +
      facet_wrap(~Chromosome, nrow = 1,
                 labeller = labeller(.cols = function(x) id)) +
      labs(y = "Proportion of Singing Males", x = "") +
      theme(
        legend.position  = "none",
        panel.border     = element_blank(),
        axis.text        = element_text(size = 14),
        axis.title       = element_text(size = 16),
        strip.background = element_blank(),
        strip.text       = element_blank()
      )
    
    png(outfile, units = "in", width = 4.5, height = 4.0, res = 300)
    print(p)
    dev.off()
  }
  
  # Sanity check: confirm expected number of plots produced
  num_plots <- length(unique_ids)
  if (num_plots != 67) {
    stop(paste("Expected 67 plots, but found", num_plots))
  } else {
    message("All 67 plots created successfully in: ", normalizePath(output_dir))
  }
}

# Run the genotype plotting function
plot_genotypes(Significanthits_Genotypes)


# Sanity check for chr5:20,957,234

# Count how many individuals of each singing category
# occur within each genotype at this SNP position.

Significanthits_Genotypes %>%
  filter(pos == 20957234) %>%
  dplyr::count(genotype, Categorical_Singing_Variable) %>%
  pivot_wider(
    names_from = Categorical_Singing_Variable,
    values_from = n,
    values_fill = 0
  )


Significanthits_Genotypes %>%
  filter(pos == 20957234) %>%
  filter(genotype=="0/0") %>%
  filter(Categorical_Singing_Variable == "0") %>%
  select(sample)
# ==============================================================================
# SECTION 12: NORTH AMERICAN POPULATION GENOTYPES
# ==============================================================================

# Load mainland (North American) genotypes at selected GEA loci
MainlandGenotypes <- read.table("MainlandGenotypesforSelectedAllles.txt", header = T)

# Match normalised proportion singing to each individual
Significanthits_Genotypes$Normalized_Prop_singing <- MetadataSilence$Normalized_Prop_singing[
  match(Significanthits_Genotypes$Sample, MetadataSilence$Samples)]

# ==============================================================================
# SECTION 13: DROSOPHILA GENE EXPRESSION (FLYATLAS2) FOR GEA LOCI
# ==============================================================================
# Use FlyAtlas2 to ask which tissues and developmental stages most highly
# express the fly orthologs of GEA candidate genes

# Load FlyBase gene IDs for GEA outlier genes and FlyAtlas2 expression data
SilenceFlybaseIDs <- read.table("SilenceOutliersGenesFlybaseIDS.txt")
FlyAtlas_DF       <- read.csv("FlyAtlas2_gene_data_2025.csv")

# Clean up IDs to ensure matching
FlyAtlas_DF$X.FBgn       <- trimws(FlyAtlas_DF$X.FBgn)
SilenceFlybaseIDs$V1      <- trimws(SilenceFlybaseIDs$V1)

# Filter FlyAtlas to only GEA candidate genes
Filtered_DF_Flyatlas <- FlyAtlas_DF[FlyAtlas_DF$X.FBgn %in% SilenceFlybaseIDs$V1, ]

# --- Tissue-level summary heatmap (mean expression across all candidate genes) ---

# Remove SD columns, whole body, and carcass; average across candidate genes
TissueStageMeans <- Filtered_DF_Flyatlas %>%
  select(-X.FBgn, -matches("\\.SD|Whole|Carcass")) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

# Pivot to long format and annotate sex/stage
TissueStageMeans_long <- TissueStageMeans %>%
  pivot_longer(cols = everything(), names_to = "Tissue_SexStage", values_to = "MeanExpression") %>%
  na.omit() %>%
  filter(!Tissue_SexStage %in% c("V.Sp.", "M.Sp.", "Accessory.Gland", "Testis", "Ovary")) %>%
  mutate(
    Sex_Stage = case_when(
      grepl("\\.M$", Tissue_SexStage)   ~ "Male",
      grepl("\\.F$", Tissue_SexStage)   ~ "Female",
      grepl("\\.L$", Tissue_SexStage)   ~ "Larval",
      grepl("Ovary", Tissue_SexStage)   ~ "Female",
      grepl("Testis", Tissue_SexStage)  ~ "Male",
      TRUE                              ~ NA_character_
    ),
    Tissue_SexStage = sub("\\..*$", "", Tissue_SexStage)   # Remove sex/stage suffix from tissue name
  )

# Order tissues by overall mean expression for cleaner heatmap
TissueStageMeans_long <- TissueStageMeans_long %>%
  group_by(Tissue_SexStage) %>%
  mutate(MeanOverall = mean(MeanExpression, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Tissue_SexStage = factor(Tissue_SexStage, levels = rev(unique(Tissue_SexStage[order(MeanOverall)]))))

# Heatmap: mean expression per tissue and sex/stage
ggplot(TissueStageMeans_long, aes(x = Sex_Stage, y = Tissue_SexStage, fill = MeanExpression)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("wheat", "orange", "red3"), na.value = "grey90") +
  theme_minimal() +
  labs(x = "", y = "", fill = "Mean Expression (FPKM)") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 15),
    panel.grid   = element_blank(),
    axis.text    = element_text(size = 15)
  )

# --- Per-gene expression heatmap across tissues ---
# 
# Filtered_DF_Flyatlas_AllGenes <- Filtered_DF_Flyatlas %>%
#   select(-X.FBgn, -matches("\\.SD|Whole|Carcass")) %>%
#   pivot_longer(cols = 4:40, names_to = "Tissue", values_to = "FPKM") %>%
#   filter(!Tissue %in% c("V.Sp.", "M.Sp.", "Accessory.Gland", "Testis", "Ovary")) %>%
#   mutate(
#     Sex_Stage = case_when(
#       grepl("\\.M$", Tissue)   ~ "Male",
#       grepl("\\.F$", Tissue)   ~ "Female",
#       grepl("\\.L$", Tissue)   ~ "Larval",
#       grepl("Ovary", Tissue)   ~ "Female",
#       grepl("Testis", Tissue)  ~ "Male",
#       TRUE                    ~ NA_character_
#     ),
#     Tissue = sub("\\..*$", "", Tissue)
#   )
# 
# # Heatmap per gene, faceted by sex/stage
# ggplot(Filtered_DF_Flyatlas_AllGenes, aes(x = Tissue, y = Symbol, fill = FPKM)) +
#   geom_tile(color = "black") +
#   scale_fill_gradientn(colors = c("wheat", "orange", "red3"), na.value = "grey90") +
#   theme_minimal() +
#   labs(x = "", y = "", fill = "Mean Expression (FPKM)") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 15),
#     panel.grid  = element_blank(),
#     axis.text   = element_text(size = 9)
#   ) +
#   facet_wrap(~Sex_Stage)

# ==============================================================================
# SECTION 14: GENE ONTOLOGY ENRICHMENT ANALYSIS
# ==============================================================================
# Test whether GEA candidate genes are enriched for particular GO functions
# relative to the genome-wide gene set, using Fisher's exact test

# Load full genome gene function counts and GEA candidate gene function counts
full_gene_set <- read.table("OrmiaCountsofAllGeneFunctions.tsv", header = T, na.strings = "NA", fill = TRUE)
full_gene_set$Function[full_gene_set$Function == ""] <- NA
full_gene_set <- na.omit(full_gene_set)

genes_of_interest <- read.table("SilenceCalibratedPval0.5.enriched.genes",
                                header = TRUE, na.strings = "NA", fill = TRUE,
                                row.names = NULL, check.names = FALSE)
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

# Clean up function names for display
Significantly_enriched_GO_terms_PosSelection$Function <- gsub(
  "_", " ", Significantly_enriched_GO_terms_PosSelection$Function)

# Top 30 enriched GO terms for plotting
Top30forPlottingPBS <- Significantly_enriched_GO_terms_PosSelection %>%
  top_n(-30, adj_p_value) %>%
  arrange(desc(enrichment))

# ==============================================================================
# SECTION 15: FST BETWEEN COLD AND HOT SPOTS (SINGING vs SILENT POPULATIONS)
# ==============================================================================
# Calculate windowed FST between a site with high silence (PonoKai) and one with
# no silence (AgStation) to identify differentiated genomic regions

# --- PonoKai vs AgStation ---
fst_ag_pono <- read.table("COLDVERSUSHOTSPOT_FST.windowed.weir.fst", header = T)

fst_ag_pono$chromosome <- as.double(gsub("chr", "", fst_ag_pono$CHROM))
fst_ag_pono <- na.omit(fst_ag_pono)

# Set 95th percentile as outlier threshold
ag_pono_threshold <- quantile(fst_ag_pono$WEIGHTED_FST, 0.95, na.rm = T)

# Manhattan plot: FST across genome
FSTPONOKAIAGSTATION <- fst_ag_pono %>%
  filter(N_VARIANTS > 100) %>%
  mutate(
    chromosome  = factor(chromosome, levels = mixedsort(unique(chromosome))),
    pos         = (BIN_START + BIN_END) / 2,
    point_color = case_when(
      WEIGHTED_FST > ag_pono_threshold           ~ "red3",
      as.numeric(chromosome) %% 2 == 0          ~ "grey60",
      TRUE                                      ~ "grey30"
    )
  ) %>%
  ggplot(aes(x = pos, y = WEIGHTED_FST, colour = point_color)) +
  geom_jitter(size = 3, alpha = 0.9) +
  geom_hline(yintercept = ag_pono_threshold, linetype = 'dashed', colour = 'black') +
  facet_grid(cols = vars(chromosome), space = "free_x", scales = "free_x", switch = "x") +
  scale_colour_identity() +
  labs(x = "Scaffold (MB)", y = expression(F[ST] * " - PonoKai vs Agstation")) +
  theme_classic() +
  theme(
    axis.text       = element_text(size = 17),
    axis.title      = element_text(size = 17),
    axis.text.x     = element_text(size = 0),
    axis.ticks.x    = element_blank(),
    strip.text      = element_text(size = 17),
    legend.position = 'none'
  )

ggarrange(SilencePlot, FSTPONOKAIAGSTATION, nrow = 2)

# Define and export FST outlier windows
fst_ag_pono <- fst_ag_pono %>%
  mutate(outlier = ifelse(WEIGHTED_FST > ag_pono_threshold, "outlier", "background"))

fst_ag_pono %>%
  filter(outlier == "outlier") %>%
  group_by(chromosome) %>%
  tally()

write.csv(fst_ag_pono %>% filter(outlier == "outlier"), "Outliers_fst_ag_pono.csv")

# --- ChurchLawn vs CommonGround ---
fst_cg_church <- read.table("COLDVERSUSHOTSPOT_FST_CHURCHLAWN_VERSUS_COMMONG.windowed.weir.fst", header = T)

fst_cg_church$chromosome <- as.double(gsub("chr", "", fst_cg_church$CHROM))
fst_cg_church <- na.omit(fst_cg_church)

fst_cg_church_threshold <- quantile(fst_cg_church$WEIGHTED_FST, 0.95, na.rm = T)

fst_cg_church_PLOT <- fst_cg_church %>%
  filter(N_VARIANTS > 100) %>%
  mutate(
    chromosome  = factor(chromosome, levels = mixedsort(unique(chromosome))),
    pos         = (BIN_START + BIN_END) / 2,
    point_color = case_when(
      WEIGHTED_FST > fst_cg_church_threshold    ~ "red3",
      as.numeric(chromosome) %% 2 == 0         ~ "grey60",
      TRUE                                     ~ "grey30"
    )
  ) %>%
  ggplot(aes(x = pos, y = WEIGHTED_FST, colour = point_color)) +
  geom_jitter(size = 3, alpha = 0.9) +
  geom_hline(yintercept = fst_cg_church_threshold, linetype = 'dashed', colour = 'black') +
  facet_grid(cols = vars(chromosome), space = "free_x", scales = "free_x", switch = "x") +
  scale_colour_identity() +
  labs(x = "Scaffold (MB)", y = expression(F[ST] * " - ChurchLawn vs CommonGround")) +
  theme_classic() +
  theme(
    axis.text       = element_text(size = 17),
    axis.title      = element_text(size = 17),
    axis.text.x     = element_text(size = 0),
    axis.ticks.x    = element_blank(),
    strip.text      = element_text(size = 17),
    legend.position = 'none'
  )

ggarrange(SilencePlot, FSTPONOKAIAGSTATION, fst_cg_church_PLOT, nrow = 3)

fst_cg_church <- fst_cg_church %>%
  mutate(outlier = ifelse(WEIGHTED_FST > fst_cg_church_threshold, "outlier", "background"))

fst_cg_church %>%
  filter(outlier == "outlier") %>%
  group_by(chromosome) %>%
  tally()

# ==============================================================================
# SECTION 16: ALLELE FREQUENCY AROUND GEA LOCI ON CHR2
# ==============================================================================
# Examine derived allele frequency in the vicinity of the top GEA hit
# near bab2 on chromosome 2

chr2_freqs <- read.table("chr2_Allelefreq.frq", header = T, fill = TRUE, comment.char = "")
chr2_freqs  <- chr2_freqs %>% filter(FREQ_Anc < 1)   # Remove fixed sites

ggplot(chr2_freqs, aes(x = POS, y = FREQ_Anc)) +
  geom_point() +
  xlim(41402551, 41573889) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  theme_bw() +
  labs(x = "Position (chr2)", y = "Allele freq (derived allele)") +
  geom_vline(xintercept = 41515212, colour = 'black') +   # bab2 start
  geom_vline(xintercept = 41573889, colour = 'black') +   # bab2 end
  geom_vline(xintercept = 41511790, colour = 'purple')    # Closest significant GEA hit

# ==============================================================================
# SECTION 17: TAJIMA'S D AT GEA CANDIDATE LOCI (PER ISLAND)
# ==============================================================================
# Check for reduced Tajima's D (negative values = recent positive selection)
# near the top GEA hit on chr2, separately for each island population

# Shared theme settings for Tajima's D plots
tajima_theme <- list(
  facet_grid(cols = vars(as.numeric(CHROM)), space = "free_x", scales = "free_x", switch = "x"),
  labs(x = "Scaffold (MB)"),
  theme_classic(),
  theme(
    legend.position = 'none',
    axis.text.y     = element_text(size = 13),
    axis.text.x     = element_text(size = 13),
    axis.ticks.x    = element_blank()
  ),
  geom_hline(yintercept = -1, linetype = 'dashed', colour = 'red', size = 1),
  geom_vline(xintercept = 41511790, colour = 'purple')   # GEA hit position
)

# Kauai
TajimasDKauai %>%
  filter(CHROM == 2) %>% na.omit() %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'grey10') +
  labs(y = "Kaua'i") + tajima_theme

# Oahu
TajimasDOahu %>%
  filter(CHROM == 2) %>% na.omit() %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'grey10') +
  labs(y = "Oahu") + tajima_theme +
  xlim(41402551, 41573889)

# Hilo
TajimasDHilo %>%
  filter(CHROM == 2) %>% na.omit() %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_jitter(size = 2.5, alpha = 0.7, colour = 'grey10') +
  labs(y = "Hilo") + tajima_theme +
  xlim(41402551, 41573889)

# ==============================================================================
# SECTION 18: ABIOTIC (ENVIRONMENTAL) LFMM
# ==============================================================================
# Test whether any SNPs are associated with environmental (climate) variation
# as a comparison to the silence phenotype GEA

# Extract population coordinates from PCA metadata
pca_all_clean <- as.data.frame(pca_all)
pca_all_clean$individuals <- as.character(pca_all_clean$individuals)

population_data <- pca_all_clean %>%
  dplyr::semi_join(MetadataSilence, by = c("individuals" = "Samples")) %>%
  dplyr::slice(match(MetadataSilence$Samples, individuals)) %>%
  dplyr::select(individuals, Location, Latitude, Longitude) %>%
  filter(Location %in% c("KAUAI", "OAHU", "HILO", "MOLOKAI")) %>%
  dplyr::select(Longitude, Latitude) %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Download WorldClim bioclimatic variables at 2.5 arcminute resolution
wclim   <- get_worldclim(coords = population_data, res = 2.5)

# PCA of environmental variables to reduce dimensionality
env_pcs   <- rasterPCA(wclim, spca = TRUE)
pc_values <- terra::extract(env_pcs$map[[1:4]], population_data)
env       <- pc_values %>% dplyr::select(PC1, PC2, PC3, PC4)

# Extract PCA loadings for biplot
loadings    <- env_pcs$model$loadings
loadings_df <- as.data.frame(loadings[, 1:2])
loadings_df$variable <- rownames(loadings_df)

# Biplot: population positions in environmental PC space with loading arrows
BiPlotPCA_EnvironmentalVariables <- ggplot(env, aes(x = PC1, y = PC2)) +
  geom_jitter(alpha = 0.6, width = 0.2) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = Comp.1 * 5, yend = Comp.2 * 5),
               arrow = arrow(length = unit(0.2, "inches")), color = "red") +
  geom_text_repel(data = loadings_df,
                  aes(x = Comp.1 * 5, y = Comp.2 * 5, label = variable),
                  color = "red", max.overlaps = Inf, size = 3) +
  theme_pubr() +
  labs(title = "PCA Biplot of WorldClim Variables",
       x = "PC1 (Temperature-Associated)", y = "PC2 (Precipitation-Associated)") +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey') +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey')

png("~/Desktop/Ochracea_POSTDOC/BiPlotPCA_EnvironmentalVariables.png",
    units = "in", width = 8.0, height = 5.0, res = 300)
BiPlotPCA_EnvironmentalVariables
dev.off()

# Replace BIO codes with full variable names in loadings
bio_variable_names <- c(
  BIO1  = "Annual Mean Temperature",       BIO2  = "Mean Diurnal Range",
  BIO3  = "Isothermality",                 BIO4  = "Temperature Seasonality",
  BIO5  = "Max Temperature of Warmest Month", BIO6 = "Min Temperature of Coldest Month",
  BIO7  = "Temperature Annual Range",      BIO8  = "Mean Temperature of Wettest Quarter",
  BIO9  = "Mean Temperature of Driest Quarter", BIO10 = "Mean Temperature of Warmest Quarter",
  BIO11 = "Mean Temperature of Coldest Quarter", BIO12 = "Annual Precipitation",
  BIO13 = "Precipitation of Wettest Month", BIO14 = "Precipitation of Driest Month",
  BIO15 = "Precipitation Seasonality",     BIO16 = "Precipitation of Wettest Quarter",
  BIO17 = "Precipitation of Driest Quarter", BIO18 = "Precipitation of Warmest Quarter",
  BIO19 = "Precipitation of Coldest Quarter"
)

replace_bio_codes_vector <- function(variable_vector) {
  modified_vector <- variable_vector
  for (bio_code in names(bio_variable_names)) {
    pattern         <- paste0("\\b", tolower(bio_code), "\\b")
    modified_vector <- gsub(pattern, bio_variable_names[bio_code], modified_vector, ignore.case = TRUE)
  }
  return(modified_vector)
}

loadings_df$variable <- replace_bio_codes_vector(loadings_df$variable)

# Run abiotic LFMM using PC1 and PC2 as environmental predictors
PC1_PC2_abiotic <- env %>% select(PC1, PC2)
imputed_vcf     <- algatr::str_impute(gen = vcf_silent, K = 2, entropy = T, repetitions = 1)
ridge_results_abiotic <- lfmm_run(imputed_vcf, PC1_PC2_abiotic, K = 2, lfmm_method = "lasso")

lfmm_table(ridge_results_abiotic$df, order = TRUE)
ridge_outliers_abiotic <- ridge_results_abiotic$df %>% filter(adjusted.pvalue < 0.05)

# ==============================================================================
# SECTION 19: ALLELE FREQUENCIES AT GEA OUTLIER LOCI
# ==============================================================================
# Characterise major allele frequencies at GEA hits on chr2 and chr5,
# and compare with SNP-wise FST between Hawaii and mainland populations

# Load allele frequencies for GEA outlier SNPs
AlleleFreqsSilenceHits <- read.table("AlleleFreqCategoricalLFMM.frq",
                                     header = TRUE, check.names = FALSE, row.names = NULL)
colnames(AlleleFreqsSilenceHits) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "A1_FREQ", "A2_FREQ")

# Compute major and minor allele frequencies
AlleleFreqsSilenceHits$MAJOR_freq <- pmax(AlleleFreqsSilenceHits$A1_FREQ, AlleleFreqsSilenceHits$A2_FREQ)
AlleleFreqsSilenceHits$MINOR_freq <- pmin(AlleleFreqsSilenceHits$A1_FREQ, AlleleFreqsSilenceHits$A2_FREQ)

# Load Hardy-Weinberg test and Hawaii-Mainland FST data, and join to allele freq data
SnpWiseFST_Silence <- read.table("Hawaii_Mainland_CategoricalSNPs.weir.fst", header = T)

AlleleFreqsSilenceHits <- AlleleFreqsSilenceHits %>%
  left_join(HardyWeinbergTest_df_filtered %>% select(POS, P_HWE), by = "POS") %>%
  left_join(SnpWiseFST_Silence %>% select(POS, WEIR_AND_COCKERHAM_FST), by = "POS")

# Set negative FST values to zero (can arise from VCFTOOLS estimation)
AlleleFreqsSilenceHits$WEIR_AND_COCKERHAM_FST[AlleleFreqsSilenceHits$WEIR_AND_COCKERHAM_FST < 0] <- 0

# Remove chr1 (not a focal region) and clean chromosome labels
AlleleFreqsSilenceHits <- AlleleFreqsSilenceHits %>% filter(CHROM != "chr1")
AlleleFreqsSilenceHits$CHROM <- gsub("chr", "", AlleleFreqsSilenceHits$CHROM)

# Chromosome colour scheme
chrom_colors <- c("2" = "olivedrab3", "5" = "plum3")

# Summarise mean and SD of major allele frequency per chromosome
chrom_summary <- AlleleFreqsSilenceHits %>%
  group_by(CHROM) %>%
  summarise(
    mean_MAJOR = mean(MAJOR_freq, na.rm = TRUE),
    sd_MAJOR   = sd(MAJOR_freq, na.rm = TRUE)
  )

# Wilcoxon test: do chr2 and chr5 have different major allele frequencies?
wilcox.test(MAJOR_freq ~ CHROM, data = AlleleFreqsSilenceHits)

# Point-range plot: mean ± 1 SD major allele frequency per chromosome
AlleleFrequenciesofSilentSNPs <- ggplot(chrom_summary, aes(x = CHROM, y = mean_MAJOR, color = CHROM)) +
  geom_pointrange(aes(ymin = mean_MAJOR - sd_MAJOR, ymax = mean_MAJOR + sd_MAJOR), size = 2) +
  scale_color_manual(values = chrom_colors) +
  labs(x = "Chromosome", y = "Major allele frequency", color = "Chromosome") +
  theme_pubr() +
  theme(
    axis.text       = element_text(size = 17),
    axis.title      = element_text(size = 17),
    axis.text.x     = element_text(size = 17),
    legend.position = 'none'
  ) +
  ylim(0.5, 1)

png("AlleleFrequenciesofSilentSNPs.png", units = "in", width = 4.0, height = 4.0, res = 300)
AlleleFrequenciesofSilentSNPs
dev.off()

# ==============================================================================
# SECTION 20: HARDY-WEINBERG EQUILIBRIUM TEST
# ==============================================================================
# Test for departures from HWE at GEA loci, which could indicate balancing
# selection (heterozygote excess) or other non-neutral processes

HardyWeinbergTest_df <- read.table("WholeGenomeHardyWeinbergTest_GP.hwe", header = T)
HardyWeinbergTest_df$CHR <- as.factor(gsub("chr", "", HardyWeinbergTest_df$CHR))

# Plot HWE p-values on chr5
HardyWeinbergTest_df %>%
  filter(CHR == "5") %>% na.omit() %>%
  ggplot(aes(x = POS, y = -log10(P_HWE))) +
  geom_jitter(size = 2.5, alpha = 0.6, colour = 'black') +
  labs(x = "Scaffold (MB)", y = "-log10 (pvalue)") +
  facet_grid(cols = vars(as.numeric(CHR)), space = "free_x", scales = "free_x", switch = "x") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', colour = 'black') +
  theme_classic() +
  theme(legend.position = 'none')

# Subset HWE results to only GEA outlier positions
HardyWeinbergTest_df_filtered <- HardyWeinbergTest_df %>%
  semi_join(ridge_outliers_silence, by = c("POS" = "chr_position"))

write.csv(HardyWeinbergTest_df_filtered, "HardyWeinbergTest_df_filtered_ForPub.csv",
          quote = F, row.names = F)

# ==============================================================================
# SECTION 21: PREPARE PUBLICATION-READY GEA OUTLIER TABLE
# ==============================================================================
# Join LFMM statistical results with gene annotation information

LFMM_gene_information <- read.table("Categorical_Var_LFMM_For_Publication.tsv", header = T)

# Merge outlier SNP stats with nearest gene information
Silence_GEA_Outliers_PublicationReady <- ridge_outliers_silence %>%
  left_join(
    LFMM_gene_information %>% select(Chromosome, End, GeneID, GeneName, DistanceFromGene),
    by = c("chromosome" = "Chromosome", "chr_position" = "End")
  ) %>%
  select(chromosome, chr_position, var, B, score, pvalue, calibrated.pvalue,
         adjusted.pvalue, GeneName, GeneID, DistanceFromGene)

write.csv(Silence_GEA_Outliers_PublicationReady, "Silence_GEA_Outliers_PublicationReady.csv",
          quote = F, row.names = F)

# ==============================================================================
# SECTION 22: GENE ONTOLOGY ENRICHMENT PLOT (CATEGORICAL VARIABLE)
# ==============================================================================

GeneOntologyBioEnrichment <- read.csv("~/Desktop/BIBTEX/SilenceLFMMCATEGORICALENRICHMENT.csv")

Top10_GeneOntologyBioEnrichment <- GeneOntologyBioEnrichment %>%
  filter(P.value < 0.05) %>%
  top_n(10, Fold.Enrichment) %>%
  arrange(desc(Fold.Enrichment))

Top10_GeneOntologyBioEnrichment_Plot <- ggplot(
  Top10_GeneOntologyBioEnrichment,
  aes(x = reorder(Gene.Set.Name, -log10(P.value)), y = -log10(P.value))
) +
  geom_col(fill = 'grey', colour = 'black') +
  coord_flip() +
  theme_pubr() +
  labs(y = "Enrichment score", x = "", title = "") +
  theme(
    axis.text  = element_text(size = 15),
    axis.title = element_text(size = 17),
    strip.text = element_text(size = 17)
  )

png("Top10_GeneOntologyBioEnrichment_Plot.png", units = "in", width = 7.0, height = 5.5, res = 300)
Top10_GeneOntologyBioEnrichment_Plot
dev.off()

# ==============================================================================
# SECTION 23: SELECTION COEFFICIENT ESTIMATION — CLUES2 LRT
# ==============================================================================
# Likelihood ratio test comparing neutral (model A) vs selection (model B)
# for GRIK1 and TLR3 candidate loci

# GRIK1
logLik_A_GRIK1 <- 573.0508
logLik_B_GRIK1 <- 623.4168
D_GRIK1        <- -2 * (logLik_A_GRIK1 - logLik_B_GRIK1)
p_value_GRIK1  <- pchisq(D_GRIK1, df = 1, lower.tail = FALSE)

cat("GRIK1 LRT statistic:", D_GRIK1, "\n")
cat("GRIK1 p-value:", p_value_GRIK1, "\n")

# TLR3
logLik_A_TLR3 <- 31.346
logLik_B_TLR3 <- 32.5445
D_TLR3        <- -2 * (logLik_A_TLR3 - logLik_B_TLR3)

cat("TLR3 LRT statistic:", D_TLR3, "\n")

# ==============================================================================
# SECTION 24: ROBUSTNESS TEST — ALTERNATE CATEGORICAL VARIABLE
# ==============================================================================
# Re-run LFMM with Princeville reassigned to high-song group to test whether
# GEA results are robust to this categorisation decision

MetadataSilenceALTERNATE <- MetadataSilence %>%
  mutate(Categorical_Singing_Variable = case_when(
    Population %in% c("AgStation")                                                              ~ 0,
    Population %in% c("CommonGround", "Hanalei", "KauaiCommunitycollege",
                      "CommunityCenter", "HumaneSociety", "Breadfruit", "Princeville")         ~ 0.25,
    Population %in% c("BYU", "PonoKai", "UH", "ChurchLawn")                                   ~ 0.5
  ))

Silence_Categorical_Variable_ALTERNATE  <- MetadataSilenceALTERNATE %>% select(Categorical_Singing_Variable)
Silence_Continuous_Variable_ALTERNATE   <- MetadataSilenceALTERNATE %>% select(Normalized_Prop_singing)

# Run LFMM with alternate categorical variable
ridge_results_silence_ALTERNATE_TESTING <- lfmm_run(
  vcf_silent, Silence_Categorical_Variable_ALTERNATE, K = 3, lfmm_method = "ridge")

lfmm_table(ridge_results_silence_ALTERNATE_TESTING$df, order = TRUE)

ridge_results_silence_ALTERNATE_TESTING_df <- ridge_results_silence_ALTERNATE_TESTING$df

# Extract significant outliers
Ridge_Silence_Outliers_DF <- ridge_results_silence_ALTERNATE_TESTING_df %>%
  filter(adjusted.pvalue < 0.05)

# Check overlap between alternate and original GEA outlier sets
intersect(Ridge_Silence_Outliers_DF$adjusted.pvalue, ridge_outliers_silence$adjusted.pvalue)
