# =============================================================================
# Ormia ochracea — Runs of Homozygosity (ROH) & Inbreeding (FROH) Analysis
#
# Biological rationale:
#   - Runs of Homozygosity (ROH): contiguous stretches of the genome where both
#     haplotypes are identical by descent (IBD). Their presence and length
#     reflect inbreeding history and demographic bottlenecks.
#       * Long ROH (> 1 Mb)  → recent inbreeding or recent bottleneck
#       * Short ROH (< 100kb) → ancient inbreeding or background relatedness
#   - Runs of Heterozygosity (ROHet): the converse — regions maintaining
#     elevated heterozygosity, potentially reflecting balancing selection
#     (e.g. at immune loci or sensory receptor genes).
#   - FROH: the proportion of the genome covered by ROH, a genome-wide
#     inbreeding coefficient analogous to the pedigree-based F.
#
# Input: LD-pruned autosomal SNPs in PLINK PED/MAP format
# Tool:  detectRUNS (Mancin et al. 2020)
# =============================================================================

# --- Load libraries ----------------------------------------------------------
library(vcfR)        # VCF import (not directly used here but often loaded alongside)
library(tidyverse)   # Data manipulation and plotting
library(ggpubr)      # Publication-ready plot helpers
library(detectRUNS)  # ROH/ROHet detection from PLINK PED files
library(ggforce)     # geom_sina for individual-level jitter in violin plots

# =============================================================================
# SECTION 1: DATA INPUT AND ROH DETECTION
# =============================================================================

# --- Input file paths --------------------------------------------------------
# LD-pruned SNP set in PLINK PED/MAP format (autosomal loci only).
# LD pruning beforehand reduces redundancy and prevents inflated ROH estimates
# caused by SNPs in strong LD appearing homozygous by chance rather than IBD.
genotype_file_roh <- "~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/admixture_plink_2/admixtools_pruned_input.ped"
mapFile_roh       <- "~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/admixture_plink_2/admixtools_pruned_input.map"

# --- Run ROH detection -------------------------------------------------------
# consecutiveRUNS.run() uses a sliding-window approach to identify ROH:
# each window of consecutive SNPs is assessed for homozygosity; adjacent
# homozygous windows are merged into runs.
#
# Key parameters:
#   ROHet = FALSE      → detect Runs of Homozygosity (not heterozygosity)
#   maxOppRun = 0      → no heterozygous SNPs tolerated within a run
#                        (strict IBD definition; relaxing this allows for
#                        genotyping error but risks false positives)
#   maxMissRun = 0     → no missing genotypes permitted within a run
#   minSNP = 15        → minimum number of SNPs to constitute a valid ROH
#                        (prevents very short spurious runs)
#   minLengthBps = 1000 → minimum physical length (1 kb) for a reported ROH
#   maxGap = 10^6      → maximum gap between consecutive SNPs within a run
#                        before the run is broken (1 Mb here)
Run_output <- consecutiveRUNS.run(
  genotype_file_roh,
  mapFile_roh,
  ROHet        = FALSE,
  maxOppRun    = 0,
  maxMissRun   = 0,
  minSNP       = 15,
  minLengthBps = 1000,
  maxGap       = 10^6
)

# --- Load pre-computed ROH results from cluster run --------------------------
# The above detection step is computationally intensive and was run on a cluster.
# Results are loaded here from a saved RDS file for downstream analysis.
Run_output <- readRDS("admixture_plink_2/ROHomo_output.rds")

# --- Summary statistics across all runs --------------------------------------
# summaryRuns() reports per-individual and per-population ROH counts, mean/total
# lengths, and SNP density within runs.
# Class = 2 defines ROH length classes for summary binning (e.g. short/long)
# snpInRuns = FALSE skips per-SNP run frequency calculation (saves time)
summaryRuns(Run_output, mapFile_roh, genotype_file_roh, Class = 2, snpInRuns = FALSE)

# --- Filter ROH output -------------------------------------------------------
# Retain only high-confidence, biologically meaningful runs:
#   lengthBps > 100,000  → focus on ROH > 100 kb (likely reflects bottleneck
#                          within the last ~50 generations; shorter ROH are more
#                          ambiguous regarding recency)
#   nSNP > 100           → require dense SNP coverage within the run
#                          (sparse runs are unreliable given genotyping error)
#   group != "ODEPLETA"  → exclude outgroup species O. depleta from analysis
Run_output <- Run_output %>%
  filter(lengthBps > 100000) %>%
  filter(nSNP > 100) %>%
  filter(group != "ODEPLETA")

# --- Exploratory plots of ROH length distribution ---------------------------
# Violin plot: shows full distribution shape per population
# Log10 y-axis spreads the distribution of lengths across orders of magnitude
ggplot(Run_output, aes(x = group, y = lengthBps)) +
  geom_violin() +
  stat_summary(fun = "mean", geom = "point", color = "blue") +
  coord_trans(y = "log10")

# Boxplot version (median + IQR; less informative about multimodality)
ggplot(Run_output, aes(x = group, y = lengthBps)) +
  geom_boxplot() +
  coord_trans(y = "log10")


# =============================================================================
# SECTION 2: PER-POPULATION ROH SUBSETTING AND SUMMARISATION
# For each population, compute the mean ROH length per genomic window
# (defined by chrom + start + end coordinates) across individuals.
# This gives a population-level ROH landscape analogous to a FST/TajD plot.
# =============================================================================

# --- Subset by population ----------------------------------------------------
Kauai_individuals_ROH   <- Run_output %>% filter(group == "KAUAI")
Oahu_individuals_ROH    <- Run_output %>% filter(group == "OAHU")
HILO_individuals_ROH    <- Run_output %>% filter(group == "HILO")
Molokai_individuals_ROH <- Run_output %>% filter(group == "MOLOKAI")
Ventura_individuals_ROH <- Run_output %>% filter(group == "VENTURA")   # Mainland outgroup
Yavapai_individuals_ROH <- Run_output %>% filter(group == "YAVAPAI")   # Mainland outgroup

# --- Summarise: mean ROH length per genomic window per population ------------
# group_by(chrom, from, to) treats each unique ROH window as a unit;
# mean length across individuals sharing that window reveals consistent
# regions of elevated homozygosity (candidate bottleneck or selection targets)

# Hawaiian populations
Kauai_individuals_ROH <- Kauai_individuals_ROH %>%
  filter(nSNP > 100) %>%
  group_by(chrom, from, to) %>%
  summarise(meanlength = mean(lengthBps))

Oahu_individuals_ROH <- Oahu_individuals_ROH %>%
  filter(nSNP > 100) %>%
  group_by(chrom, from, to) %>%
  summarise(meanlength = mean(lengthBps))

HILO_individuals_ROH <- HILO_individuals_ROH %>%
  filter(nSNP > 100) %>%
  group_by(chrom, from, to) %>%
  summarise(meanlength = mean(lengthBps))

# Mainland outgroup populations (Arizona and California)
Ventura_individuals_ROH <- Ventura_individuals_ROH %>%
  filter(nSNP > 100) %>%
  group_by(chrom, from, to) %>%
  summarise(meanlength = mean(lengthBps))

Yavapai_individuals_ROH <- Yavapai_individuals_ROH %>%
  filter(nSNP > 100) %>%
  group_by(chrom, from, to) %>%
  summarise(meanlength = mean(lengthBps))


# =============================================================================
# SECTION 3: GENOME-WIDE ROH VISUALISATION
# Plot mean ROH length along each chromosome (Manhattan-style) for each
# population, to identify chromosomal regions with consistently long ROH.
# =============================================================================

# Diagnostic plot: frequency of SNPs found within runs across chromosome 3
# Useful for checking whether ROH are evenly distributed or cluster at
# specific loci (e.g. near centromeres or low-recombination regions)
plot_SnpsInRuns(
  runs        = Run_output[Run_output$chrom == 3, ],
  genotypeFile = genotype_file_roh,
  mapFile      = mapFile_roh
)

# --- Reusable plotting function for genome-wide ROH length -------------------
# Takes a per-population summarised ROH dataframe (output of Section 2) and
# produces a LOESS-smoothed genome-wide plot faceted by chromosome.
# The midpoint of each ROH window ((from + to) / 2) is used as the x position.
plot_roh_length <- function(df) {
  df %>%
    na.omit() %>%
    mutate(chr_position = (from + to) / 2) %>%   # Window midpoint in bp
    ggplot(aes(x = chr_position, y = meanlength)) +
    geom_smooth(
      method  = "loess",
      alpha   = 0.6,
      span    = 0.1,   # Narrow span = local smoothing; increase for broader trends
      se      = FALSE,
      colour  = "red"
    ) +
    labs(
      x = "Scaffold (MB)",
      y = "Mean Length - ROH"
    ) +
    # Facet by chromosome; free x-scale accounts for different chromosome sizes
    facet_grid(
      cols   = vars(as.numeric(chrom)),
      space  = "free_x",
      scales = "free_x",
      switch = "x"         # Move chromosome labels to bottom of plot
    ) +
    theme_classic() +
    theme(
      axis.text    = element_text(size = 14),
      axis.title   = element_text(size = 13, face = "bold"),
      legend.position = "none",
      axis.text.y  = element_text(size = 13),
      axis.ticks.x = element_blank(),   # Remove x tick marks (positions not informative)
      axis.text.x  = element_blank()    # Remove x axis labels (chromosome label in strip)
    )
}

# Generate per-island plots
KauaiROHplot <- plot_roh_length(Kauai_individuals_ROH)
OahuROHplot  <- plot_roh_length(Oahu_individuals_ROH)
HiloROHplot  <- plot_roh_length(HILO_individuals_ROH)

# Stacked panel: allows visual comparison of ROH landscape across islands
ggarrange(KauaiROHplot, OahuROHplot, HiloROHplot, nrow = 3)

# Mainland outgroup comparisons (plotted separately; not arranged into panel here)
plot_roh_length(Ventura_individuals_ROH)
plot_roh_length(Yavapai_individuals_ROH)

# --- Histograms of mean ROH length per population ----------------------------
# Gives a sense of the distribution of ROH lengths (short vs long tail)
# Note: Hilo and Ventura use the correct object names but Hilo_individuals_ROH
# (capital H) should be HILO_individuals_ROH to match the object defined above
hist(Kauai_individuals_ROH$meanlength)
hist(HILO_individuals_ROH$meanlength)   # Original used Hilo_ (lowercase) — check object name
hist(Yavapai_individuals_ROH$meanlength)
hist(Ventura_individuals_ROH$meanlength)


# =============================================================================
# SECTION 4: FROH — GENOME-WIDE INBREEDING COEFFICIENT
#
# FROH = total length of ROH > threshold / total autosomal genome length
# Analogous to the pedigree inbreeding coefficient F but estimated directly
# from genotype data without requiring known pedigrees.
# Higher FROH → more genomic IBD → more inbreeding / stronger bottleneck.
# Island populations that have undergone founder events are expected to have
# significantly higher FROH than mainland source populations.
# =============================================================================

# Compute FROH per individual, genome-wide
# genome_wide = TRUE uses total autosomal length from the map file as denominator
Froh_output <- Froh_inbreeding(Run_output, mapFile_roh, genome_wide = TRUE)
Froh_output  # Preview: one row per individual with their FROH value

# --- Define population order for plotting ------------------------------------
# Ordered geographically: Hawaiian islands first (from NW to SE in the chain),
# then California mainland, then Arizona mainland
level_order <- c(
  "KAUAI", "OAHU", "HILO", "MOLOKAI",          # Hawaiian islands
  "VENTURA", "LOSANGELES", "SANTABARBARA",       # California
  "SANTACRUZ", "COCHISE", "YAVAPAI"             # Arizona
)

# --- FROH violin + sina plot -------------------------------------------------
# Violin: shows full distribution shape per population
# geom_sina: jitters individual points within the violin width
#            (better than simple jitter as it respects the density shape)
# stat_summary crossbar: marks the population mean FROH
ROH_plot <- ggplot(
  Froh_output,
  aes(x = factor(group, level = level_order), y = Froh_genome)
) +
  geom_violin() +
  geom_sina(alpha = 0.25, colour = "black") +    # Individual data points
  stat_summary(
    fun    = "mean",
    geom   = "crossbar",
    width  = 0.5,
    colour = "red3"                              # Mean marked in red
  ) +
  theme_pubr() +
  labs(
    x = "",
    y = "FROH (Homozygous runs > 100KB)"
  ) +
  coord_flip()   # Horizontal layout; easier to read population labels

# Save figure
png("~/Desktop/Ochracea_POSTDOC/FROH_plot.png",
    units = "in", width = 5.0, height = 5.0, res = 300)
ROH_plot
dev.off()
