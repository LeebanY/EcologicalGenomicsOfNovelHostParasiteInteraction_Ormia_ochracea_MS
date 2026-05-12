# ==============================================================================
# COMBINING GWAS + SELECTION SCANS
# ==============================================================================
# This script:
#   1. Combines PBS, XP-EHH, and SweeD genome scan plots
#   2. Creates a Venn diagram showing overlap between tests
#   3. Focuses on three candidate behavioural regions:
#        - GRIK1 (Chr2)
#        - APP (Chr3)
#        - GRIN2B (Chr4)
#   4. Plots sitewise FST values for candidate genes
#
# Required packages:
#   tidyverse
#   ggplot2
#   ggpubr
#   VennDiagram
#   grid
#   gtools
# ==============================================================================


# ==============================================================================
# SECTION 1: COMBINE GENOME-WIDE SELECTION SCANS
# ==============================================================================

# Save combined PBS, XP-EHH, and SweeD plots
png(
  "GWAS_Selection_SCANS.png",
  units = "in",
  width = 13,
  height = 4.5,
  res = 300
)

ggarrange(
  PBS_threeway_plot,
  XPEHH_plot,
  SweedPlot,
  nrow = 3,
  ncol = 1
)

dev.off()

# Display PBS plot in viewer
PBS_threeway_plot


# ==============================================================================
# SECTION 2: VENN DIAGRAM OF OVERLAPPING OUTLIERS
# ==============================================================================

# Install package if necessary
# install.packages("VennDiagram")

library(VennDiagram)

# Draw overlap between PBS, XP-EHH, and SweeD outliers
draw.triple.venn(
  area1 = 298,   # PBS outliers
  area2 = 355,   # XP-EHH outliers
  area3 = 183,   # SweeD outliers
  
  n12   = 178,   # PBS + XP-EHH overlap
  n23   = 124,   # XP-EHH + SweeD overlap
  n13   = 120,   # PBS + SweeD overlap
  n123  = 67,    # Shared by all three tests
  
  category = c("PBS", "XP-EHH", "SweeD"),
  fill = c("red", "blue", "green")
)

dev.off()


# ==============================================================================
# SECTION 3: PREPARE GENE ANNOTATION DATA
# ==============================================================================

# Convert GFF object to dataframe
gff_df <- as.data.frame(gff)

# Extract only gene entries and retain relevant columns
genes <- gff_df %>%
  filter(type == "gene") %>%
  select(
    chromosome = seqnames,
    start,
    end,
    strand,
    Name
  ) %>%
  mutate(
    chromosome = factor(
      chromosome,
      levels = levels(PBS_threeway_plot$data$chromosome)
    )
  )


# ==============================================================================
# SECTION 4: FUNCTION TO GENERATE CANDIDATE REGION PBS PLOTS
# ==============================================================================

# This function generates a PBS plot for a focal chromosome region
plot_candidate_region <- function(data,
                                  chromosome_id,
                                  genes_df,
                                  x_limits,
                                  gene_name,
                                  vline_start,
                                  vline_end,
                                  threshold) {
  
  data %>%
    filter(chromosome == chromosome_id) %>%
    na.omit() %>%
    
    mutate(
      chr_position = ((start + end) / 2) / 1e6,
      
      chromosome = factor(
        chromosome,
        levels = gtools::mixedsort(unique(chromosome))
      ),
      
      # Highlight significant windows
      point_color = case_when(
        PBSn1A > threshold ~ "red3",
        as.numeric(factor(chromosome)) %% 2 == 0 ~ "grey60",
        TRUE ~ "grey30"
      )
    ) %>%
    
    ggplot(aes(x = chr_position,
               y = PBSn1A,
               colour = point_color)) +
    
    geom_jitter(size = 2.5, alpha = 1) +
    
    facet_grid(
      cols = vars(chromosome),
      space = "free_x",
      scales = "free_x",
      switch = "x"
    ) +
    
    # Significance threshold
    geom_hline(
      yintercept = threshold,
      linetype = "dashed"
    ) +
    
    # Candidate gene boundaries
    geom_vline(xintercept = vline_start, colour = "black") +
    geom_vline(xintercept = vline_end, colour = "black") +
    
    # Gene density rug
    geom_rug(
      data = genes_df,
      aes(x = ((start + end) / 2) / 1e6),
      sides = "t",
      inherit.aes = FALSE,
      colour = "darkred",
      length = unit(0.05, "npc")
    ) +
    
    scale_colour_identity() +
    
    xlim(x_limits[1], x_limits[2]) +
    
    labs(
      x = "Position (Mb)",
      y = expression(PBSn1[Hawaii]),
      title = gene_name
    ) +
    
    theme_classic() +
    
    theme(
      axis.text = element_text(size = 17),
      axis.title = element_text(size = 17, face = "bold"),
      strip.text = element_text(size = 17),
      strip.background = element_blank(),
      legend.position = "none"
    )
}


# ==============================================================================
# SECTION 5: CHROMOSOME 2 — GRIK1
# ==============================================================================

chr2_genes <- genes %>%
  filter(chromosome == "2")

CHR2_PBS <- plot_candidate_region(
  data         = PBS_threeway,
  chromosome_id = "2",
  genes_df     = chr2_genes,
  x_limits     = c(50.557906, 59.85),
  gene_name    = "GRIK1",
  vline_start  = 58.557906,
  vline_end    = 58.565988,
  threshold    = threshold
)


# ==============================================================================
# SECTION 6: CHROMOSOME 3 — APP
# ==============================================================================

chr3_genes <- genes %>%
  filter(chromosome == "3")

CHR3_PBS <- plot_candidate_region(
  data          = PBS_threeway,
  chromosome_id = "3",
  genes_df      = chr3_genes,
  x_limits      = c(5, 10),
  gene_name     = "APP",
  vline_start   = 7.859926,
  vline_end     = 7.869537,
  threshold     = threshold
)


# ==============================================================================
# SECTION 7: CHROMOSOME 4 — GRIN2B
# ==============================================================================

chr4_genes <- genes %>%
  filter(chromosome == "4")

CHR4_PBS <- plot_candidate_region(
  data          = PBS_threeway,
  chromosome_id = "4",
  genes_df      = chr4_genes,
  x_limits      = c(48, 52),
  gene_name     = "GRIN2B",
  vline_start   = 50.600002,
  vline_end     = 50.603907,
  threshold     = threshold
)


# ==============================================================================
# SECTION 8: SAVE CANDIDATE PBS PLOTS
# ==============================================================================

png(
  "pbs_top_hits_behaviour.png",
  units = "in",
  width = 10,
  height = 3,
  res = 300
)

ggarrange(
  CHR2_PBS,
  CHR3_PBS,
  CHR4_PBS,
  nrow = 1,
  ncol = 3
)

dev.off()


# ==============================================================================
# SECTION 9: COMBINED CHROMOSOME 2 MULTI-SCAN PANEL
# ==============================================================================

png(
  "GWAS_Selection_SCANS_CHR2.png",
  units = "in",
  width = 4,
  height = 6.5,
  res = 300
)

ggarrange(
  CHR2_PBS,
  CHR2_xpehh,
  CHR2_SWEED,
  nrow = 3,
  ncol = 1
)

dev.off()


# ==============================================================================
# SECTION 10: SITEWISE FST PLOTS FOR CANDIDATE GENES
# ==============================================================================

# ------------------------------------------------------------------------------
# Function to plot sitewise FST
# ------------------------------------------------------------------------------

plot_fst <- function(file_name, gene_title) {
  
  fst_data <- read.table(file_name, header = TRUE) %>%
    na.omit()
  
  ggplot(
    fst_data,
    aes(x = POS, y = WEIR_AND_COCKERHAM_FST)
  ) +
    
    geom_point(
      colour = "black",
      alpha = 1,
      size = 3
    ) +
    
    labs(
      x = "Position",
      y = "FST (Hawaii - Mainland)",
      title = gene_title
    ) +
    
    theme_pubr()
}


# ------------------------------------------------------------------------------
# GRIK1
# ------------------------------------------------------------------------------

GRIK1_FST <- plot_fst(
  "GRIK1.sitewise_fst.weir.fst",
  "GRIK1"
)

GRIK1_FST


# ------------------------------------------------------------------------------
# APP
# ------------------------------------------------------------------------------

APP_FST <- plot_fst(
  "APP.sitewise_fst.weir.fst",
  "APP"
)

APP_FST


# ------------------------------------------------------------------------------
# GRIN2B
# ------------------------------------------------------------------------------

GRIN2B_FST <- plot_fst(
  "GRIN2B.sitewise_fst.weir.fst",
  "GRIN2B"
)

GRIN2B_FST


# ==============================================================================
# END OF SCRIPT
# ==============================================================================