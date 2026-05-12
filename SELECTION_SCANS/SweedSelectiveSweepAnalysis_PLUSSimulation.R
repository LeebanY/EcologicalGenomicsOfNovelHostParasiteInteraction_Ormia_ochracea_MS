# ==============================================================================
# LOAD REQUIRED PACKAGES
# ==============================================================================

library(tidyverse)

# Define colours for background vs outlier points
outlier_colours <- c("gray10", "orangered2")


# ==============================================================================
# SWEEED ANALYSIS: HAWAIIAN O. OCHRACEA
# ==============================================================================

# Load chromosome-specific SweeD output files
Chr1_Sweed <- read.table(
  "sweed-results/SweeD_Report.Chr1_Hawaiian_Ormia",
  header = TRUE
)

Chr2_Sweed <- read.table(
  "sweed-results/SweeD_Report.Chr2_Hawaiian_Ormia",
  header = TRUE
)

Chr3_Sweed <- read.table(
  "sweed-results/SweeD_Report.Chr3_Hawaiian_Ormia",
  header = TRUE
)

Chr4_Sweed <- read.table(
  "sweed-results/SweeD_Report.Chr4_Hawaiian_Ormia",
  header = TRUE
)

Chr5_Sweed <- read.table(
  "sweed-results/SweeD_Report.Chr5_Hawaiian_Ormia",
  header = TRUE
)


# ==============================================================================
# ADD CHROMOSOME IDENTIFIERS
# ==============================================================================

Chr1_Sweed$chr <- 1
Chr2_Sweed$chr <- 2
Chr3_Sweed$chr <- 3
Chr4_Sweed$chr <- 4
Chr5_Sweed$chr <- 5


# ==============================================================================
# COMBINE ALL CHROMOSOMES INTO A SINGLE DATAFRAME
# ==============================================================================

combined_sweed_df <- rbind(
  Chr1_Sweed,
  Chr2_Sweed,
  Chr3_Sweed,
  Chr4_Sweed,
  Chr5_Sweed
)


# ==============================================================================
# EXPORT OUTLIER WINDOWS
# ==============================================================================

OutlierSweedTest <- combined_sweed_df %>%
  filter(outlier == "outlier")

write.csv(
  OutlierSweedTest,
  "OutlierSweedTest.csv",
  quote = FALSE,
  row.names = FALSE
)


# ==============================================================================
# CHROMOSOME 2 SWEEP PLOT
# ==============================================================================

Chr2_Sweep_Plot <- combined_sweed_df %>%
  
  filter(chr == "2") %>%
  
  ggplot(
    aes(
      x = Position,
      y = Likelihood,
      colour = outlier
    )
  ) +
  
  geom_jitter(size = 4, alpha = 1) +
  
  scale_colour_manual(values = outlier_colours) +
  
  facet_grid(
    cols = vars(as.numeric(chr)),
    space = "free_x",
    scales = "free_x",
    switch = "x"
  ) +
  
  labs(
    x = "Scaffold (MB)",
    y = "CLR"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  
  # Focus on candidate sweep region
  xlim(50000000, 60000000)


# ==============================================================================
# CHROMOSOME 3 SWEEP PLOT
# ==============================================================================

Chr3_Sweep_Plot <- combined_sweed_df %>%
  
  filter(chr == "3") %>%
  
  ggplot(
    aes(
      x = Position,
      y = Likelihood,
      colour = outlier
    )
  ) +
  
  geom_jitter(size = 4, alpha = 1) +
  
  scale_colour_manual(values = outlier_colours) +
  
  facet_grid(
    cols = vars(as.numeric(chr)),
    space = "free_x",
    scales = "free_x",
    switch = "x"
  ) +
  
  labs(
    x = "Scaffold (MB)",
    y = "CLR"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  
  # Focus on candidate sweep region
  xlim(6000000, 10000000)


# ==============================================================================
# CHROMOSOME 4 SWEEP PLOT
# ==============================================================================

Chr4_Sweep_Plot <- combined_sweed_df %>%
  
  filter(chr == "4") %>%
  
  ggplot(
    aes(
      x = Position,
      y = Likelihood,
      colour = outlier
    )
  ) +
  
  geom_jitter(size = 4, alpha = 1) +
  
  scale_colour_manual(values = outlier_colours) +
  
  facet_grid(
    cols = vars(as.numeric(chr)),
    space = "free_x",
    scales = "free_x",
    switch = "x"
  ) +
  
  labs(
    x = "Scaffold (MB)",
    y = "CLR"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  
  # Focus on candidate sweep region
  xlim(50000000, 58000000)


# ==============================================================================
# IDENTIFY GENOME-WIDE SWEEP OUTLIERS
# ==============================================================================

# Calculate 95th percentile CLR threshold
quantile(
  combined_sweed_df$Likelihood,
  0.95,
  na.rm = TRUE
)

# Define threshold
my_threshold <- quantile(
  combined_sweed_df$Likelihood,
  0.95,
  na.rm = TRUE
)

# Label outlier windows
combined_sweed_df <- combined_sweed_df %>%
  mutate(
    outlier = ifelse(
      Likelihood > my_threshold,
      "outlier",
      "background"
    )
  )

# Count outliers per chromosome
OutlierSweedTest <- combined_sweed_df %>%
  filter(outlier == "outlier") %>%
  group_by(chr) %>%
  tally()


# ==============================================================================
# GENOME-WIDE SWEEP PLOT
# ==============================================================================

# CLR threshold based on bottleneck simulations
threshold_SWEED <- 6.348

SweedPlot <- combined_sweed_df %>%
  
  mutate(
    chr = factor(
      chr,
      levels = mixedsort(unique(chr))
    ),
    
    # Define significant sweep windows
    is_outlier = Likelihood > threshold_SWEED
  ) %>%
  
  ggplot(aes(x = Position, y = Likelihood)) +
  
  # Background windows
  geom_jitter(
    data = . %>% filter(!is_outlier),
    aes(colour = as.factor(as.numeric(factor(chr)) %% 2)),
    size = 2.5,
    alpha = 1
  ) +
  
  # Significant sweep windows
  geom_jitter(
    data = . %>% filter(is_outlier),
    colour = "red3",
    size = 2.5,
    alpha = 1
  ) +
  
  # Significance threshold
  geom_hline(
    yintercept = threshold_SWEED,
    linetype = "dashed"
  ) +
  
  facet_grid(
    cols = vars(chr),
    space = "free_x",
    scales = "free_x",
    switch = "x"
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "grey60",
      "1" = "grey80"
    )
  ) +
  
  labs(
    x = "",
    y = "CLR"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 17),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 17),
    strip.background = element_blank(),
    legend.position = "none"
  )


# ==============================================================================
# MAINLAND SWEEP ANALYSIS: ARIZONA
# ==============================================================================

# Load chromosome-specific Arizona SweeD output
Chr1_Sweed_A <- read.table(
  "sweed-results_Arizona/SweeD_Report.Chr1_Arizona_Ormia",
  header = TRUE
)

Chr2_Sweed_A <- read.table(
  "sweed-results_Arizona/SweeD_Report.Chr2_Arizona_Ormia",
  header = TRUE
)

Chr3_Sweed_A <- read.table(
  "sweed-results_Arizona/SweeD_Report.Chr3_Arizona_Ormia",
  header = TRUE
)

Chr4_Sweed_A <- read.table(
  "sweed-results_Arizona/SweeD_Report.Chr4_Arizona_Ormia",
  header = TRUE
)

Chr5_Sweed_A <- read.table(
  "sweed-results_Arizona/SweeD_Report.Chr5_Arizona_Ormia",
  header = TRUE
)


# ==============================================================================
# ADD CHROMOSOME IDENTIFIERS
# ==============================================================================

Chr1_Sweed_A$chr <- 1
Chr2_Sweed_A$chr <- 2
Chr3_Sweed_A$chr <- 3
Chr4_Sweed_A$chr <- 4
Chr5_Sweed_A$chr <- 5


# ==============================================================================
# COMBINE ARIZONA DATA
# ==============================================================================

combined_sweed_df_Arizona <- rbind(
  Chr1_Sweed_A,
  Chr2_Sweed_A,
  Chr3_Sweed_A,
  Chr4_Sweed_A,
  Chr5_Sweed_A
)


# ==============================================================================
# IDENTIFY ARIZONA OUTLIERS
# ==============================================================================

quantile(
  combined_sweed_df_Arizona$Likelihood,
  0.99,
  na.rm = TRUE
)

my_threshold <- quantile(
  combined_sweed_df_Arizona$Likelihood,
  0.99,
  na.rm = TRUE
)

combined_sweed_df_Arizona <- combined_sweed_df_Arizona %>%
  mutate(
    outlier = ifelse(
      Likelihood > my_threshold,
      "outlier",
      "background"
    )
  )

# Count outliers
combined_sweed_df %>%
  filter(outlier == "outlier") %>%
  group_by(chr) %>%
  tally()

# Extract outlier windows
OutlierSweedTest_Arizona <- combined_sweed_df_Arizona %>%
  filter(outlier == "outlier")


# ==============================================================================
# GENOME-WIDE SWEEP PLOT: ARIZONA
# ==============================================================================

combined_sweed_df_Arizona %>%
  
  mutate(
    chr = factor(
      chr,
      levels = mixedsort(unique(chr))
    ),
    
    is_outlier = Likelihood > 4.683955
  ) %>%
  
  ggplot(aes(x = Position, y = Likelihood)) +
  
  geom_jitter(
    data = . %>% filter(!is_outlier),
    aes(colour = as.factor(as.numeric(factor(chr)) %% 2)),
    size = 2.5,
    alpha = 1
  ) +
  
  geom_jitter(
    data = . %>% filter(is_outlier),
    colour = "red3",
    size = 2.5,
    alpha = 1
  ) +
  
  geom_hline(
    yintercept = 4.683955,
    linetype = "dashed"
  ) +
  
  facet_grid(
    cols = vars(chr),
    space = "free_x",
    scales = "free_x",
    switch = "x"
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "grey60",
      "1" = "grey30"
    )
  ) +
  
  labs(
    x = "",
    y = "CLR (Arizona)",
    title = ""
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 17),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 17),
    legend.position = "none"
  )


# ==============================================================================
# MAINLAND SWEEP ANALYSIS: CALIFORNIA
# ==============================================================================

# Load chromosome-specific California SweeD output
Chr1_Sweed_C <- read.table(
  "sweed-results_California/SweeD_Report.Chr1_California_Ormia",
  header = TRUE
)

Chr2_Sweed_C <- read.table(
  "sweed-results_California/SweeD_Report.Chr2_California_Ormia",
  header = TRUE
)

Chr3_Sweed_C <- read.table(
  "sweed-results_California/SweeD_Report.Chr3_California_Ormia",
  header = TRUE
)

Chr4_Sweed_C <- read.table(
  "sweed-results_California/SweeD_Report.Chr4_California_Ormia",
  header = TRUE
)

Chr5_Sweed_C <- read.table(
  "sweed-results_California/SweeD_Report.Chr5_California_Ormia",
  header = TRUE
)


# ==============================================================================
# ADD CHROMOSOME IDENTIFIERS
# ==============================================================================

Chr1_Sweed_C$chr <- 1
Chr2_Sweed_C$chr <- 2
Chr3_Sweed_C$chr <- 3
Chr4_Sweed_C$chr <- 4
Chr5_Sweed_C$chr <- 5


# ==============================================================================
# COMBINE CALIFORNIA DATA
# ==============================================================================

combined_sweed_df_California <- rbind(
  Chr1_Sweed_C,
  Chr2_Sweed_C,
  Chr3_Sweed_C,
  Chr4_Sweed_C,
  Chr5_Sweed_C
)


# ==============================================================================
# IDENTIFY CALIFORNIA OUTLIERS
# ==============================================================================

quantile(
  combined_sweed_df_California$Likelihood,
  0.99,
  na.rm = TRUE
)

my_threshold <- quantile(
  combined_sweed_df_California$Likelihood,
  0.99,
  na.rm = TRUE
)

combined_sweed_df_California <- combined_sweed_df_California %>%
  mutate(
    outlier = ifelse(
      Likelihood > my_threshold,
      "outlier",
      "background"
    )
  )

# Count outliers
combined_sweed_df %>%
  filter(outlier == "outlier") %>%
  group_by(chr) %>%
  tally()

# Extract outlier windows
OutlierSweedTest_California <- combined_sweed_df_California %>%
  filter(outlier == "outlier")


# ==============================================================================
# GENOME-WIDE SWEEP PLOT: CALIFORNIA
# ==============================================================================

combined_sweed_df_California %>%
  
  mutate(
    chr = factor(
      chr,
      levels = mixedsort(unique(chr))
    ),
    
    is_outlier = Likelihood > 3.660814
  ) %>%
  
  ggplot(aes(x = Position, y = Likelihood)) +
  
  geom_jitter(
    data = . %>% filter(!is_outlier),
    aes(colour = as.factor(as.numeric(factor(chr)) %% 2)),
    size = 2.5,
    alpha = 1
  ) +
  
  geom_jitter(
    data = . %>% filter(is_outlier),
    colour = "red3",
    size = 2.5,
    alpha = 1
  ) +
  
  geom_hline(
    yintercept = 3.660814,
    linetype = "dashed"
  ) +
  
  facet_grid(
    cols = vars(chr),
    space = "free_x",
    scales = "free_x",
    switch = "x"
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "grey60",
      "1" = "grey30"
    )
  ) +
  
  labs(
    x = "",
    y = "CLR (California)",
    title = ""
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 17),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 17),
    legend.position = "none"
  )


# ==============================================================================
# BOTTLENECK SIMULATIONS
# ==============================================================================

# Simulated CLR values under bottleneck demographic model
SimulationResult_Sweed <- read.table(
  "SweeD_Report.SimulationBottleneck_Ormia",
  header = TRUE
)


# ==============================================================================
# MEAN OBSERVED CLR VALUES BY CHROMOSOME
# ==============================================================================

OutlierSweedTest %>%
  group_by(chr) %>%
  summarise(
    mean_Likelihood = mean(Likelihood, na.rm = TRUE)
  )


# ==============================================================================
# COMPARE SIMULATED VS OBSERVED CLR VALUES
# ==============================================================================

SimulationVsObservedValues <- ggplot(
  SimulationResult_Sweed,
  aes(x = Likelihood)
) +
  
  geom_histogram(
    fill = "grey80",
    colour = "black"
  ) +
  
  # Observed mean CLR values
  geom_vline(
    xintercept = 10.4,
    linetype = "dashed",
    colour = "tomato1",
    size = 2
  ) +
  
  geom_vline(
    xintercept = 11.1,
    linetype = "dashed",
    colour = "cyan3",
    size = 2
  ) +
  
  geom_vline(
    xintercept = 8.36,
    linetype = "dashed",
    colour = "darkgreen",
    size = 2
  ) +
  
  theme_pubr() +
  
  labs(
    x = "Composite Likelihood Ratio (Bottleneck Simulation)",
    y = "count (SNPs)"
  ) +
  
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 17),
    strip.text = element_text(size = 17),
    legend.position = "none"
  )


# ==============================================================================
# SAVE SIMULATION HISTOGRAM
# ==============================================================================

png(
  "SweedSimulation_HistogramOfSimulatedValuesWithObservedMeans.png",
  units = "in",
  width = 8,
  height = 3,
  res = 300
)

SimulationVsObservedValues

dev.off()
