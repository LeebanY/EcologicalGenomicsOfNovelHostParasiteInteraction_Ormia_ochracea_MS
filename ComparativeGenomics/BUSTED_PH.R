library(tidyverse)
library(viridis)
library(ggpubr)
library(data.table)
library(ggExtra)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(ggbeeswarm)

################################ BUSTED PH Analysis ###########################
# Which PBS genes are under selection, which have different distributions compared to background genes? #

busted_ph_results<-read.csv("busted_ph_PBS_genes.csv")

# now sort through for genes which are positively selected in Ormia.
busted_ph_results %>%
  filter(LRT_DiversifyingPositiveSelection_Ormia < 0.05)

p.adjust(busted_ph_results$LRT_DiversifyingPositiveSelection_Ormia, method = BH, n = length(busted_ph_results$LRT_DiversifyingPositiveSelection_Ormia))

# so of all the pbs genes (348 in the analyses, 143 genes show evidence of selection in Ormia)

# Now what proportion also show evidence of differences in distribution between Ormia and other Tachinids.
busted_ph_results %>%
  filter(LRT_DiversifyingPositiveSelection_Ormia < 0.05) %>%
  filter(LRT_DifferenceDistributions_Ormia_Tachinids < 0.05)
  
# 81 genes show evidence of selection in Ormia and difference in selective regime between Ormia and other Tachinids.

# Finally identify the best candidates: genes which show (a) evidence of positive selection in Ormia, but (b) not in Tachinids and (c) a statistically significant difference between Ormia and Tachinids.
BestCandidatesForAnalysisOverMacroTime<-busted_ph_results %>%
  filter(LRT_DiversifyingPositiveSelection_Ormia < 0.05) %>%
  filter(LRT_DifferenceDistributions_Ormia_Tachinids < 0.05) %>%
  filter(LRT_DiversifyingPositiveSelection_Tachinids > 0.05)
  
# 15 genes show this pattern. I've checked and they don' t seem to be enriched in any particular function.


# Now look for the number of genes that conversely are under positive selection in Tachinids but not Ormia.
busted_ph_results %>%
  filter(LRT_DifferenceDistributions_Ormia_Tachinids < 0.05)
# There are 18 genes which are under positively selection in Tachinids, but not Ormia.


# Now we ask how many are under positive selection in both and with a difference in dn/ds.
GoodCandidatesConsistentSelectionAcrossTachinidsAndDifference<-busted_ph_results %>%
  filter(LRT_DiversifyingPositiveSelection_Ormia < 0.05) %>%
  filter(LRT_DifferenceDistributions_Ormia_Tachinids < 0.05) %>%
  filter(LRT_DiversifyingPositiveSelection_Tachinids < 0.05)

# 66 genes show positive selection across the tree and show a difference in distribution of positive selection.
# This strongly suggests that PBS genes aren't uniquely under selection in Ormia, instead, they appear to have diversified across Ormia and in Ormia. 
# Ultimately, this supports the view of recurrent selection on PBS genes. 

# Now to find a way of plotting this succincintly.

plot_df <- busted_ph_results %>%
  mutate(
    x = -log10(LRT_DiversifyingPositiveSelection_Ormia),
    y = -log10(LRT_DiversifyingPositiveSelection_Tachinids),
    # cap Inf (and any >30) at 30 so they stay on the panel
    x = if_else(is.infinite(x) | x > 30, 30, x),
    y = if_else(is.infinite(y) | y > 30, 30, y),
    # color by the third variable (red if < 0.05, else grey)
    dist_sig = if_else(LRT_DifferenceDistributions_Ormia_Tachinids < 0.05,
                       "p < 0.05", "p ≥ 0.05")
  )

# Ok this plot isn't that informative but it plots the difference in log-likelihood between Ormia positive selection and Tachinid selection.
ggplot(plot_df, aes(x = x, y = y, colour = dist_sig)) +
  geom_jitter(size = 3, alpha = 0.6, width = 0.5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  scale_colour_manual(values = c("p < 0.05" = "tomato2", "p ≥ 0.05" = "grey50"),
                      name = "Distributions test") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 30))



# Perhaps the easiest thing to do is just plot a barplot of the number of genes which 
# (A) are not under selection in Ormia or Tachinids. n = 61
# (B) selection in Ormia only. n = 29.
# (C) seleciton in Tachinids only but not ormia. n= 144
# (D) selection in both. n= 114
# (E) of those under selection in both, how many show a difference in selective regime. n=66.

# This is fairly strong evidence that genes under local selection in Ormia,
# come from a subset of genes which show evidence of positive selection across Tachinidae, more than half of which seem to have also shown a change in effect size of positive selection.

# Put another way, of all genes under selection in Ormia (n= 143), most are genes also under selection in Tachinidae.

OrmiaSpecificGenesBUSTEDPH<-busted_ph_results %>%
  filter(LRT_DiversifyingPositiveSelection_Ormia < 0.05)

busted_ph_results<-busted_ph_results %>%
  mutate(
    SelectiveRegime = case_when(
      LRT_DiversifyingPositiveSelection_Ormia > 0.05 & LRT_DiversifyingPositiveSelection_Tachinids > 0.05 ~ "No Positive Selection",
      LRT_DiversifyingPositiveSelection_Ormia < 0.05 & LRT_DiversifyingPositiveSelection_Tachinids > 0.05 ~ "Positive Selection in Ormia",
      LRT_DiversifyingPositiveSelection_Ormia > 0.05 & LRT_DiversifyingPositiveSelection_Tachinids < 0.05 ~ "Positive Selection in Background branches",
      LRT_DiversifyingPositiveSelection_Ormia < 0.05 & LRT_DiversifyingPositiveSelection_Tachinids < 0.05 & LRT_DifferenceDistributions_Ormia_Tachinids > 0.05 ~ "Selection Across Tree - No difference between Ormia and background",
      LRT_DiversifyingPositiveSelection_Ormia < 0.05 & LRT_DiversifyingPositiveSelection_Tachinids < 0.05 & LRT_DifferenceDistributions_Ormia_Tachinids < 0.05 ~ "Selection Across Tree - Difference between Ormia and background"
    )
  )
  

# Bar plot of SelectiveRegime counts
BarPlotBustedPH<-busted_ph_results %>%
  mutate(
    SelectiveRegime = factor(
      SelectiveRegime,
      levels = c(
        "No Positive Selection",
        "Positive Selection in Background branches",
        "Positive Selection in Ormia",
        "Selection Across Tree - No difference between Ormia and background",
        "Selection Across Tree - Difference between Ormia and background"
      )
    )
  ) %>%
  count(SelectiveRegime) %>%
  ggplot(aes(x = SelectiveRegime, y = n, fill = SelectiveRegime)) +
  geom_col(width = 0.7, colour = "black") +
  geom_text(aes(label = n), hjust = -0.2, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL,
    y = "Number of Genes",
    title = ""
  ) +
  coord_flip() +
  theme_pubr(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )

library(dplyr)
library(ggplot2)
library(ggpubr)

BarPlotBustedPH <- busted_ph_results %>%
  mutate(
    SelectiveRegime = factor(
      SelectiveRegime,
      levels = c(
        "No Positive Selection",
        "Positive Selection in Background branches",
        "Positive Selection in Ormia",
        "Selection Across Tree - No difference between Ormia and background",
        "Selection Across Tree - Difference between Ormia and background"
      )
    ),
    # Add manual line breaks to long labels
    SelectiveRegime = recode(
      SelectiveRegime,
      "Selection Across Tree - No difference between Ormia and background" =
        "Selection Across Tree:\nNo difference between Ormia and background",
      "Selection Across Tree - Difference between Ormia and background" =
        "Selection Across Tree:\nDifference between Ormia and background"
    )
  ) %>%
  count(SelectiveRegime) %>%
  ggplot(aes(x = SelectiveRegime, y = n, fill = SelectiveRegime)) +
  geom_col(width = 0.7, colour = "black") +
  geom_text(aes(label = n), hjust = -0.2, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Number of Genes", title = "") +
  coord_flip() +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )+
  scale_x_discrete(expand = expansion(add = c(0.3, 0.3)))



png("BarPlotBustedPH.png", units="in", width=12, height=5, res=300)
BarPlotBustedPH
dev.off()

## lollipop plot
# Compute percentages and create lollipop plot
LollipopBustedPH <- busted_ph_results %>%
  mutate(
    Regime = case_when(
      SelectiveRegime == "No Positive Selection" ~ "No selection",
      SelectiveRegime == "Positive Selection in Ormia" ~ "Ormia only",
      SelectiveRegime == "Positive Selection in Background branches" ~ "Background only",
      SelectiveRegime %in% c(
        "Selection Across Tree - No difference between Ormia and background",
        "Selection Across Tree - Difference between Ormia and background"
      ) ~ "Across tree (inc. Ormia)",
      TRUE ~ NA_character_
    ),
    Regime = factor(
      Regime,
      levels = c(
        "No selection",
        "Ormia only",
        "Background only",
        "Across tree (inc. Ormia)"
      )
    ),
    Highlight = ifelse(
      Regime %in% c("No selection", "Ormia only"),
      "No historical selection",
      "Historical selection"
    )
  ) %>%
  count(Regime, Highlight) %>%
  # Calculate percentages
  mutate(percent = n / sum(n) * 100) %>%
  ggplot(aes(x = Regime, y = percent, colour = Highlight)) +
  
  # thick sticks
  geom_segment(aes(x = Regime, xend = Regime, y = 0, yend = percent),
               linewidth = 1.4) +
  
  # large lollipops
  geom_point(size = 13) +
  
  # percentage labels inside points
  geom_text(aes(label = paste0(round(percent, 1), "")),
            colour = "white",
            size = 5,
            fontface = "bold") +
  
  scale_colour_manual(
    values = c(
      "No historical selection" = "grey40",
      "Historical selection"    = "tomato3"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Percentage (%) of genes under episodic selection."
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 12),  # tilted for readability
    axis.title = element_text(size = 14)
  )+scale_y_continuous(expand = expansion(mult = c(0, 0.05)))



LollipopBustedPH

png("LollipopPlotBustedPH.png", units="in", width=6.7, height=4, res=300)
LollipopBustedPH+theme(
  plot.margin = margin(1, 1, 1, 1, "cm") # top, right, bottom, left
)
dev.off()



####b trying boxplots.
OrmiaPSG<-ggplot(
  busted_ph_results,
  aes(
    x = "",
    y = -log10(LRT_DiversifyingPositiveSelection_Ormia)
  )
) +
  geom_violin(fill = "skyblue", alpha = 0.5, colour='grey3') +
  geom_quasirandom(size=3, alpha=0.5, colour="royalblue1") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  coord_cartesian(ylim = c(0, 15.3), clip = "off") +
  theme_classic() +
  labs(x = NULL, y = expression(-log[10](P)))+
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14),  # tilted for readability
    axis.title = element_text(size = 14))

# ------------------------------------------------------------- # 
BackgroundPSG<-ggplot(
  busted_ph_results,
  aes(
    x = "",
    y = -log10(LRT_DiversifyingPositiveSelection_Tachinids)
  )
) +
  geom_violin(fill = "plum2", alpha = 0.5, colour='grey3') +
  geom_quasirandom(size=3, alpha=0.5, colour="purple") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  coord_cartesian(ylim = c(0, 15.3), clip = "off") +
  theme_classic() +
  labs(x = NULL, y = expression(-log[10](P)))+
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14),  # tilted for readability
    axis.title = element_text(size = 14))



# ------------------------------------------------------------- #
PSG_DifferenceinSP <- ggplot(
  busted_ph_results,
  aes(
    x = "",
    y = -log10(LRT_DifferenceDistributions_Ormia_Tachinids)
  )
) +
  geom_violin(fill = "orange", alpha = 0.5, colour = "grey3") +
  geom_quasirandom(size = 3, alpha = 0.5, colour = "orangered") +
  
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.4,
    fatten = 0,
    colour = "black"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  coord_cartesian(ylim = c(0, 15.3), clip = "off") +
  theme_classic() +
  labs(x = NULL, y = expression(-log[10](P))) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14)
  )


ViolinPlotsforPSGs<-ggarrange(BackgroundPSG, OrmiaPSG, PSG_DifferenceinSP, nrow=1,
          ncol=3)


png("ViolinPlotPSGs.png", units="in", width=9, height=4, res=300)
ViolinPlotsforPSGs+theme(plot.margin = margin(1, 1, 1, 1, "cm"))
dev.off()
                         