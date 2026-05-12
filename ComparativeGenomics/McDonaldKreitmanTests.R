##### MK TEST - calculation for Hawaii and Mainland -- for PBS genes -----
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)
custom_colors <- c(
  "PBS outlier_Significant" = "#b30000",         # Vibrant dark red
  "PBS outlier_Not significant" = "#fcaeae",     # Light red
  "background_Significant" = "#666666",          # Dark grey
  "background_Not significant" = "#cccccc"       # Light grey
)


custom_colors2 <- c(
  "PBS outlier_Significant" = "#007f7f",         # Dark turquoise blue
  "PBS outlier_Not significant" = "#add8e6",     # Light blue
  "background_Significant" = "#666666",          # Dark grey
  "background_Not significant" = "#cccccc"       # Light grey
)

# retrieve the dataset -- MK Hawaii
MK.Hawaii<-read_tsv("mk.Hawaii.tsv")

PBS_transcripts <- read.table("PBS_hits_transcripts.txt", header=T)

# identify matches in full dataset with pbs genes in dn/ds data

MK.Hawaii <- MK.Hawaii %>%
  mutate(PBS = case_when(
    transcript  %in% PBS_transcripts$TranscriptsPBS ~ "PBS outlier",
    TRUE ~ "background"
  ))

# remove transcripts where there is little information about divergence or polymorphism -- these are probably unreliable.
MK.Hawaii<-na.omit(MK.Hawaii)

#perform multiple testing correction on raw p-value.
MK.Hawaii$adjusted_pval<-p.adjust(MK.Hawaii$pval, method = "fdr", n = length(MK.Hawaii$pval))

# remove sites where neutrality index is infinite because dS is probably sky high.
MK.Hawaii<-MK.Hawaii %>%
  filter(odds_ni != "Inf")

ggplot(MK.Hawaii, aes(x = dS)) +
  geom_histogram()+
  geom_vline(xintercept = 100)

# now test whether the proportion of genes in PBS under adaptive evolution is higher than expected.
prop_data_stacked <- MK.Hawaii %>%
  mutate(significant = ifelse(adjusted_pval < 0.05, "Significant", "Not significant")) %>%
  group_by(PBS, significant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PBS) %>%
  mutate(
    proportion = n / sum(n),
    fill_group = paste(PBS, significant, sep = "_")
  )

ProportionAdaptiveGenesHawaiiStackedBar<-ggplot(prop_data_stacked, aes(x = PBS, y = proportion, fill = fill_group)) +
  geom_col(color = "black", width=0.95) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  labs(
    title = "",
    x = "",
    y = "Proportion of genes under adaptive evolution"
  ) +
  theme_classic()+
  theme(
    axis.text = element_text(size = 17),         # Increase overall axis text size
    axis.text.x = element_blank()                # Remove x-axis text only
  )+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())


MK.Hawaii.Outlier<-MK.Hawaii %>%
  filter(adjusted_pval < 0.05)

# 32 outliers in Hawaii. 

### But now I want to explore the DoS selection of PBS genes and the Neutrality index.

DoS_Hawaii<-ggplot(MK.Hawaii, aes(x = PBS, y = dos, colour = PBS, shape = PBS)) +
  geom_violin(aes(fill = PBS), width = 0.9, alpha = 0.3, color = NA) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, colour = "black", fatten = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_blank()
  ) +
  scale_fill_manual(values = c("grey", "tomato")) +         # Custom fill colors
  labs(y = "Direction of Selection (DoS)", x = "")+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
  
  
# Direction of selection appears to be negative for most genes and PBS genes appear to be pretty much the same.
# This implies long-term negative selection on these genes -- they're conserved and probably pretty important.

ggplot(MK.Mainland, aes(x=dos, y=-log10(adjusted_pval)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))+
  theme_pubr()
  

# What about the neutrality index.



###--------------------- what about the mainland population. -----------------
#------------------------------------------------------------------------------
  
#------------------------------------------------------------------------------
# retrieve the dataset -- MK Hawaii
MK.Mainland<-read_tsv("mk.Mainland.tsv")

# identify matches in full dataset with pbs genes in dn/ds data

MK.Mainland<- MK.Mainland %>%
  mutate(PBS = case_when(
    transcript  %in% PBS_transcripts$TranscriptsPBS ~ "PBS outlier",
    TRUE ~ "background"
  ))

# remove transcripts where there is little information about divergence or polymorphism -- these are probably unreliable.
MK.Mainland<-na.omit(MK.Mainland)

#perform multiple testing correction on raw p-value.
MK.Mainland$adjusted_pval<-p.adjust(MK.Mainland$pval, method = "fdr", n = length(MK.Mainland$pval))

# remove sites where neutrality index is infinite because dS is probably sky high.
MK.Mainland<-MK.Mainland %>%
  filter(odds_ni != "Inf")

# now test whether the proportion of genes in PBS under adaptive evolution is higher than expected.
prop_data_stacked2 <- MK.Mainland %>%
  mutate(significant = ifelse(adjusted_pval < 0.05, "Significant", "Not significant")) %>%
  group_by(PBS, significant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PBS) %>%
  mutate(
    proportion = n / sum(n),
    fill_group = paste(PBS, significant, sep = "_")
  )

ProportionAdaptiveGenesMainlandStackedBar<-ggplot(prop_data_stacked2, aes(x = PBS, y = proportion, fill = fill_group)) +
  geom_col(color = "black", width=0.95) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors2, guide = "none") +
  labs(
    title = "",
    x = "",
    y = "Proportion of genes under adaptive evolution"
  ) +
  theme_classic()+
  theme(
    axis.text = element_text(size = 17),         # Increase overall axis text size
    axis.text.x = element_blank()                # Remove x-axis text only
  )+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())


MK.Mainland %>%
  filter(adjusted_pval < 0.05) %>%
  filter(PBS == "PBS outlier") %>%
  ggplot(aes(x=PBS, y=dos))+
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, dotsize=1.2)


MK.Mainland.Outlier<-MK.Mainland %>%
  filter(adjusted_pval < 0.05)
# In the mainland, there are fewer adaptive genes than Hawaii but this difference is probably not significant.
# 21 adaptive genes in mainland.
# 32 adaptive genes in Hawaii.


# Ok now you've calculated the proportion


# Combined together what genes are under long-term adaptive evolution.

SignificantMKOutliersinHawaiiMainland<-intersect(MK.Hawaii.Outlier$transcript, MK.Mainland.Outlier$transcript)

write.table(SignificantMKOutliersinHawaiiMainland, "SignificantMKOutliersinHawaiiMainland.txt", row.names = F, quote=F)



## What about the direction of selection for these genes: 

Dos_Selection_Mainland<-ggplot(MK.Mainland, aes(x = PBS, y = dos, colour = PBS, shape = PBS)) +
  geom_violin(aes(fill = PBS), width = 0.9, alpha = 0.3, color = NA) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, colour = "black", fatten = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_blank()
  ) +
  scale_fill_manual(values = c("grey", "#007f7f")) +         # Custom fill colors
  labs(y = "Direction of Selection (DoS)", x = "")+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

t.test(dos ~ PBS, data = MK.Hawaii, var.equal=F)
t.test(dos ~ PBS, data = MK.Mainland, var.equal=F)

DoS_Final<-ggarrange(DoS_Hawaii, Dos_Selection_Mainland, nrow=2, ncol=1)

# Plot figures.
png("DirectionOfSelectionHawaiiMainland.png", units="in", width=4.0, height=6.0, res=300)
DoS_Final
dev.off()

Alpha_StackedBars<-ggarrange(ProportionAdaptiveGenesHawaiiStackedBar, ProportionAdaptiveGenesMainlandStackedBar, nrow=2, ncol=1)

png("AlphaBarChartHawaiiMainland.png", units="in", width=4.0, height=6.0, res=300)
Alpha_StackedBars
dev.off()


# 2x2 contingency table test.
# Hawaii populations.
MK_test_Hawaii <- xtabs(n ~ PBS + significant, data = prop_data_stacked)

# View the table
MK_test_Hawaii

# Run chi-squared test
chisq.test(MK_test_Hawaii)


# Mainland.
# Reshape into 2x2 contingency table
MK_test_mainland <- xtabs(n ~ PBS + significant, data = prop_data_stacked2)

# View the table
MK_test_mainland

# Run chi-squared test
chisq.test(MK_test_mainland)
prop.table(MK_test_mainland, margin = 1)
prop.table(MK_test_Hawaii, margin = 1)



### identify genes that are under adaptive evolution in Ormia. ochracea
# DoS on the x axis and then p-value of McDonald Kreitman.
