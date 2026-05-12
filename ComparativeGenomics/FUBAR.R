##### FUBAR -- SITE LEVEL POSITIVE SELECTION ACROSS TACHINIDS -- for PBS genes -----
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)



###### ANALYSE FUBAR RESULTS FOR PBS GENES - PRODUCE ANOTHER STACKED BAR CHART.
FUBAR_results<-read.table("FUBAR_results.txt", header=T)


FUBAR_results <- FUBAR_results %>%
  mutate(EvidenceforSelection = if_else(
    NoOfSitesUnderDiversifyingPositiveSelection > 0,
    "Yes",
    "No"
  ))



# Step 1: Summarise counts and proportions
selection_summary <- FUBAR_results %>%
  count(EvidenceforSelection) %>%
  mutate(Proportion = n / sum(n))

# Step 2: Plot as stacked bar
FUBAR_test_StackedBar <- ggplot(selection_summary, 
                                aes(x = "Selection", y = Proportion, fill = EvidenceforSelection)) +
  geom_col(width = 0.6, colour = "black") +
  scale_fill_manual(values = c("No" = "#fcaeae", "Yes" = "#b30000")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = NULL, fill = "Evidence for Selection") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 19),
    legend.position = "none"
  )
png("FUBAR_test_StackedBar.png", units="in", width=2.0, height=5.0, res=300)
FUBAR_test_StackedBar+theme(plot.margin = margin(1, 1, 1, 1, "cm"))
dev.off()



# histogram of the number of sites under pervasive selection.
FUBAR_histogram<-FUBAR_results %>%
  filter(NoOfSitesUnderDiversifyingPositiveSelection > 0) %>%
  ggplot(aes(x = NoOfSitesUnderDiversifyingPositiveSelection)) +
  geom_histogram(
    binwidth = 1,                 # Adjust binwidth if needed
    fill = "#b30000",             # Dark green
    colour = "black",             # Outline
    size = 0.4
  ) +
  theme_classic() +
  labs(
    x = "Sites under positive episodic selection",
    y = "count (locally-adaptive genes)"
  ) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )


mean(FUBAR_results$NoOfSitesUnderDiversifyingPositiveSelection)

png("FUBAR_test_Histogram.png", units="in", width=5.0, height=5.0, res=300)
FUBAR_histogram+theme(plot.margin = margin(1, 1, 1, 1, "cm"))
dev.off()