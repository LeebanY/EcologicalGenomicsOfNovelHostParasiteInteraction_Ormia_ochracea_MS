library(tidyverse)
library(ggpubr)


# ==============================================================================
# GONE2 -- RECENT POPULATION SIZE FLUCTUATIONS.
# Using the software GONE2 to estimate how Ormia populations have changed over the last 150 generations. 
# ==============================================================================

RecentNe<-read.table("Recent_Ne_Rec1.5.2.txt", header = T)

my_palette <- c(
  "Kauai" = "tomato1",    # Orange
  "Oahu" = "deepskyblue",     # Blue
  "Hilo" = "forestgreen",   # Yellow
  "Arizona" = "orange",   # Dark Blue
  "California" = "lightpink1" # Red
)

# pivot the data to a wider format so you can plot. 
RecentNe_long<-RecentNe %>%
  pivot_longer(
    cols = ends_with("diploids"),
    names_to = "Island",
    values_to = "Ne_diploids") %>%
  filter(Generation > 9) 

# formatting. 
RecentNe_long$Island<-gsub("_Ne_diploids", "", RecentNe_long$Island)


# Mutate to add the region column. 
RecentNe_long <- RecentNe_long %>%
  mutate(Region = case_when(
    Island %in% c("California", "Arizona") ~ "Mainland",
    Island %in% c("Kauai", "Oahu", "Hilo") ~ "Hawaii",
    TRUE ~ NA_character_   # catch anything else
  ))

# Plot the data 
NethroughtimePlot<-ggplot(RecentNe_long, aes(x = Generation, y = Ne_diploids, colour = Island)) +
  geom_step(linewidth = 1) +
  theme_classic() +
  labs(y = "Effective population size (Ne)") +
  theme(legend.position = "none") +
  scale_color_manual(values = my_palette) +
  facet_wrap(~ Region, scales = "free_y")+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  )+scale_y_log10(breaks = 10^(1:6),
                    labels = scales::comma)



# Save :-)
png("NeThroughTimeRecent.png", units="in", width=7.0, height=3.0, res=300)
NethroughtimePlot
dev.off()


