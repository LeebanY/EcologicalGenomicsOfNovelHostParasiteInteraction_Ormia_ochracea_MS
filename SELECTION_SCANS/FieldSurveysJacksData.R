# Library
library(tidyverse)
library(ggpubr)


# ==============================================================================
# FIELD SURVEY DATA -- FOR FIGURE 4. 
# Mine and Jack's combined field notes, which illustrate how cricket populations have changed (in terms of adaptive v.s. singing-capable morphs)
# ==============================================================================

SurveyInformation<-read.csv("JacksSurveyData.csv", header=T, fill = TRUE, sep = ",")
# replace na's with zeros which works for this particular use case.
SurveyInformation[is.na(SurveyInformation)] <- 0


# now I want to create an average of singing males across islands -- this is a useful measure because there seems to be little structure within islands.
SurveySummary <- SurveyInformation %>%
  filter(Year != "2022.5") %>%
  filter(Site!= "VacantLot" & Site != "AstronomyCenter") %>%
  mutate(NonSingingMales = TotalMale - Nw_male) %>%
  pivot_longer(cols = c(Nw_male, NonSingingMales),
               names_to = "MaleType",
               values_to = "Count") %>%
  mutate(Proportion = Count / TotalMale) 
SurveySummary$Year<-as.character(SurveySummary$Year)
SurveySummary$Year<-gsub("2022", "2022 - January", SurveySummary$Year)

# Ok but I also want to add my data for 2025 June to this dataframe.

SelectionPressDF_JoinJacksData<-SelectionPressDF %>%
  select(Site, Island, Nw_male, total_males) %>%
  mutate(NonSingingMales = total_males - Nw_male) %>%
  pivot_longer(cols = c(Nw_male, NonSingingMales),
               names_to = "MaleType",
               values_to = "Count") %>%
  mutate(Proportion = Count / total_males) %>%
  mutate(Year = "2022 - June")

# bind my data and Jacks long-term data.
SurveySummary <- bind_rows(SurveySummary, SelectionPressDF_JoinJacksData)
# Make labels more readable.
SurveySummary$MaleType<-gsub("Nw_male", "Singing-capable", SurveySummary$MaleType)
SurveySummary$MaleType<-gsub("NonSingingMales", "Silent morphs", SurveySummary$MaleType)

# reorder sites by Island
SurveySummary <- SurveySummary %>%
  mutate(Site = factor(Site, levels = unique(Site[order(Island, Site)])))

SurveySummary %>%
  group_by(Year, MaleType) %>%
  summarise(ProportionSinging=mean(Proportion)*100)

# Rename the sites so they are more easily seen on the graph.
SurveySummary$Site<-gsub("ChurchLawn", "Hawaii - CL", SurveySummary$Site)
SurveySummary$Site<-gsub("UH", "Hawaii - UH", SurveySummary$Site)
SurveySummary$Site<-gsub("BYU", "Oahu - BYU", SurveySummary$Site)
SurveySummary$Site<-gsub("CommunityCenter", "Oahu - CC", SurveySummary$Site)
SurveySummary$Site<-gsub("AgStation", "Kauai - AS", SurveySummary$Site)
SurveySummary$Site<-gsub("CommonGround", "Kauai - CG", SurveySummary$Site)
SurveySummary$Site<-gsub("PonoKai", "Kauai - PK", SurveySummary$Site)
SurveySummary$Site<-gsub("Kamilo", "Oahu - KP", SurveySummary$Site)
SurveySummary$Site<-gsub("Breadfruit", "Kauai - BI", SurveySummary$Site)
SurveySummary$Site<-gsub("Hanalei", "Kauai - WC", SurveySummary$Site)
SurveySummary$Site<-gsub("Princeville", "Kauai - PV", SurveySummary$Site)
SurveySummary$Site<-gsub("HumaneSociety", "Kauai - HS", SurveySummary$Site)
SurveySummary$Site<-gsub("KauaiCommunitycollege", "Kauai - KCC", SurveySummary$Site)

June2022<-SurveySummary %>%
  filter(Year == "2022 - June") %>%
  filter(MaleType == "Silent morphs")

mean(June2022$Proportion)

EcologicalInfoPlot<-ggplot(SurveySummary, aes(x = Site, y = Proportion, fill = MaleType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Singing-capable" = "lightblue",
                               "Silent morphs" = "maroon")) +
  labs(y = "Percentage of total males", x = "", fill="") +
  theme_pubr()+
  facet_wrap(~ Year, nrow = 1, scales = "free_x", drop = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # smaller x-axis text
    axis.text.y = element_text(size = 12)                         # larger y-axis text
  )

EcologicalInfoPlot <- ggplot(SurveySummary, aes(x = Site, y = Proportion, fill = MaleType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Singing-capable" = "turquoise3",
                               "Silent morphs" = "violetred3")) +
  labs(y = "Percentage of total males", x = "", fill="") +
  theme_pubr() +
  facet_wrap(~ Year, nrow = 1, scales = "free_x", drop = TRUE) +
  theme(
    # remove x-axis text and ticks
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    
    # increase y-axis text size
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    
    # facet strip text larger + remove background
    strip.text = element_text(size = 16),
    strip.background = element_blank(),
    
    # optional: remove facet box line if present
    panel.spacing = unit(1, "lines")
  )+theme(legend.position="bottom")


png("FieldSurveyAcrossYearsPlot.png", units="in", width=14.0, height=4.0, res=300)
EcologicalInfoPlot
dev.off()

# png("FieldSurveyAcrossYearsPlot.png", units="in", width=4.0, height=3.0, res=300)
# FieldSurveyPlotJack
# dev.off()
