# filtering decisions. 
setwd("~/Desktop/Ochracea_POSTDOC/Ochracea_POSTDOC/filtering_decisions1/")

library(tidyverse)
library(viridis)
library(ggpubr)


var_qual <- read_delim("./Ormia_subset.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

var_qual %>%
  filter(qual > 40) # most quality scores are way above 30, let's go with 40 just to be sure we are keeping the best variants.


# missingness
var_miss <- read.table("Ormia_subset.lmiss", header = T)

a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

summary(var_miss$fmiss)

var_miss %>%
  filter(fmiss > 0.01) # tolerate max missingness of 5-10% of individuals. 


# minor allele frequency filter
var_freq <- read_delim("./Ormia_subset.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)


var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))


a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

summary(var_freq$maf)

# MAF looks correct, but MAF filtering is probably not a good idea given that we've sequenced at high depth and retain only high qual sites.
# Remaining SNPs are likely to be true rare variants and I don't want to bias demographic inference.


# examine the heterozygosity. 

ind_het <- read_delim("./Ormia_subset.het", delim = "\t",
                     col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)


ind_het$Fam_ID<-pca$Location[match(ind_het$ind, pca$individuals)]
ind_het<-na.omit(ind_het)

ind_het <- ind_het %>%
  mutate(state = case_when(
    Fam_ID %in% c("KAUAI", "HILO", "MOLOKAI", "OAHU") ~ "Hawaii",
    Fam_ID %in% c("VENTURA", "SANTABARBARA", "LOSANGELES") ~ "California",
    Fam_ID %in% c("COCHISE", "YAVAPAI", "SANTACRUZ") ~ "Arizona",
    TRUE ~ "OTHER"  # Optional: Handles unexpected values
  ))

Heterozygosity_plot<-ggplot(ind_het, aes(x=state, y=f))+
  geom_boxplot(width=0.5, alpha=0.4)+
  geom_jitter(width=0.1, alpha=0.3, colour='darkgrey')+
  theme_pubr()+
  geom_hline(yintercept=0, linetype='dashed')+
  theme(legend.position="none")+
  labs(x='',
       y='Inbreeding coefficient (F)')+
  theme(text = element_text(size = 15))


png("~/Desktop/Ochracea_POSTDOC/Inbreeding_coefficient_plot.png", units="in", width=5.0, height=5.0, res=300)
Heterozygosity_plot
dev.off()


ind_missingess<-read.table("Ormia_subset.imiss", header = T)
