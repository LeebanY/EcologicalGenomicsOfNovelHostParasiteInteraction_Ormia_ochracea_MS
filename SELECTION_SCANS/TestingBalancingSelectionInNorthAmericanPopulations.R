# This script is to test the selective regime in genes which appear important in host signal detection.
# Descrptive stuff: Intermediate haplotypes are high near the centromeres and sequence diversity (pi) decreases in centromeres probably because of filtering.
# Hypotheses:
# 1. Null: These genes are evolving neutrally and are not under any particular kind of selection. Tajima's D will not distinguish significantly from the background. 
# 2. These genes could be evolving under negative selection (if they are important they could be incredibly consdered) I could look at the non-synonymous polymorphism in these genes and compare to the background.
# 3. These genes could be under directional selection, in which case Tajima's D will be low, the ratio of non-synonymous/synonymous polymorphism will be high, and pi will be low too.
# 4. Under balancing selection, Tajima's D and pi will be high relative to the background, non-synonymous polymorphism will also be high. 
################################################################################


### CALIFORNIA
# Import California data.
TajimasDCalifornia<-read_table("California.TajimasD")
TajimasDCalifornia$CHROM<-gsub("chr", "", TajimasDCalifornia$CHROM)
# Remove windows with poor SNP representation.
TajimasDCalifornia<-TajimasDCalifornia %>%
  filter(N_SNPS > 100)

# Make sure Tajima's D is actually a double.
TajimasDCalifornia$TajimaD<-as.double(TajimasDCalifornia$TajimaD)
# Mean Tajima's D in California: 0.3081467.
mean(TajimasDCalifornia$TajimaD)

# Plot histogram for Cali population. Tajima's D almost centered around 0.
HistogramTajimasD_California<-TajimasDCalifornia %>%
  ggplot(aes(x=TajimaD))+
  geom_histogram(fill='grey', colour='black')+
  theme_pubr()+
  labs(x=expression(paste("Tajima's D - California", italic(" D"))),
       y="count (10kb windows)")+
  geom_vline(xintercept = 0, linetype='dashed', colour='black')

# Genome-wide plot of Tajima's D in Californian populations.
CalifornianTajimasDLandscape<-TajimasDCalifornia %>%
  na.omit() %>%
  ggplot(aes(x=BIN_START, y=TajimaD))+
  geom_jitter(size=2.5,alpha=0.7, colour='lightgrey')+
  geom_smooth(method='loess', alpha=0.6, span=0.1, se=F, colour='pink')+
  labs(y=expression(paste("Tajima's", italic(" D"))))+
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Scaffold (MB)", title="California")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=13,face="bold"))+
  theme_classic()+
  theme(legend.position = 'none')+
  theme(axis.text.y=element_text(size=13))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
  geom_hline(yintercept = 2, linetype='dashed', colour='red', size=1)+
  geom_hline(yintercept = -2, linetype='dashed', colour='red', size=1)



### Import Arizona data.
TajimasDArizona<-read_table("Arizona.TajimasD")
TajimasDArizona$CHROM<-gsub("chr", "", TajimasDArizona$CHROM)
# Remove windows with poor SNP representation.
TajimasDArizona<-TajimasDArizona %>%
  filter(N_SNPS > 100)

TajimasDArizona$TajimaD<-as.double(TajimasDArizona$TajimaD)
# Mean D in Arizona: 0.5286888. Again, fairly neutral expectation met.
mean(TajimasDArizona$TajimaD)

# Hisotgram shows normal distribution of Tajima's D values.
HistogramTajimasDArizona<-TajimasDArizona %>%
  ggplot(aes(x=TajimaD))+
  geom_histogram(fill='grey', colour='black')+
  theme_pubr()+
  labs(x=expression(paste("Tajima's D - Arizona", italic(" D"))),
       y="count (10kb windows)")+
  geom_vline(xintercept = 0, linetype='dashed', colour='black')

# Genome-wide plot of Tajima's D values.
ArizonanTajimasDLandscape<-TajimasDArizona %>%
  na.omit() %>%
  ggplot(aes(x=BIN_START/1000000, y=TajimaD))+
  geom_jitter(size=2.5,alpha=0.7, colour='lightgrey')+
  geom_smooth(method='loess', alpha=0.6, span=0.1, se=F, colour='orange')+
  labs(y=expression(paste("Tajima's", italic(" D"))))+
  facet_grid(cols = vars(as.numeric(CHROM)),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Scaffold (MB)", title="Arizona")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=13,face="bold"))+
  theme_classic()+
  theme(legend.position = 'none')+
  theme(axis.text.y=element_text(size=13))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.x = element_blank())+
  geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
  geom_hline(yintercept = 2, linetype='dashed', colour='red', size=1)+
  geom_hline(yintercept = -2, linetype='dashed', colour='red', size=1)


MainlandTajDPlots<-ggarrange(ArizonanTajimasDLandscape, CalifornianTajimasDLandscape, nrow=2)
  
png("MainlandTajDPlots.png", units="in", width=8, height=6, res=300)
MainlandTajDPlots
dev.off()


# Ok now I need to find coordinates for significantly DE genes in the genome and then overlap them with Tajima's D.
library(GenomicRanges)
library(tidyverse)
library(ggpubr)

# DEG positions. 
genes_of_interest <- read.delim("EvidenceOfSelectiveSweeps.SharedAcrossTests.genes.gff", header = FALSE, comment.char = "#")
colnames(genes_of_interest) <- c("chr", "source", "feature", "start", "end", "score", "strand", "phase", "attribute")
interest_gr <- GRanges(seqnames = genes_of_interest$chr,
                       ranges = IRanges(start = genes_of_interest$start, end = genes_of_interest$end))

# Now convert Tajimas D data for California and Arziona to the same sort of data.
# Arizona first.
TajimasDArizona<-TajimasDArizona %>%
  mutate(BIN_END=BIN_START + 9999)

TajD_Arizona_gr <- GRanges(seqnames = TajimasDArizona$CHROM,
                         ranges = IRanges(start = TajimasDArizona$BIN_START, end = TajimasDArizona$BIN_END),
                         D = TajimasDArizona$TajimaD)

TajD_Arizona_gr_dataset <- as.data.frame(TajD_Arizona_gr)

# Now process California.
# Arizona first.
TajimasDCalifornia<-TajimasDCalifornia %>%
  mutate(BIN_END=BIN_START + 9999)

TajD_Cali_gr <- GRanges(seqnames = TajimasDCalifornia$CHROM,
                           ranges = IRanges(start = TajimasDCalifornia$BIN_START, end = TajimasDCalifornia$BIN_END),
                           D = TajimasDCalifornia$TajimaD)

TajD_Cali_gr_dataset <- as.data.frame(TajD_Cali_gr)


# Combine the datasets.
DE_Arizona <- findOverlaps(interest_gr, TajD_Arizona_gr)
DE_California <- findOverlaps(interest_gr, TajD_Cali_gr)

# California match overlap and add new DE column for t.test.
Taj_overlaps_California <- findOverlaps(TajD_Cali_gr, interest_gr)
TajD_Cali_gr_dataset$DE_genes <- FALSE
TajD_Cali_gr_dataset$DE_genes[queryHits(Taj_overlaps_California)] <- TRUE

# Arizona match overlap
Taj_overlaps_Arizona <- findOverlaps(TajD_Arizona_gr, interest_gr)
TajD_Arizona_gr_dataset$DE_genes <- FALSE
TajD_Arizona_gr_dataset$DE_genes[queryHits(Taj_overlaps_Arizona)] <- TRUE

t.test(D ~ DE_genes, data = TajD_Cali_gr_dataset, var.equal = F)
t.test(D ~ DE_genes, data = TajD_Arizona_gr_dataset, var.equal = F)


# how many genes show low genetic diversity.
DE_Arizona_TajD_df %>%
  filter(TajD < -1) 

DE_California_TajD_df %>%
  filter(TajD < -1)

# how many genes show high intermediate haplotypes.
DE_Arizona_TajD_df %>%
  filter(TajD > 1) 

DE_California_TajD_df %>%
  filter(TajD > 1)





#### --------- Measuring nucleotide diversity 
# ok now test for 
NucleotideDiversity<-read.table("Altogether_PI_pixy.txt", header=T)

NucleotideDiversity_California<-NucleotideDiversity %>%
  filter(pop == "California")

NucleotideDiversity_Arizona<-NucleotideDiversity %>%
  filter(pop == "Arizona")

# California first.

Pi_California_gr <- GRanges(seqnames = NucleotideDiversity_California$chromosome,
                           ranges = IRanges(start = NucleotideDiversity_California$window_pos_1, end = NucleotideDiversity_California$window_pos_2),
                           Pi = NucleotideDiversity_California$avg_pi)

Pi_California_gr_dataset <- as.data.frame(Pi_California_gr)

# Arizona first.

Pi_Arizona_gr <- GRanges(seqnames = NucleotideDiversity_Arizona$chromosome,
                            ranges = IRanges(start = NucleotideDiversity_Arizona$window_pos_1, end = NucleotideDiversity_Arizona$window_pos_2),
                            Pi = NucleotideDiversity_Arizona$avg_pi)

Pi_Arizona_gr_dataset <- as.data.frame(Pi_Arizona_gr)


# Combine the datasets.
Pi_DE_Arizona <- findOverlaps(interest_gr, Pi_Arizona_gr)
Pi_DE_Cali <- findOverlaps(interest_gr, Pi_California_gr)

# California match overlap and add new DE column for t.test.
overlaps_California <- findOverlaps(Pi_California_gr, interest_gr)
Pi_California_gr_dataset$DE_genes <- FALSE
Pi_California_gr_dataset$DE_genes[queryHits(overlaps_California)] <- TRUE

# Arizona match overlap
overlaps_Arizona <- findOverlaps(Pi_Arizona_gr, interest_gr)
Pi_Arizona_gr_dataset$DE_genes <- FALSE
Pi_Arizona_gr_dataset$DE_genes[queryHits(overlaps_Arizona)] <- TRUE

t.test(Pi ~ DE_genes, data = Pi_California_gr_dataset, var.equal = F)
t.test(Pi ~ DE_genes, data = Pi_Arizona_gr_dataset, var.equal = F)


### now do the same with 
