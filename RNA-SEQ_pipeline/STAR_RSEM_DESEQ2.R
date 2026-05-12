# Load packages ----- 
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
#BiocManager::install("DESeq2")
library(DESeq2)
library(gt)
library(DT)
library(plotly)
library(ashr)
library(devtools)
#BiocManager::install("rhdf5")
library(rhdf5)
#devtools::install_github("pachterlab/sleuth")
#install.packages("tximport")
library(babelgene)
library(sleuth)
library(limma)
library(clusterProfiler)
library(ggrepel)
# BiocManager::install("tximport")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# library(BiocManager)
BiocManager::install("tximport")
library(tximport)
# the essentials ----
# this chunk contains the minimal essential code from this script. Simply uncomment the lines below and run the code.
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot)# allows you to combine multiple plots in one figure
library(ggpubr)
library(sva)


targets <- read_csv("sample_metadata_RNAexperiment.csv")# read in your study design

path <- file.path("RSEM",targets$sample, paste0(targets$sample, ".genes.results.gz")) # set file paths to your mapped data

all(file.exists(path))

Txi_gene <- tximport(path, 
                     type = "rsem", 
                     txIn = FALSE,
                     txOut = FALSE)
Txi_gene$length[Txi_gene$length <= 0] <- 1

dds <- DESeqDataSetFromTximport(Txi_gene, targets, ~condition)
# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)

# capture sample labels from the study design file that you worked with and saved as 'targets' in step 1
targets
sampleLabels <- targets$sample
##### ---------------- OLD SCRIPT -- IGNORE --------------
# # Generate summary stats for your data ----
# # 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# # then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
# # then we use the 'rowSds', 'rowMeans' and 'rowMedians' functions from the matrixStats package
# myTPM.stats <- transform(myTPM, 
#                          SD=rowSds(myTPM), 
#                          AVG=rowMeans(myTPM),
#                          MED=rowMedians(myTPM))
# 
# 
# # Create your first plot using ggplot2 ----
# # produce a scatter plot of the transformed data
# ggplot(myTPM.stats) + 
#   aes(x = SD, y = MED) +
#   geom_point(shape=25, size=3)
# 
# 
# # Make a DGElist from your counts, and plot ----
# myDGEList <- DGEList(myCounts)
# # take a look at the DGEList object 
# myDGEList
# #DEGList objects are a good R data file to consider saving to you working directory
# save(myDGEList, file = "myDGEList_RSEM")
# #Saved DGEList objects can be easily shared and loaded into an R environment
# load(file = "myDGEList_RSEM")
# 
# # use the 'cpm' function from EdgeR to get counts per million
# cpm <- cpm(myDGEList) 
# colSums(cpm)
# log2.cpm <- cpm(myDGEList, log=TRUE)
# 
# # 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
# log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
# log2.cpm.df
# # add your sample names to this dataframe (we lost these when we read our data in with tximport)
# #colnames(log2.cpm.df) <- c("gene_name", sampleLabels)
# #log2.cpm.df$geneID<-gsub("OrmOch_", "", log2.cpm.df$geneID)
# 
# #first, take a look at how many genes or transcripts have no read counts at all
# table(rowSums(myDGEList$counts==0)==16)
# keepers <- rowSums(cpm>1)>=7
# myDGEList.filtered <- myDGEList[keepers,]
# dim(myDGEList.filtered)
# 
# log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
# log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
# colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
# 
# log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
#                                            cols = AS2:DC8, # column names to be stored as a SINGLE variable
#                                            names_to = "samples", # name of that new variable (column)
#                                            values_to = "expression") # name of new variable (column) storing all the values (data)
# 
# ggplot(log2.cpm.filtered.df.pivot) +
#   aes(x=samples, y=expression, fill=samples) +
#   geom_violin(trim = FALSE, show.legend = FALSE) +
#   stat_summary(fun = "median", 
#                geom = "point", 
#                shape = 95, 
#                size = 10, 
#                color = "black", 
#                show.legend = FALSE) +
#   labs(y="log2 expression", x = "sample",
#        title="Log2 Counts per Million (CPM)",
#        subtitle="filtered, non-normalized",
#        caption=paste0("produced on ", Sys.time())) +
#   theme_bw()
# 
# 
# # check the PCA before performing normalisation.
# log2.cpm.matrix <- as.matrix(log2.cpm.filtered.df[, -1])
# log2.cpm.matrix <- apply(log2.cpm.matrix, 2, as.numeric)
# 
# # PCA
# pca.res <- prcomp(t(log2.cpm.matrix), scale. = FALSE, retx = TRUE)
# pc.var <- pca.res$sdev^2
# pc.per <- round(pc.var / sum(pc.var) * 100, 1)
# 
# # Combine PCA results with metadata
# pca.res.df <- as_tibble(pca.res$x)
# group <- targets$condition
# pca.res.df$Group <- group
# pca.res.df$Sample <- sample_names
# 
# ggplot(pca.res.df, aes(x = PC1, y = PC2, color = Group, label = paste(Sample, Group, sep = " - "))) +
#   geom_point(size = 4) +
#   stat_ellipse() +
#   geom_text_repel(size = 3, max.overlaps = 20) +
#   xlab(paste0("PC1 (", pc.per[1], "%)")) + 
#   ylab(paste0("PC2 (", pc.per[2], "%)")) +
#   labs(title = "PCA Plot",
#        caption = paste0("produced on ", Sys.time())) +
#   coord_fixed() +
#   theme_bw()
# 
# 
# 
# # Normalize your data ----
# myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
# # take a look at this new DGEList object...how has it changed?
# 
# # use the 'cpm' function from EdgeR to get counts per million from your normalized data
# log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
# log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
# colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
# 
# #rownames(countDataMatrix) <- log2.cpm.filtered.norm.df[ , 1]
# 
# # PCA ----
# library(tidyverse)
# library(DT)
# library(gt)
# library(plotly)
# library(ggrepel)
# 
# group <- targets$condition
# sample_names<-targets$sample
# group <- factor(group)
# 
# # Remove the gene ID column and convert to a numeric matrix
# log2.cpm.matrix <- as.matrix(log2.cpm.filtered.norm.df[, -1])
# log2.cpm.matrix <- apply(log2.cpm.matrix, 2, as.numeric)
# 
# # PCA
# pca.res <- prcomp(t(log2.cpm.matrix), scale. = FALSE, retx = TRUE)
# pc.var <- pca.res$sdev^2
# pc.per <- round(pc.var / sum(pc.var) * 100, 1)
# 
# # Combine PCA results with metadata
# pca.res.df <- as_tibble(pca.res$x)
# pca.res.df$Group <- group
# pca.res.df$Sample <- sample_names
# 
# pca.plot <- ggplot(pca.res.df, aes(x = PC1, y = PC2, color = Group, label = paste(Sample, Group, sep = " - "))) +
#   geom_point(size = 4) +
#   stat_ellipse() +
#   geom_text_repel(size = 3, max.overlaps = 20) +
#   xlab(paste0("PC1 (", pc.per[1], "%)")) + 
#   ylab(paste0("PC2 (", pc.per[2], "%)")) +
#   labs(title = "PCA Plot",
#        caption = paste0("produced on ", Sys.time())) +
#   coord_fixed() +
#   theme_bw()
# 
# pca.plot
# 
# # Looks like CS3 is a weird sample, might be worth removing it to see whether the analyses changes. 
# colnames(myDGEList.filtered.norm$counts) <- sampleLabels
# rownames(myDGEList.filtered.norm$samples) <- sampleLabels
# 
# samples_to_keep <- colnames(myDGEList.filtered.norm$counts) != "CS3"
# 
# # Subset DGEList and metadata
# myDGEList.filtered.norm <- myDGEList.filtered.norm[, samples_to_keep]
# targets <- targets[samples_to_keep, ]
# 
# ### -
# group <- factor(targets$condition)
# mod <- model.matrix(~ group, data = targets)  # Full model
# mod0 <- model.matrix(~ 1, data = targets)# Null model
# 
# Temp_myDGEList  <- myDGEList.filtered.norm$counts
# idxmyDGEList  <- rowMeans(Temp_myDGEList) > 7
# Temp_dataset_myDGEList  <- Temp_myDGEList[idxmyDGEList, ]
# colnames(Temp_dataset_myDGEList) <- NULL
# 
# # Estimate the number of significant surrogate variables.
# svseq <- svaseq(Temp_dataset_myDGEList, mod, mod0)
# 
# # add them to the design.
# targets$SV1 <- svseq$sv[,1]
# 
# # Rebuild design matrix and proceed with voom/limma
# group <- factor(targets$condition)
# design <- model.matrix(~ 0 + group + SV1, data = targets)
# colnames(design)[1:length(levels(group))] <- levels(group)
# # 
# # #Set up your design matrix ----
# # group <- factor(targets$condition)
# # design <- model.matrix(~0 + group)
# # colnames(design) <- levels(group)
# 
# #Model mean-variance trend and fit linear model to data ----
# # Use VOOM function from Limma package to model the mean-variance relationship
# # Model mean-variance trend
# v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# 
# # Fit model
# fit <- lmFit(v.DEGList.filtered.norm, design)
# 
# # Define contrasts for your comparison of interest
# contrast.matrix <- makeContrasts(Perception = Control - Song, levels = design)
# 
# # Apply contrasts and compute statistics
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# 
# # View top differentially expressed genes
# topTable(fit2, coef = "Perception")
# myTopHits <- topTable(fit2, adjust ="fdr", coef = "Perception", number=10, sort.by="logFC")
# 
# 
# myTopHits.df <- myTopHits %>%
#   as_tibble(rownames = "geneID")
# 
# #gt(myTopHits.df)
# 
# vplot <- ggplot(myTopHits.df) +
#   aes(y=-log10(P.Value), x=logFC, text = paste("Symbol:", geneID)) +
#   geom_point(size=2.5, alpha=0.75) +
#   geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
#   geom_vline(xintercept = 2, linetype="longdash", colour="#BE684D", size=1) +
#   geom_vline(xintercept = -2, linetype="longdash", colour="#2C467A", size=1) +
#   annotate("rect", xmin = 2, xmax = 12, ymin = -log10(0.05), ymax = 7.5, alpha=.2, fill="#BE684D") +
#   annotate("rect", xmin = -2, xmax = -12, ymin = -log10(0.05), ymax = 7.5, alpha=.2, fill="#2C467A") +
#   labs(title="",
#        subtitle = "")+
#   ylab("-log10(p-value)")+
#   theme_pubr()+
#   theme(axis.text = element_text(size=14))
# 
# vplot
# 
# Significant_hits_DE<-myTopHits.df %>%
#   filter(P.Value < 0.05)
# 
# 
# write.csv(Significant_hits_DE, "Significant_hits_DE_LIMMA_0.05Pval.csv", quote = F, row.names = F)
# 

#--------- ------ THIS IS THE NEW SCRIPT ----------------
# Load packages ----- 
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
BiocManager::install("DESeq2")
library(DESeq2)
library(gt)
library(DT)
library(plotly)
library(ashr)
BiocManager::install("sva")
library(sva)

dds <- DESeqDataSetFromTximport(Txi_gene, targets, ~condition)
dds$condition <- relevel(dds$condition, ref = "Control")

# Step 3: Run DESeq -----
dds <- DESeq(dds)

# samples_to_remove <- c("CS2", "CS6")
# 
# dds <- dds[, !colnames(dds) %in% samples_to_remove]

res <- results(dds)

# Leeban has added this in to try get rid of the other effects.

Temp_dataset  <- counts(dds, normalized = T)
idx  <- rowMeans(Temp_dataset) > 7
Temp_dataset  <- Temp_dataset[idx, ]
mod  <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~  1, colData(dds))
svseq <- svaseq(Temp_dataset, mod, mod0) # now I have some surrogate variables which I can regress out.
n.sv <- num.sv(Temp_dataset, mod, method="leek")

targets$SV1 <- svseq$sv[,1]
targets$SV2 <- svseq$sv[,2]
targets$SV3 <- svseq$sv[,3]
targets$SV4 <- svseq$sv[,4]


targets$Batch<-as.factor(targets$Batch)

designSV <- ~ condition +SV1 +SV2 +SV3+ SV4

DDS_Altered <- DESeqDataSetFromTximport(Txi_gene,
                                                  targets,
                                                  design =  designSV)

# run de seq
DDS_Altered <- DESeq(DDS_Altered)
res_DDS_Altered <- results(DDS_Altered)

# PCA
vsd <- vst(DDS_Altered, blind = FALSE)
plotPCA(vsd, intgroup="SV1")
# Extract PCA data
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Add labels (assuming rownames are sample names)
pcaData$sample <- rownames(pcaData)

# Plot with labels
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +  # adjust vertical position of labels
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

# LogFC Correction 

contrast <- c("condition", "Control", "Song")
ddlfc <- lfcShrink(DDS_Altered, contrast =contrast, res=res_DDS_Altered, type="ashr")
plotMA(ddlfc, ylim=c(-5,5))

# -- summarise results
summary(res_DDS_Altered)

# # set the padj as < 0.01
# res0.01 <- results(dds, alpha = 0.01)
# summary(res0.01)
# sum(res0.01$padj <0.05, na.rm =TRUE)
# # set the padj as < 0.05
# res0.05 <- results(dds, alpha = 0.05)
# summary(res0.05)
# 
# sum(res0.01$padj <0.05, na.rm =TRUE)
# 
# # MA plot
# plotMA(res0.01)
# plotMA(res0.05)
# 
# plotMA(resLFC)
# 
# 
# # plot counts 
# 
# plotCounts(dds, gene=which.min(res$padj), intgroup = "condition")
# 
# order(res0.01$padj)

# export result of DEGs
DEG <- ddlfc %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")



# # filter genes with lfc >2 and p<0.01
# DEG_filtered_2_0.05 <- DEG %>%
#   dplyr::filter((log2FoldChange>=2 | log2FoldChange <= (-2)) & padj <0.05)

# what about significant genes that are not based on some arbitrary threshold. Given that we are talking about neural genes, there's maybe not a huge DE.
DEG_filtered_0.05 <- DEG %>%
  dplyr::filter(padj <0.05)


# Upregulated in female heads when hearing song.
UpregulatedDEGs<-DEG_filtered_0.05 %>%
  filter(log2FoldChange > 0)

# Downregulated in female heads when hearing song.
DownregulatedDEGs<-DEG_filtered_0.05 %>%
  filter(log2FoldChange < 0)

# filter genes with lfc >1 and p<0.05 --- 
DEG_filtered_1_0.05 <- DEG %>%  
  dplyr::filter((log2FoldChange>=1 | log2FoldChange <= (-1)) & padj <0.05)

DEG_filtered_1_0.05_up <- DEG %>%  
  dplyr::filter(log2FoldChange>=1 & padj <0.05)

DEG_filtered_1_0.05_down <- DEG %>%  
  dplyr::filter(log2FoldChange<=-1 & padj <0.05)

# Volcano plot
DEG2 <- na.omit(DEG)

ggplot(data=DEG_filtered_0.05, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.4, size=3) + 
  labs(x="log2 fold change") +
  ylab("-log10 pvalue") +
  ggtitle("DEG between Control and Song") + 
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(size=15, hjust=0.5),) +
  scale_color_manual(values=c('#a121f0','#bebebe', '#ffad21'))
  

plotCounts(DDS_Altered, gene="gene:ENSEGHG00000004131", intgroup="condition")


write.csv(DEG_filtered_2_0.05, "DEG_ORMIA_HEADS_logfc2_PADJ0.05.csv")

write.csv(DEG_filtered_0.05, "DEG_ORMIA_HEADS_PADJ0.05.csv")
write.csv(DEG_filtered_1_0.05, "DEG_filtered_1_0.05.csv")


# Annotations file -- Eggnog mapped proteins for the DEGs mapped to Insecta, then extracted one-to-one orthologs of D.melanogaster only via Flybase.
annotations_file<-read.csv("DEG_ORMIA_HEADS_Annotations.csv")

DEG2$gene_id<-gsub("gene:", "", DEG2$gene_id)
DEG2$gene_id<-sapply(strsplit(DEG2$gene_id,"_"), `[`, 1)
annotations_file$X.query <- gsub("\\.1", "", annotations_file$X.query)
annotations_file$X.query<-gsub("GHP", "GHG", annotations_file$X.query)
# Re-do DEG plot but with annotations. Add annotations to DEG2
DEG2$Annotations<-annotations_file$ASSOCIATED_GENE[match(DEG2$gene_id, annotations_file$X.query)]
DEG2$Interesting<-annotations_file$GO_Interesting[match(DEG2$gene_id, annotations_file$X.query)]



# Add a category column
DEG$colour_group <- "NS"  # default: not significant
DEG$colour_group[DEG$padj < 0.05 & DEG$log2FoldChange < 0] <- "Down"
DEG$colour_group[DEG$padj < 0.05 & DEG$log2FoldChange > 0] <- "Up"

# Updated plot
DEG_plot <- ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), colour=colour_group)) +
  geom_point(alpha=0.3, size=4) +
  labs(x="log2 fold change",
       title="") +
  ylab("-log10 pvalue") +
  geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black', size=1) +
  geom_vline(xintercept = 0, linetype='dashed', colour='black', size=1) +
  theme_pubr(base_size = 20) + 
  theme(plot.title = element_text(size=15, hjust=0.5)) +
  scale_color_manual(values=c('Down'='blue', 'NS'='grey50', 'Up'='red')) +
  theme(legend.position = 'none')


png("DEG_plot.png", units="in", width=5, height=4, res=300)
DEG_plot
dev.off()


# save this and keep 
DEG_filtered_0.05_SV4<- DEG_filtered_0.05
intersect(DEG_filtered_0.05$gene_id, DEG_filtered_0.05_SV4$gene_id)

###------------------ Overrpresentation test between DEG and PBS ---------------------
A_total <- 298
B_total <- 1328
overlap <- 84
background_size <- 11049

# Build contingency table
#        In B   Not in B
# In A     15      298 - 15
# Not A  238-15  rest of universe

a <- overlap
b <- A_total - overlap
c <- B_total - overlap
d <- background_size - (a + b + c)

contingency_table <- matrix(c(a, b, c, d), nrow = 2,
                            dimnames = list("A" = c("In A", "Not in A"),
                                            "B" = c("In B", "Not in B")))
contingency_table

fisher.test(contingency_table, alternative = "greater")  # for overrepresentation


###------------------ Overrpresentation test between DEG and XPEHH ---------------------
A_total <- 355
B_total <- 1327
overlap <- 119
background_size <- 11049

# Build contingency table
#        In B   Not in B
# In A     15      298 - 15
# Not A  238-15  rest of universe

a <- overlap
b <- A_total - overlap
c <- B_total - overlap
d <- background_size - (a + b + c)

contingency_table <- matrix(c(a, b, c, d), nrow = 2,
                            dimnames = list("A" = c("In A", "Not in A"),
                                            "B" = c("In B", "Not in B")))
contingency_table

fisher.test(contingency_table, alternative = "greater")  # for overrepresentation

###------------------ Overrpresentation test between DEG and SweeD ---------------------
A_total <- 183
B_total <- 1327
overlap <- 65
background_size <- 11049

# Build contingency table
#        In B   Not in B
# In A     15      298 - 15
# Not A  238-15  rest of universe

a <- overlap
b <- A_total - overlap
c <- B_total - overlap
d <- background_size - (a + b + c)

contingency_table <- matrix(c(a, b, c, d), nrow = 2,
                            dimnames = list("A" = c("In A", "Not in A"),
                                            "B" = c("In B", "Not in B")))
contingency_table

fisher.test(contingency_table, alternative = "greater")  # for overrepresentation



# export count data and deseq results.
colnames(DDS_Altered) <- targets$sample

count_df <- as.data.frame(counts(DDS_Altered))
count_df$gene <- rownames(count_df)
count_df <- count_df[, c("gene", setdiff(colnames(count_df), "gene"))]

# export the DEG
DEG2

write.csv(DEG2, "DifferentialExpressionDataframe_Oochracea_SongVsNoSongExperiment.csv", quote = F, row.names = F)
write.csv(count_df, "NormalisedCountDataframe_Oochracea_SongVsNoSongExperiment.csv", quote = F, row.names = F)
