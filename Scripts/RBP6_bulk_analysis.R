library(phateR)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(phateR)
library(dplyr)
library(Seurat)
library(rgl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(viridis)
library(scran)
library(RColorBrewer)
library(SingleCellExperiment)
library(reshape2)
library(slingshot)
library(stats)
library(stringi)
library(harmony)
library(stringr)
library(gplots)
library(EnhancedVolcano)
library(DESeq2)

library(unixtools)


files_all <- list.files("/data/Ross/lucy_rbp6_analysis/paper/bulk_zikova/bulk_input")

day <- c()
replicate <- c()

for (i in files_all){
  
  split_i <- strsplit(i, "_")[[1]]
  
  day <- append(day, split_i[12])
  replicate <- append(replicate, split_i[13])
  
}


day_replicate <- paste0(day, "_", replicate)

first <- read.delim("/data/Ross/lucy_rbp6_analysis/paper/bulk_zikova/bulk_input/GSM4185419_cutadapt_imb_butter_2017_03_falk_rnaseq_trypanosomes_tc_1_0d_1_S1_R1_001.readcounts.tsv", 
                    sep ="\t", header = F)

first$V1 <- str_replace_all(first$V1, "-t26_1", "")

#Convert 427 genes to 927 orthologs
orthologs <- read.csv("927_to_427_syntenic_orthologs.csv", header = T)

#Remove genes with multiple inputs
genes_keep <- row.names(orthologs)[str_detect(orthologs$Input.Ortholog.s., ",", negate = T)]
orthologs <- orthologs[genes_keep,]
orthologs <- orthologs[, c("Gene.ID", "Input.Ortholog.s.")]

#Remove genes which map to more than one Gene ID
single_map_genes <- names(table(orthologs$Input.Ortholog.s.))[ which(table(orthologs$Input.Ortholog.s.) == 1) ]
orthologs <- subset(orthologs, Input.Ortholog.s. %in% single_map_genes)

#Filter orthologs to remove genes that aren't in the PCF object
orthologs <- subset(orthologs,Input.Ortholog.s. %in% first$V1 )

first <- subset(first, V1 %in% orthologs$Input.Ortholog.s.)

orthologs_ordered <- orthologs[order(match(orthologs$Input.Ortholog.s., first$V1)), ]

first <- first[,-1]
names(first)


exp_mtx <- matrix(0, ncol = length(day_replicate), nrow = length(first))

path <- "/data/Ross/lucy_rbp6_analysis/paper/bulk_zikova/bulk_input/"

for(i in 1:length(day_replicate)){
  x <- read.delim(paste0(path,files_all[i]), sep = "\t", header = F)
  
  x$V1 <- str_replace_all(x$V1, "-t26_1", "")
  
  x <- subset(x, V1 %in% orthologs$Input.Ortholog.s.)
  
  orthologs_ordered <- orthologs[order(match(orthologs$Input.Ortholog.s., x$V1)), ]
  
  exp_mtx[,i] <- as.numeric(x$V2)
}

colnames(exp_mtx) <- day_replicate
row.names(exp_mtx) <- orthologs_ordered$Gene.ID

metaData_obj <- data.frame("day" = day,
                           "dayReplicate" = day_replicate)

metaData_obj$day <- factor(metaData_obj$day, levels = c("0d", "2d", "3d",
                                                                      "4d", "6d", "8d"))

metaData_obj$dayReplicate <- factor(metaData_obj$dayReplicate, levels = c(day_replicate))


library(DESeq2)  

dds_cruzi <- DESeqDataSetFromMatrix(exp_mtx,
                                    colData = metaData_obj,
                                    design = ~day)

smallestGroupSize <- 4
keep <- rowSums(counts(dds_cruzi) >= 10) >= smallestGroupSize
dds_cruzi <- dds_cruzi[keep,]

dds_cruzi <- DESeq(dds_cruzi)
res <- results(dds_cruzi)

vsdata <- vst(dds_cruzi, blind=FALSE)

#trans-sialidase, group V
pdf("/datastore/Ross/cruzi_paper/bulk/plots/BulkRNA_transSialidse_groupV_plotCounts_allSamples.pdf",
    width = 8)
plotCounts(dds_cruzi, gene="Tb927.6.510", intgroup="day") 
dev.off()

#trans-sialidase, group IV
pdf("/datastore/Ross/cruzi_paper/bulk/plots/BulkRNA_transSialidse_groupIV_plotCounts_allSamples.pdf",
    width = 8.5)
plotCounts(dds_cruzi, gene="C4B63_75g82", intgroup="shortSample") 
dev.off()

plotCounts(dds_cruzi, gene="C4B63_11g29", intgroup="stage_fine") 
plotCounts(dds_cruzi, gene="C4B63_113g41", intgroup="stage_fine") 


pdf("/datastore/Ross/cruzi_paper/bulk/plots/BulkRNA_deseq2_PCA_stages.pdf")
plotPCA(vsdata, intgroup="day")
dev.off()







