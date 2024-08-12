library(stringr)
library(Seurat)
library(ggplot2)
library(data.table)
library(tidyr)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(dplyr)

setwd("/data/Ross/lucy_rbp6_analysis/paper")

pcf_obj <- readRDS("emma_freshFrozen_927/PCF_intergrated.rds")

pcf_obj <- subset(pcf_obj, sample == "PCF fresh")

#remove "Tbrucei-" prefix
pcf_mtx <- as.matrix(pcf_obj@assays$RNA@counts)
pcf_meta <- pcf_obj@meta.data

row.names(pcf_mtx) <- str_replace_all(row.names(pcf_mtx), "Tbrucei-", "")


pcf_obj <- CreateSeuratObject(pcf_mtx, meta.data = pcf_meta)

scale_factor <- median(pcf_obj$nCount_RNA)
pcf_obj <- NormalizeData(pcf_obj, scale.factor = scale_factor)

all927 <- read.csv("all927_genes.csv")

#remove genes ID which occur multiple times

single_genes <- names(table(all927$Gene.ID))[which(table(all927$Gene.ID) == 1)]

all927 <- subset(all927, Gene.ID %in% single_genes)

row.names(all927) <- all927$Gene.ID

ribosomal_genes <- all927$Gene.ID[str_detect(all927$Product.Description, "ribosomal protein|ribosomal subunit protein")]
ribosomal_genes <- intersect(ribosomal_genes, row.names(pcf_obj))

#Find ribosomal proteins that are expressed in more than 10% of cells
exp_mtx <- pcf_obj@assays$RNA@layers$counts
row.names(exp_mtx) <- row.names(pcf_obj)

riboProtein_percent <- rowSums(exp_mtx[ ribosomal_genes,] > 0) / dim(pcf_obj)[2]

riboProtein_percent <- riboProtein_percent[which(riboProtein_percent >= 0.1)]

pcf_obj <- MetaFeature(pcf_obj, features = names(riboProtein_percent), meta.name = "RiboProtein")

Cell_cycle_regulated_genes <- read.delim("Cell_cycle_regulated_genes.text")

#Cell cycle analysis
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@layers$counts
  row.names(counts) <- row.names(object)
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

### part one ### 

cellcyclemarkers <- c(Cell_cycle_regulated_genes$Early.G1, Cell_cycle_regulated_genes$Late.G1, Cell_cycle_regulated_genes$S.phase, Cell_cycle_regulated_genes$G2.M.phase) 
head(cellcyclemarkers)

cellcyclemarkers <- intersect(cellcyclemarkers, row.names(pcf_obj))

## calculate the percentage of cells expressing each marker genes
cell_prct <- PrctCellExpringGene(pcf_obj, genes = cellcyclemarkers, group.by = "all")

genes_10prct <- subset(cell_prct, subset = cell_prct$Cell_proportion > 0.05)
genes_10 <- genes_10prct$Markers
head(genes_10)

# Get list of marker genes present in at least 10% cells for each phase
s.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$S.phase %in% genes_10)
s.genes <- s.genes$S.phase

g2m.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$G2.M.phase %in% genes_10)
g2m.genes <- g2m.genes$G2.M.phase

g1e.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$Early.G1 %in% genes_10)
g1e.genes <- g1e.genes$Early.G1

g1l.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$Late.G1 %in% genes_10)
g1l.genes <- g1l.genes$Late.G1

#Old cell cycle analysis
# Calculate an expression score for each phase and save it to the seurat object
pcf_obj <- MetaFeature(pcf_obj, features = s.genes, meta.name = "S.aggregate")
pcf_obj <- MetaFeature(pcf_obj, features = g2m.genes, meta.name = "G2M.aggregate")
pcf_obj <- MetaFeature(pcf_obj, features = g1e.genes, meta.name = "G1e.aggregate")
pcf_obj <- MetaFeature(pcf_obj, features = g1l.genes, meta.name = "G1l.aggregate")


# Create and dataframe with the expression score of each cell and each phase
dfLog <- data.frame(pcf_obj@meta.data[["S.aggregate"]], 
                    pcf_obj@meta.data[["G2M.aggregate"]], 
                    pcf_obj@meta.data[["G1e.aggregate"]],
                    pcf_obj@meta.data[["G1l.aggregate"]])

colnames(dfLog) <- c("S", "G2M", "G1e", "G1l")
rownames(dfLog) <- colnames(pcf_obj)
# Find ratio between the score and average
dfLog$S.ratio <- dfLog$S / mean(dfLog$S)
dfLog$G2M.ratio <- dfLog$G2M / mean(dfLog$G2M)
dfLog$G1e.ratio <- dfLog$G1e / mean(dfLog$G1e)
dfLog$G1l.ratio <- dfLog$G1l / mean(dfLog$G1l)

#to decide on the FC value, Emma suggests plotting the ratios as histograms first

p1 <- ggplot(dfLog, aes(x=G1e.ratio)) + 
  geom_histogram(binwidth=0.01, fill = "#ABD857", colour = "#ABD857", alpha=0.8) +
  geom_vline(xintercept = 1.05, colour = "red", linetype="dashed") + theme_bw() + xlim(0, 2) + ylim(0, 100)

pdf("RBP6_analysis/plots/emmaPCF_G1e_histogram.pdf")
p1
dev.off()

p2 <- ggplot(dfLog, aes(x=G1l.ratio)) + 
  geom_histogram(binwidth=0.01, fill = "#ABD857", colour = "#ABD857", alpha=0.8) +
  geom_vline(xintercept = 0.9, colour = "red", linetype="dashed") + theme_bw() + xlim(0, 2) + ylim(0, 100)

pdf("RBP6_analysis/plots/emmaPCF_G1l_histogram.pdf")
p2
dev.off()

p3 <- ggplot(dfLog, aes(x=S.ratio)) + 
  geom_histogram(binwidth=0.01, fill = "#F9E855", colour = "#F9E855", alpha=0.8) + 
  geom_vline(xintercept = 1, colour = "red", linetype="dashed") + theme_bw() + xlim(0, 2) + ylim(0, 100)

pdf("RBP6_analysis/plots/emmaPCF_S_histogram.pdf")
p3
dev.off()

p4 <- ggplot(dfLog, aes(x=G2M.ratio, )) + 
  geom_histogram(binwidth=0.01, fill = "#3D0750", colour = "#3D0750", alpha=0.8) +
  geom_vline(xintercept = 1.05, colour = "red", linetype="dashed") + theme_bw()  + ylim(0, 100)

pdf("RBP6_analysis/plots/emmaPCF_G2M_histogram.pdf")
p4
dev.off()

# Find the top scoring phase of each cell, with FC > 1.05 
assignmentsLog <- apply(
  X = dfLog[, 5:8],
  MARGIN = 1,
  FUN = function(scores, first = 'S', second = 'G2M', third = "G1e", fourth = "G1l", null = 'Unlabelled') {
    if (all(scores < 1)) {
      return(null)
    } else {
      if (length(which(x = scores == max(scores))) > 1.05) {
        return('Undecided')
      } else {
        return(c(first, second, third, fourth)[which(x = scores == max(scores))])
      }
    }
  }    
)

dfLog$Phase <- assignmentsLog
dfLog$Cluster <- pcf_obj@active.ident

pcf_obj$Phase <- dfLog$Phase


saveRDS(pcf_obj,"/data/Ross/lucy_rbp6_analysis/paper/integration_analysis/objects/briggs_processed.rds")


