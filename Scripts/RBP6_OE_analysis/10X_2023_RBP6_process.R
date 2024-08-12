#### Load libraries ####
library(stringr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(dplyr)
library(phateR)
library(clustree)
#### load data from umi_tools count ####

#This is not mapped against VSGS

setwd("/data/Ross/lucy_rbp6_analysis/paper")

Cell_cycle_regulated_genes <- read.delim("Cell_cycle_regulated_genes.text")

cellranger_obj <- Read10X("RBP6_analysis/objects/RBP6_mVSG_filtered_feature_bc_matrix")


VSG_list <- read.delim("VSG_list", header = F)
VSG_list$V1 <- str_replace_all(VSG_list$V1, "Tbrucei---", "")

rrna_genes <- read.delim("ribosomal_RNA.txt", header = F)
rrna_genes$V1 <- str_replace(rrna_genes$V1, "_", "-")

tryp_obj<- CreateSeuratObject(as.matrix(cellranger_obj), min.cells = 10, project = "lucy")

mt_genes <- row.names(tryp_obj)[str_detect(row.names(tryp_obj), "^Tb", negate = TRUE)]
mt_genes <- mt_genes[str_detect(mt_genes, "tmp", negate = TRUE)]
mt_genes <- mt_genes[str_detect(mt_genes, "rRNA", negate = TRUE)]

rrna_genes <- intersect(rrna_genes$V1, row.names(tryp_obj))


tryp_obj[["percent.rRNA"]] <- PercentageFeatureSet(tryp_obj, features = rrna_genes)
tryp_obj[["percent.mt"]] <- PercentageFeatureSet(tryp_obj, features = mt_genes)

#Print Median values of percent rRNA and kDNA
median(tryp_obj$percent.mt)
median(tryp_obj$percent.rRNA)

pdf("plots/regressUMI/percent_mt_outlier.pdf")
VlnPlot(tryp_obj, features = "percent.mt") 
dev.off()

# tryp_obj <- subset(tryp_obj, percent.mt < 20)

pdf("RBP6_analysis/plots/violin_nFeature_RNA.pdf")
VlnPlot(tryp_obj, features = "nFeature_RNA") + 
  geom_hline(aes(yintercept = 400)) +
  geom_hline(aes(yintercept = 2500))
dev.off()

pdf("plots/regressUMI/violin_percent_mt.pdf")
VlnPlot(tryp_obj, features = "percent.mt") + geom_hline(aes(yintercept = 5))
dev.off()

pdf("plots/regressUMI/violin_percent_rRNA.pdf")
VlnPlot(tryp_obj, features = "percent.rRNA") + geom_hline(aes(yintercept = 4.5))
dev.off()

pdf("plots/regressUMI/scatter_nFeature_nCount.pdf")
FeatureScatter(tryp_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(aes(yintercept = 400)) +
  geom_hline(aes(yintercept = 2500))
  
dev.off()

FeatureScatter(tryp_obj, feature1 = "nCount_RNA", feature2 = "percent.rRNA")
FeatureScatter(tryp_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(tryp_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

tryp_obj <- subset(tryp_obj, subset = nFeature_RNA > 400 & nFeature_RNA < 2500 & percent.rRNA < 4.5 & percent.mt < 5)

scale_factor <- median(tryp_obj$nCount_RNA)

tryp_obj <- NormalizeData(tryp_obj, scale.factor = scale_factor)

tryp_obj <- FindVariableFeatures(tryp_obj, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10_var <- head(VariableFeatures(tryp_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(tryp_obj)

pdf("RBP6_analysis/plots/var_feature_plot_intersect.pdf")
LabelPoints(plot = plot1, points = top10_var, repel = TRUE)
dev.off()

tryp_obj <- ScaleData(tryp_obj, vars.to.regress = "nCount_RNA")

tryp_obj <- RunPCA(tryp_obj, npcs = 50)

pdf("plots/regressUMI/pca_UMIregress_plot.pdf")
DimPlot(tryp_obj)
dev.off()


pdf("plots/regressUMI/elbowPlot_regressUMI.pdf")
ElbowPlot(tryp_obj, ndims = 50)
dev.off()

tryp_obj <- FindNeighbors(tryp_obj, dims = 1:12)
tryp_obj <- FindClusters(tryp_obj, resolution = c(0.1, 0.2, 0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

pdf("plots/regressUMI/clustree_regressUMI.pdf", height = 14)
clustree(tryp_obj, prefix = "RNA_snn_res.")
dev.off()

tryp_obj <- FindClusters(tryp_obj, resolution = 0.6)


tryp_obj <- RunUMAP(tryp_obj, dims = 1:12)

pdf("plots/regressUMI/UMAP_regressUMI_plotClusters.pdf")
DimPlot(tryp_obj, label = T, group.by = "seurat_clusters")
dev.off()

saveRDS(tryp_obj, "RBP6_analysis/objects/2023_10X_sample_QC_mVSG.rds")

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

## calculate the percentage of cells expressing each marker genes
cell_prct <- PrctCellExpringGene(tryp_obj, genes = cellcyclemarkers, group.by = "all")

genes_10prct <- subset(cell_prct, subset = cell_prct$Cell_proportion > 0.1)
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
tryp_obj <- MetaFeature(tryp_obj, features = s.genes, meta.name = "S.aggregate")
tryp_obj <- MetaFeature(tryp_obj, features = g2m.genes, meta.name = "G2M.aggregate")
tryp_obj <- MetaFeature(tryp_obj, features = g1e.genes, meta.name = "G1e.aggregate")
tryp_obj <- MetaFeature(tryp_obj, features = g1l.genes, meta.name = "G1l.aggregate")


# Create and dataframe with the expression score of each cell and each phase
dfLog <- data.frame(tryp_obj@meta.data[["S.aggregate"]], 
                    tryp_obj@meta.data[["G2M.aggregate"]], 
                    tryp_obj@meta.data[["G1e.aggregate"]],
                    tryp_obj@meta.data[["G1l.aggregate"]])

colnames(dfLog) <- c("S", "G2M", "G1e", "G1l")
rownames(dfLog) <- colnames(tryp_obj)
# Find ratio between the score and average
dfLog$S.ratio <- dfLog$S / mean(dfLog$S)
dfLog$G2M.ratio <- dfLog$G2M / mean(dfLog$G2M)
dfLog$G1e.ratio <- dfLog$G1e / mean(dfLog$G1e)
dfLog$G1l.ratio <- dfLog$G1l / mean(dfLog$G1l)

#to decide on the FC value, Emma suggests plotting the ratios as histograms first

p1 <- ggplot(dfLog, aes(x=G1e.ratio)) + 
  geom_histogram(binwidth=0.01, fill = "#ABD857", colour = "#ABD857", alpha=0.8) +
  geom_vline(xintercept = 1, colour = "red", linetype="dashed") + theme_bw() + xlim(0, max(dfLog$G1e.ratio)) + ylim(0, 250)

pdf("plots/regressUMI/G1e_histogram.pdf")
p1
dev.off()

p2 <- ggplot(dfLog, aes(x=G1l.ratio)) + 
  geom_histogram(binwidth=0.01, fill = "#ABD857", colour = "#ABD857", alpha=0.8) +
  geom_vline(xintercept = 0.90, colour = "red", linetype="dashed") + theme_bw() + xlim(0, max(dfLog$G1l.ratio)) + ylim(0, 250)

pdf("plots/regressUMI/G1l_histogram.pdf")
p2
dev.off()

p3 <- ggplot(dfLog, aes(x=S.ratio)) + 
  geom_histogram(binwidth=0.01, fill = "#F9E855", colour = "#F9E855", alpha=0.8) + 
  geom_vline(xintercept = 0.85, colour = "red", linetype="dashed") + theme_bw() + xlim(0, max(dfLog$S.ratio)) + ylim(0, 250)

pdf("RBP6_analysis/plots/S_histogram.pdf")
p3
dev.off()

p4 <- ggplot(dfLog, aes(x=G2M.ratio, )) + 
  geom_histogram(binwidth=0.01, fill = "#3D0750", colour = "#3D0750", alpha=0.8) +
  geom_vline(xintercept = 1, colour = "red", linetype="dashed") + theme_bw() + xlim(0, max(dfLog$G2M.ratio))  + ylim(0, 250)

pdf("plots/regressUMI/G2M_histogram.pdf")
p4
dev.off()

# Find the top scoring phase of each cell, with FC > 1.05 
assignmentsLog <- apply(
  X = dfLog[, 5:8],
  MARGIN = 1,
  FUN = function(scores, first = 'S', second = 'G2M', third = "G1e", fourth = "G1l", null = 'Unlabelled') {
    if (all(scores < 0.9)) {
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
dfLog$Cluster <- tryp_obj@active.ident

tryp_obj$Phase <- dfLog$Phase

pdf("plots/regressUMI/UMAP_regressUMI_plotCellCyclePhase.pdf")
DimPlot(tryp_obj, group.by = "Phase")
dev.off()

pdf("plots/regressUMI/PCA_regressUMI_plotCellCyclePhase.pdf")
DimPlot(tryp_obj, group.by = "Phase", reduction = "pca")
dev.off()

library(RColorBrewer)
pt <- table(tryp_obj$Phase, Idents(tryp_obj))
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/regressUMI/proportionPlot_regressUMI_plotCellCyclePhase.pdf")
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()

saveRDS(tryp_obj, "RBP6_analysis/objects/2023_10X_sample_phase_QC_mVSG.rds")



