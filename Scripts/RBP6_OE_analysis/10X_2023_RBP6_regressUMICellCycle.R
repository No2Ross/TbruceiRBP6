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
library(stringr)
library(smplot2)
#### load data from umi_tools count ####

setwd("/data/Ross/lucy_rbp6_analysis/paper")

all927 <- read.csv("all927_genes.csv")

#remove genes ID which occur multiple times

single_genes <- names(table(all927$Gene.ID))[which(table(all927$Gene.ID) == 1)]

all927 <- subset(all927, Gene.ID %in% single_genes)

row.names(all927) <- all927$Gene.ID

Cell_cycle_regulated_genes <- read.delim("Cell_cycle_regulated_genes.text")

cellcycle_all <- unique(c(Cell_cycle_regulated_genes$Early.G1, Cell_cycle_regulated_genes$Late.G1, Cell_cycle_regulated_genes$S.phase,
                          Cell_cycle_regulated_genes$G2.M.phase))

tryp_obj <- readRDS("RBP6_analysis/objects/2023_10X_sample_phase_QC_mVSG.rds")

#Intersection of variable features and cell cycle genes.
length(intersect(tryp_obj@assays[["RNA"]]@meta.data[["var.features"]][is.na(tryp_obj@assays[["RNA"]]@meta.data[["var.features"]]) == F], 
                 cellcycle_all))

var_features <- tryp_obj@assays[["RNA"]]@meta.data[["var.features"]][is.na(tryp_obj@assays[["RNA"]]@meta.data[["var.features"]]) == F]

var_features <- setdiff(var_features, cellcycle_all)


tryp_obj <- ScaleData(tryp_obj, vars.to.regress = c("nCount_RNA", "G2M.aggregate", 
                                                    "G1e.aggregate", "S.aggregate", "G1l.aggregate"),
                      features = var_features)

tryp_obj <- RunPCA(tryp_obj, npcs = 50, features =var_features)

pdf("plots/regressUMI_CC/elbowPlot_regressUMICC.pdf")
ElbowPlot(tryp_obj, ndims = 50)
dev.off()

tryp_obj <- FindNeighbors(tryp_obj, dims = 1:10)
tryp_obj <- FindClusters(tryp_obj, resolution = c(0.1, 0.2, 0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

pdf("plots/regressUMI_CC/clustree_regressUMICellCycle.pdf", height = 14)
clustree(tryp_obj, prefix = "RNA_snn_res.")
dev.off()

tryp_obj <- FindClusters(tryp_obj, resolution =  0.6)

tryp_obj <- RunUMAP(tryp_obj, dims = 1:10)

pdf("plots/regressUMI_CC/UMAP_regressUMICellCycle_plotClusters.pdf")
DimPlot(tryp_obj, label = T, group.by = "seurat_clusters")
dev.off()

pdf("plots/regressUMI_CC/UMAP_regressUMICellCycle_plotCellCyclePhase.pdf")
DimPlot(tryp_obj, group.by = "Phase")
dev.off()

pdf("plots/regressUMI_CC/PCA_regressUMICellCycle_plotClusters.pdf")
DimPlot(tryp_obj, label = T, group.by = "seurat_clusters", reduction = "pca")
dev.off()

pdf("plots/regressUMI_CC/PCA_regressUMICellCycle_plotPhase.pdf")
DimPlot(tryp_obj, label = T, group.by = "Phase", reduction = "pca")
dev.off()

library(RColorBrewer)
pt <- table(tryp_obj$Phase, Idents(tryp_obj))
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/regressUMI_CC/proportionPlot_regressUMICellCycle_plotCellCyclePhase.pdf")
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()



markers <- FindAllMarkers(tryp_obj, only.pos = T, logfc.threshold = 0.5, min.pct = 0.05)

markers_subset <- subset(markers, p_val_adj < 0.05)

markers_subset$pct_diff <- abs(markers_subset$pct.1 - markers_subset$pct.2)

markers_subset %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  ungroup() -> top15_rbp6

markers_subset %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> top10

DimPlot(tryp_obj, label = T, group.by = "seurat_clusters")

pdf("plots/regressUMI_CC/heatmap_withRibosomalProteins_regressUMICellCycle.pdf")
DoHeatmap(tryp_obj, features = top10$gene)
dev.off()


#Check markers of clusters 0 and 6 to see if any of the markers are not ribosomal proteins
cluster0_markers <- subset(markers_subset, cluster == "0")
cluster0_markers <- subset(cluster0_markers, gene %in% row.names(all927))
cluster0_markers_description <- all927[cluster0_markers$gene,]
percent_ribo_C0 <- sum(str_detect(cluster0_markers_description$Product.Description, "ribosomal protein")) / dim(cluster0_markers_description)[1]

cluster6_markers <- subset(markers_subset, cluster == "6")
cluster6_markers <- subset(cluster6_markers, gene %in% row.names(all927))
cluster6_markers_description <- all927[cluster6_markers$gene,]
percent_ribo_C6 <- sum(str_detect(cluster6_markers_description$Product.Description, "ribosomal protein")) / dim(cluster6_markers_description)[1]

nonribosomal_genes <- all927$Gene.ID[str_detect(all927$Product.Description, "ribosomal protein|ribosomal subunit protein", negate = T)]

#Remove ribosomal protein genes
markers_subset <- subset(markers_subset, gene %in% nonribosomal_genes)

markers_subset$pct_diff <- abs(markers_subset$pct.1 - markers_subset$pct.2)

markers_subset %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> top10


# #Cluster 3 markers
# #Hexokinase 1, GPEET, ZC3H37
# pdf("plots/regressUMI_CC/VlnPlot_C3_markers_regressUMICellCycle.pdf")
# VlnPlot(tryp_obj, features = c("Tb927.10.2010", "Tb927.6.510", "Tb927.10.12780"))
# dev.off()
# 
# #Cluster 2 markers
# #ZC3H36, PAD7, hypothetical protein, conserved
# pdf("plots/regressUMI_CC/VlnPlot_C2_markers_regressUMICellCycle.pdf")
# VlnPlot(tryp_obj, features = c("Tb927.10.12760", "Tb927.7.5990", "Tb927.11.1280"))
# dev.off()
# 
# #ZC3H35, ZC3H36, ZC3H37
# FeaturePlot(tryp_obj, features = c("Tb927.10.12740", "Tb927.10.12760", "Tb927.10.12780"))
# 
# #Cluster 4 markers
# #PAG4, PAG2, PAG5
# pdf("plots/regressUMI_CC/VlnPlot_C4_markers_regressUMICellCycle.pdf")
# VlnPlot(tryp_obj, features = c("Tb927.10.10210", "Tb927.10.10220", "Tb927.10.10230"))
# dev.off()

pdf("plots/regressUMI_CC/heatmap_withoutRibosomalProteins_regressUMICellCycle.pdf")
DoHeatmap(tryp_obj, features = top10$gene)
dev.off()

ribosomal_genes <- all927$Gene.ID[str_detect(all927$Product.Description, "ribosomal protein|ribosomal subunit protein")]
ribosomal_genes <- intersect(ribosomal_genes, row.names(tryp_obj))


#Find ribosomal proteins that are expressed in more than 10% of cells

exp_mtx <- tryp_obj@assays$RNA@layers$counts
row.names(exp_mtx) <- row.names(tryp_obj)

riboProtein_percent <- rowSums(exp_mtx[ ribosomal_genes,] > 0) / dim(tryp_obj)[2]

test <- riboProtein_percent

names(test) <- ribosomal_genes

riboProtein_percent <- riboProtein_percent[which(riboProtein_percent >= 0.1)]

tryp_obj <- MetaFeature(tryp_obj, features = names(riboProtein_percent), meta.name = "RiboProtein")

tryp_obj$RiboProtein <- tryp_obj$RiboProtein / mean(tryp_obj$RiboProtein)


pdf("plots/regressUMI_CC/UMAP_RiboProteinMetaFeature_regress_UMICellCycle.pdf")
FeaturePlot(tryp_obj, features = "RiboProtein")
dev.off()

pdf("plots/regressUMI_CC/PCA_RiboProteinMetaFeature_regress_UMICellCycle.pdf")
FeaturePlot(tryp_obj, features = "RiboProtein", reduction = "pca")
dev.off()

#Relationship between cell PC value and ribosomal protein metafeature
for (i in 1:4){
  plot_df <- data.frame("Ribosomal_protein_feature" = tryp_obj$RiboProtein,
                        "PC" = tryp_obj@reductions$pca@cell.embeddings[,i])
  
  pdf(paste0("plots/regressUMI_CC/RB_meta_vs_PC",i,"_regress_UMICellCycle.pdf"))
  print(ggplot(data = plot_df, mapping = aes(x = Ribosomal_protein_feature,
                                       y = PC)) + geom_point()+
    sm_statCorr(alternative = "less") + ylab(paste0("PC",i)) + xlab("Ribosomal protein feature value"))
  dev.off()
  
  print(ggplot(data = plot_df, mapping = aes(x = Ribosomal_protein_feature,
                                             y = PC)) + geom_point()+
          sm_statCorr(alternative = "less") + ylab(paste0("PC",i)) + xlab("Ribosomal protein feature value"))
  
}

pdf("plots/regressUMI_CC/VizDimLoading_PC1_regress_UMICellCycle.pdf")
VizDimLoadings(tryp_obj, dims = 1:1, nfeatures = 50, balanced = T)
dev.off()

DimPlot(tryp_obj, group.by = "seurat_clusters", reduction = "pca")

saveRDS(tryp_obj, "RBP6_analysis/objects/RBP6_regressUMICellCycle.obj")


