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

#This is not mapped against VSGS

setwd("/data/Ross/lucy_rbp6_analysis/paper")

Cell_cycle_regulated_genes <- read.delim("Cell_cycle_regulated_genes.text")

cellcycle_all <- unique(c(Cell_cycle_regulated_genes$Early.G1, Cell_cycle_regulated_genes$Late.G1, Cell_cycle_regulated_genes$S.phase,
                          Cell_cycle_regulated_genes$G2.M.phase))

tryp_obj <- readRDS("RBP6_analysis/objects/2023_10X_sample_phase_QC_mVSG.rds")

all927 <- read.csv("all927_genes.csv")

#remove genes ID which occur multiple times

single_genes <- names(table(all927$Gene.ID))[which(table(all927$Gene.ID) == 1)]

all927 <- subset(all927, Gene.ID %in% single_genes)

row.names(all927) <- all927$Gene.ID


ribosomal_genes <- all927$Gene.ID[str_detect(all927$Product.Description, "ribosomal protein|ribosomal subunit protein")]
ribosomal_genes <- intersect(ribosomal_genes, row.names(tryp_obj))

#Find ribosomal proteins that are expressed in more than 10% of cells

exp_mtx <- tryp_obj@assays$RNA@layers$counts
row.names(exp_mtx) <- row.names(tryp_obj)

riboProtein_percent <- rowSums(exp_mtx[ ribosomal_genes,] > 0) / dim(tryp_obj)[2]

riboProtein_percent <- riboProtein_percent[which(riboProtein_percent >= 0.1)]

tryp_obj <- MetaFeature(tryp_obj, features = names(riboProtein_percent), meta.name = "RiboProtein")

tryp_obj$RiboProtein <- tryp_obj$RiboProtein / mean(tryp_obj$RiboProtein)

var_features <- tryp_obj@assays[["RNA"]]@meta.data[["var.features"]][is.na(tryp_obj@assays[["RNA"]]@meta.data[["var.features"]]) == F]

var_features <- setdiff(var_features, cellcycle_all)
var_features <- setdiff(var_features, ribosomal_genes)


tryp_obj <- ScaleData(tryp_obj, vars.to.regress = c("nCount_RNA", "G2M.aggregate", "G1e.aggregate",
                                                    "S.aggregate", "G1l.aggregate", "RiboProtein"), feature = var_features)

tryp_obj <- RunPCA(tryp_obj, npcs = 50, feature = var_features)

ElbowPlot(tryp_obj, ndims = 50)

pdf("plots/regressUMI_CC_ribo/elbowPlot_regressRiboUMICellCycle.pdf")
ElbowPlot(tryp_obj, ndims = 50)
dev.off()

tryp_obj <- FindNeighbors(tryp_obj, dims = 1:15)
tryp_obj <- FindClusters(tryp_obj, resolution = c(0.1, 0.2, 0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

pdf("plots/regressUMI_CC_ribo/clustree_regressRiboUMICellCycle.pdf", height = 14)
clustree(tryp_obj, prefix = "RNA_snn_res.")
dev.off()

tryp_obj <- FindClusters(tryp_obj, resolution =  0.5)

tryp_obj <- RunUMAP(tryp_obj, dims = 1:15)

DimPlot(tryp_obj, label = T, group.by = "seurat_clusters")


pdf("plots/regressUMI_CC_ribo/UMAP_regressRiboUMICellCycle_plotClusters.pdf")
DimPlot(tryp_obj, label = T, group.by = "seurat_clusters")
dev.off()

pdf("plots/regressUMI_CC_ribo/UMAP_regressRiboUMICellCycle_plotCellCyclePhase.pdf")
DimPlot(tryp_obj, group.by = "Phase")
dev.off()

pdf("plots/regressUMI_CC_ribo/PCA_regressRiboUMICellCycle_plotClusters.pdf")
DimPlot(tryp_obj, label = T, group.by = "seurat_clusters", reduction = "pca")
dev.off()

pdf("plots/regressUMI_CC_ribo/PCA_regressRiboUMICellCycle_plotCellCyclePhase.pdf")
DimPlot(tryp_obj, label = T, group.by = "Phase", reduction = "pca")
dev.off()

pdf("plots/regressUMI_CC_ribo/UMAP_RiboProtein_featurePlot_regressRiboUMICellCycle.pdf")
FeaturePlot(tryp_obj, features = "RiboProtein")
dev.off()

pdf("plots/regressUMI_CC_ribo/PCA_RiboProtein_featurePlot_regressRiboUMICellCycle.pdf")
FeaturePlot(tryp_obj, features = "RiboProtein", reduction = "pca")
dev.off()


tryp_obj <- RenameIdents(tryp_obj, "0" = "ZC3H36_EP_hi_GPEET_mid", "1" = "RHS_EP_hi_GPEET_mid", "2" = "histone_EP_hi_GPEET_mid", "3" = "EP_hi_GPEET_mid", 
                         "4" = "late_RBP6_PAG", "5" = "early_Procyclics", "6" = "late_possibleMetacyclics")

tryp_obj$cellType <- Idents(tryp_obj)

pdf("plots/regressUMI_CC_ribo/Vln_gpeet_ep_barp_regressRiboUMICellCycle.pdf", width = 30, height = 14)
VlnPlot(tryp_obj, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"), group.by = "cellType")
dev.off()

DimPlot(tryp_obj, label = T, group.by = "seurat_clusters")

plot_df <- data.frame("Ribosomal_protein_feature" = tryp_obj$RiboProtein,
                      "PC" = tryp_obj@reductions$pca@cell.embeddings[,1])

pdf(paste0("plots/regressUMI_CC_ribo/RB_meta_vs_PC1_regress_UMICellCycle.pdf"))
print(ggplot(data = plot_df, mapping = aes(x = Ribosomal_protein_feature,
                                           y = PC)) + geom_point()+
        sm_statCorr(alternative = "less") + ylab(paste0("PC","1")) + xlab("Ribosomal protein feature value"))
dev.off()


markers <- FindAllMarkers(tryp_obj, only.pos = T, logfc.threshold = 0.5, min.pct = 0.05)

markers_subset <- subset(markers, p_val_adj < 0.05)

markers_subset %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  ungroup() -> top15_rbp6

markers_subset %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> top10

avgexp = AggregateExpression(tryp_obj, return.seurat = T, group.by = 'cellType')

pdf("plots/regressUMI_CC_ribo/MarkerHeatmap.pdf", width = 7, height = 10)
DoHeatmap(avgexp, features = top10$gene)
dev.off()

write.csv(markers_subset, "RBP6_analysis/DEG_tables/clusterMarkers_riboCCUMIregress.csv")

write.csv(top10, "RBP6_analysis/DEG_tables/top10_clusterMarkers_riboCCUMIregress.csv")



pdf("plots/regressUMI_CC_ribo/UMAP_regressUMIRiboCC_plotCellType.pdf", width = 10)
DimPlot(tryp_obj, label = T, reduction = "umap")
dev.off()

markers <- FindAllMarkers(tryp_obj, only.pos = T, logfc.threshold = 0.25, min.pct = 0.05)

markers_subset <- subset(markers, p_val_adj < 0.05)

write.csv(markers_subset, "RBP6_analysis/DEG_tables/clusterMarkers_riboCCUMIregress.csv")

saveRDS(tryp_obj, "RBP6_analysis/objects/RBP6_regressRiboUMICellCycle.obj")

library(RColorBrewer)
pt <- table(tryp_obj$Phase,Idents(tryp_obj))
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/regressUMI_CC_ribo/barplot_cellType_phase.pdf")
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster ID") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()

#Cluster 4
FeaturePlot(tryp_obj, features = "Tb927.7.6600")

#Cluster 1 (ZC3H36)
FeaturePlot(tryp_obj, features = "Tb927.10.12760")

#Compare PAGs and metacyclics
pag_vs_meta <- FindMarkers(tryp_obj, ident.1 = "late_Metacyclics", ident.2 = "late_RBP6_PAG", logfc.threshold = 0.5)
pag_vs_meta$pct_diff <- abs(pag_vs_meta$pct.1 - pag_vs_meta$pct.2)
pag_vs_meta <- subset(pag_vs_meta, p_val_adj < 0.05)

#RBP6 point mutation overexpression associated genes (compared with WT RBP6 overexpression)
#ZC3H36, RBP5, 
#MCP12, STE/STE11 serine/threonine-protein kinase, 
#GRESAG 4, PSP1 C-terminal conserved region
#Hypothetical, hypothetical
VlnPlot(tryp_obj, features = c("Tb927.10.12760", "Tb927.11.12100",
                               "Tb927.10.12840", "Tb927.6.2030",
                               "Tb927.6.430", "Tb927.9.9370",
                               "Tb927.4.4940", "Tb927.10.4340"))

cluster4_products <- all927[subset(top15_rbp6, cluster == 4)$gene,]

FeaturePlot(tryp_obj, features = c("Tb09.v4.0119", "Tb927.11.17870"))

#Clusters markers
#RHS for cluster 0, ZC3H36 for Cluster 1, PAG4 for cluster 2, GPEET for cluster 3, SGM1.7 for cluster 4
FeaturePlot(tryp_obj, features = c("Tb927.2.240",
                                   "Tb927.10.12760",
                                   "Tb927.10.10210",
                                   "Tb927.6.510",
                                   "Tb927.7.6600"))

VlnPlot(tryp_obj, features = c("Tb927.2.240",
                                   "Tb927.10.12760",
                                   "Tb927.10.10210",
                                   "Tb927.6.510",
                                   "Tb927.7.6600"))

#ZC3H36, ZC3H37, ZC3H38
VlnPlot(tryp_obj, features = c("Tb927.10.12760", "Tb927.10.12780", "Tb927.10.12800"))

#GPEET, EP1, EP2 and BARP 
pdf("plots/regressUMI_CC_ribo/VlnPlot_GPEET_EPs_BARP_regressRiboUMICellCycle.pdf")
VlnPlot(tryp_obj, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"))
dev.off()

pdf("plots/regressUMI_CC_ribo/FeaturePlot_GPEET_EPs_BARP_regressRiboUMICellCycle.pdf")
FeaturePlot(tryp_obj, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"))
dev.off()




#Plot marker genes
#EP1, BARP, RBP6
VlnPlot(tryp_obj, features = c("Tb927.10.10260", "Tb927.9.15640", "Tb927.3.2930"))
FeaturePlot(tryp_obj, features = c("Tb927.10.10260", "Tb927.9.15640", "Tb927.3.2930"))



