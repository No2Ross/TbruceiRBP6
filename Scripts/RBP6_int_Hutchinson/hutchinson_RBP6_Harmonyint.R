#### Load libraries ####
library(stringr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(dplyr)
library(harmony)
library(phateR)
library(clustree)

setwd("/data/Ross/lucy_rbp6_analysis/paper")

tryp_obj <- readRDS("RBP6_analysis/objects/RBP6_regressRiboUMICellCycle.obj")
hutchinson_obj <- readRDS("/data/Ross/lucy_rbp6_analysis/paper/integration_analysis/objects/seb_lib1_processed.rds")

hutchinson_obj$orig.ident <- "WT"
tryp_obj$orig.ident <- "RBP6"

tryp_obj$condition <- "in_vitro"
hutchinson_obj$condition <- "in_vivo"

tryp_obj$origin <- "Laidlaw"
hutchinson_obj$origin <- "Hutchinson"

tryp_obj$cellType <- tryp_obj@active.ident
hutchinson_obj$cellType <- hutchinson_obj@active.ident

shared_genes <- intersect(row.names(tryp_obj), row.names(hutchinson_obj))

merge_obj <- merge(tryp_obj, hutchinson_obj)

merge_obj <- merge_obj[shared_genes,]

norm_scale_factor <- median(c(tryp_obj$nCount_RNA, hutchinson_obj$nCount_RNA))

merge_obj <- NormalizeData(merge_obj, scale.factor = norm_scale_factor)


merge_obj <- FindVariableFeatures(merge_obj)


merge_obj <- ScaleData(merge_obj, vars.to.regress = c("nCount_RNA", "G2M.aggregate", "G1e.aggregate",
                                                      "S.aggregate", "G1l.aggregate", "RiboProtein"))
merge_obj <- RunPCA(merge_obj)

pdf("plots/hutchinsonRBP6_harmony/elbowPlot_HutchRBP6_Harmonyint.pdf")
ElbowPlot(merge_obj)
dev.off()

merge_obj <- IntegrateLayers(object = merge_obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",
                             verbose = FALSE)

# re-join layers after integration
merge_obj[["RNA"]] <- JoinLayers(merge_obj[["RNA"]])

merge_obj <- FindNeighbors(merge_obj, reduction = "harmony", dims = 1:12)

merge_obj <- FindClusters(merge_obj, resolution = c(0.1, 0.2, 0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

pdf("plots/hutchinsonRBP6_harmony/clustree_HutchRBP6_Harmonyint.pdf", height = 14)
clustree(merge_obj, prefix = "RNA_snn_res.")
dev.off()

merge_obj <- FindClusters(merge_obj, resolution = 0.3)

merge_obj <- RunUMAP(merge_obj, dims = 1:12, reduction = "harmony")

DimPlot(merge_obj)

DimPlot(merge_obj, group.by = "orig.ident", reduction = "pca")

DimPlot(merge_obj, group.by = "orig.ident")

DimPlot(merge_obj, group.by = "cellType", label = T)

#Order the cell type names
merge_obj$cellType <- factor(merge_obj$cellType, levels = c("Midgut forms", "Attached Epimastigote", "Gametes", "Pre-metacyclic", "Metacyclics",
                                                   "early_Procyclics", "EP_hi_GPEET_mid", "histone_EP_hi_GPEET_mid", "RHS_EP_hi_GPEET_mid",
                                                   "ZC3H36_EP_hi_GPEET_mid","late_RBP6_PAG", "late_possibleMetacyclics"))

DimPlot(merge_obj, group.by = "cellType", label = T)

pdf("plots/hutchinsonRBP6_harmony/UMAP_cluster_splitOrigin_HutchRBP6_Harmonyint.pdf", width = 14)
DimPlot(merge_obj, group.by = "seurat_clusters", split.by = "origin", label = T)
dev.off()

pdf("plots/hutchinsonRBP6_harmony/UMAP_cellType_splitOrigin_HutchRBP6_Harmonyint.pdf", width = 14)
DimPlot(merge_obj, group.by = "cellType", split.by = "origin")
dev.off()

pdf("plots/hutchinsonRBP6_harmony/UMAP_origin_HutchRBP6_Harmonyint.pdf", width = 7)
DimPlot(merge_obj,group.by = "origin")
dev.off()

DimPlot(merge_obj, group.by = "seurat_clusters", label = T)

pdf("plots/hutchinsonRBP6_harmony/FeaturePlot_GPEET_EPs_BARP_RBP6_Hutchinson_Harmonyint.pdf")
FeaturePlot(merge_obj, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"))
dev.off()


pdf("plots/hutchinsonRBP6_harmony/VlnPlot_GPEET_EPs_BARP_RBP6_Hutchinson_Harmonyint.pdf", width = 30, height = 14)
VlnPlot(merge_obj, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"), group.by = "cellType", ncol = 2)
dev.off()

VlnPlot(merge_obj, features = "Tb927.10.12760", group.by = "cellType")


library(RColorBrewer)
pt <- table(Idents(merge_obj), merge_obj$origin)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/hutchinsonRBP6_harmony/Harmonyint_proportionPlot_sample_cluster.pdf")
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()


pt <- table(merge_obj$origin, Idents(merge_obj))
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/hutchinsonRBP6_harmony/Harmonyint_proportionPlot_cluster_sample.pdf", width = 14)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()


pt <- table(merge_obj$cellType, Idents(merge_obj))
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/hutchinsonRBP6_harmony/Harmonyint_proportionPlot_cluster_cellType.pdf", width = 14)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()


pt <- table(Idents(merge_obj),merge_obj$cellType)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/hutchinsonRBP6_harmony/Harmonyint_proportionPlot_cellType_cluster.pdf", width = 30)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()


avgexp = AggregateExpression(merge_obj, return.seurat = T, group.by = 'cellType')

pdf("plots/hutchinsonRBP6_harmony/GPEET_EP_BARP_Heatmap.pdf", width = 7, height = 10)
DoHeatmap(avgexp, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"))
dev.off()



#Compare RBP6 vs Hutchinson

Idents(merge_obj) <- "origin"
rbp6_markers <- FindMarkers(merge_obj, ident.1 = "Laidlaw", ident.2 = "Hutchinson", logfc.threshold = 0.5, min.pct = 0.1, only.pos = T)
invivo_markers <- FindMarkers(merge_obj, ident.1 = "Hutchinson", ident.2 = "Laidlaw", logfc.threshold = 0.5, min.pct = 0.1, only.pos = T)

rbp6_markers <- subset(rbp6_markers, pct.2 > 0)
invivo_markers <- subset(invivo_markers, pct.2 > 0)

rbp6_markers %>%
  slice_max(n = 25, order_by = avg_log2FC) %>%
  ungroup() -> top25_rbp6

invivo_markers %>%
  slice_max(n = 25, order_by = avg_log2FC) %>%
  ungroup() -> top25_invivo

avgexp = AggregateExpression(merge_obj, return.seurat = T, group.by = 'origin')

pdf("plots/hutchinsonRBP6_harmony/Hutchinson_vs_laidlaw_heatmap_top25Each.pdf", width = 7, height = 10)
DoHeatmap(ScaleData(merge_obj, features = row.names(merge_obj)), features = unique( c( row.names(top25_rbp6), row.names(top25_invivo) ) ) )
dev.off()

write.csv(rbp6_markers, "/data/Ross/lucy_rbp6_analysis/paper/integration_analysis/DEG_tables/RBP6_vs_inVivo_markers.csv")
write.csv(invivo_markers, "/data/Ross/lucy_rbp6_analysis/paper/integration_analysis/DEG_tables/InVivo_vs_RBP6_markers.csv")

#PAG VlnPlot - PAG1, PAG2, PAG4, PAG5

Idents(merge_obj) <- "cellType"
VlnPlot(merge_obj, features = c("Tb927.10.10240", "Tb927.10.10220","Tb927.10.10210", "Tb927.10.10230"), ncol = 2)


