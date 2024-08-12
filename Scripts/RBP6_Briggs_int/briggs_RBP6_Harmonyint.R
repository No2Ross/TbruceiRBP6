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

setwd("/data/Ross/lucy_rbp6_analysis/paper")

tryp_obj <- readRDS("RBP6_analysis/objects/RBP6_regressRiboUMICellCycle.obj")
pcf_obj <- readRDS("/data/Ross/lucy_rbp6_analysis/paper/integration_analysis/objects/briggs_processed.rds")

pcf_obj$orig.ident <- "WT"
tryp_obj$orig.ident <- "RBP6"

tryp_obj$condition <- "in_vitro"
pcf_obj$condition <- "in_vivo"

tryp_obj$origin <- "Laidlaw"
pcf_obj$origin <- "Briggs"

tryp_obj$cellType <- tryp_obj@active.ident
pcf_obj$cellType <- pcf_obj@active.ident

shared_genes <- intersect(row.names(tryp_obj), row.names(pcf_obj))


merge_obj <- merge(tryp_obj, pcf_obj)

merge_obj <- merge_obj[shared_genes,]

norm_scale_factor <- median(c(tryp_obj$nCount_RNA, pcf_obj$nCount_RNA))


merge_obj <- NormalizeData(merge_obj, scale.factor = norm_scale_factor)


merge_obj <- FindVariableFeatures(merge_obj)



merge_obj <- ScaleData(merge_obj, vars.to.regress = c("nCount_RNA", "G2M.aggregate", "G1e.aggregate",
                                                      "S.aggregate", "G1l.aggregate", "RiboProtein"))
merge_obj <- RunPCA(merge_obj)

ElbowPlot(merge_obj)

merge_obj <- IntegrateLayers(object = merge_obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",
                             verbose = FALSE)

# re-join layers after integration
merge_obj[["RNA"]] <- JoinLayers(merge_obj[["RNA"]])

merge_obj <- FindNeighbors(merge_obj, reduction = "harmony", dims = 1:11)

merge_obj <- FindClusters(merge_obj, resolution = c(0.1, 0.2, 0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

pdf("plots/briggsRBP6_harmony/clustree_BriggsRBP6_Harmonyint.pdf", height = 14)
clustree(merge_obj, prefix = "RNA_snn_res.")
dev.off()

merge_obj <- FindClusters(merge_obj, resolution = 0.4)

merge_obj <- RunUMAP(merge_obj, dims = 1:11, reduction = "harmony")

DimPlot(merge_obj, group.by = "orig.ident", reduction = "pca")

DimPlot(merge_obj, group.by = "orig.ident")

DimPlot(merge_obj, group.by = "orig.ident", reduction = "harmony")

DimPlot(merge_obj, group.by = "cellType", label = T)

DimPlot(merge_obj, group.by = "seurat_clusters", label = T)

#Order the cell type names
merge_obj$cellType <- factor(merge_obj$cellType, levels = c("Fresh",
                                                            "early_Procyclics", "EP_hi_GPEET_mid", "histone_EP_hi_GPEET_mid", "RHS_EP_hi_GPEET_mid",
                                                            "ZC3H36_EP_hi_GPEET_mid","late_RBP6_PAG", "late_Metacyclics"))

pdf("plots/briggsRBP6_harmony/UMAP_cluster_splitOrigin_BriggsRBP6_Harmonyint.pdf", width = 14)
DimPlot(merge_obj, group.by = "seurat_clusters", split.by = "origin", label = T)
dev.off()

pdf("plots/briggsRBP6_harmony/UMAP_cellType_splitOrigin_BriggsRBP6_Harmonyint.pdf", width = 14)
DimPlot(merge_obj, group.by = "cellType", split.by = "origin")
dev.off()

pdf("plots/briggsRBP6_harmony/VlnPlot_GPEET_EPs_BARP_RBP6_pcf_Harmonyint.pdf", width = 20, height = 14)
VlnPlot(merge_obj, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"), group.by = "cellType", ncol = 2)
dev.off()

pdf("plots/briggsRBP6_harmony/FeaturePlot_GPEET_EPs_BARP_RBP6_pcf_Harmonyint.pdf")
FeaturePlot(merge_obj, features = c("Tb927.6.510","Tb927.10.10260", "Tb927.10.10250", "Tb927.9.15640"))
dev.off()

pdf("plots/briggsRBP6_harmony/UMAP_CC_pcf_Harmonyint.pdf")
DimPlot(merge_obj, group.by  = c("Phase"))
dev.off()

pdf("plots/briggsRBP6_harmony/UMAP_Ribo_pcf_Harmonyint.pdf")
FeaturePlot(merge_obj, features = c("RiboProtein"))
dev.off()

#Cluster markers
Idents(merge_obj) <- "cellType"

markers <- FindAllMarkers(merge_obj, only.pos = T, logfc.threshold = 0.5, min.pct = 0.05)

markers_subset <- subset(markers, p_val_adj < 0.05)

markers_subset %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  ungroup() -> top15_rbp6

#markers plot
#ZC3H36, #Histone H4
#RHS5, VSG interrupted (pseudogene)
#PAG4, Hexokinase
#GPEET
VlnPlot(merge_obj, features = c("Tb927.10.12760", "Tb927.5.4180",
                                "Tb927.2.240", "Tb927.11.17870",
                                "Tb927.10.10210", "Tb927.10.2010",
                                "Tb927.6.510"))


#early_Procyclics
VlnPlot(merge_obj, features = c("Tb927.6.2790"))

Idents(merge_obj) <- "seurat_clusters"

library(RColorBrewer)
pt <- table(Idents(merge_obj), merge_obj$origin)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf("plots/briggsRBP6_harmony/PCF_RBP6_Harmonyint_proportionPlot_sample_cluster.pdf")
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

pdf("plots/briggsRBP6_harmony/PCF_RBP6_Harmonyint_proportionPlot_cluster_sample.pdf", width = 14)
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

pdf("plots/briggsRBP6_harmony/PCF_RBP6_Harmonyint_proportionPlot_cluster_cellType.pdf", width = 14)
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

pdf("plots/briggsRBP6_harmony/PCF_RBP6_Harmonyint_proportionPlot_cellType_cluster.pdf", width = 20)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()

deg_path <- "/data/Ross/lucy_rbp6_analysis/paper/DEG/RBP6Briggs/"


Idents(merge_obj) <- "cellType"
stats_all <- FindMarkers(merge_obj, ident.1 = "early_Procyclics", ident.2 = "Fresh", logfc.threshold = 0, min.pct = 0)


markers <- FindMarkers(merge_obj, ident.1 = "early_Procyclics", ident.2 = "Fresh", logfc.threshold = 0.5, min.pct = 0.1)
markers <- subset(markers, pct.2 > 0 & pct.1 > 0 & p_val_adj < 0.05)
write.csv(markers, paste0(deg_path, "PCF_briggsVSrbp6_markers.csv"))

rbp6_markers <- subset(markers, avg_log2FC > 0)
briggs_markers <- subset(markers, avg_log2FC < 0)

VlnPlot(merge_obj, features = "Tb927.7.5940")

rbp6_markers %>%
  slice_max(n = 25, order_by = avg_log2FC) %>%
  ungroup() -> top25_rbp6

briggs_markers %>%
  slice_max(n = 25, order_by = avg_log2FC) %>%
  ungroup() -> top25_briggs

briggs_markers$pct_diff <- abs(briggs_markers$pct.1 - briggs_markers$pct.2)
rbp6_markers$pct_diff <- abs(rbp6_markers$pct.1 - rbp6_markers$pct.2)

write.csv(rbp6_markers, "/data/Ross/lucy_rbp6_analysis/paper/integration_analysis/DEG_tables/RBP6_vs_briggs_markers.csv")
write.csv(briggs_markers, "/data/Ross/lucy_rbp6_analysis/paper/integration_analysis/DEG_tables/briggs_vs_RBP6_markers.csv")

avgexp = AggregateExpression(merge_obj, return.seurat = T, group.by = 'origin')

Idents(merge_obj) <- "origin"
bulk_comparison <- FindMarkers(merge_obj, ident.1 = "Laidlaw", ident.2 = "Briggs", features = row.names(markers), 
                               logfc.threshold = 0, min.pct = 0)

bulk_comparison <- bulk_comparison[order(match(row.names(bulk_comparison), row.names(markers))),]
write.csv(bulk_comparison, paste0(deg_path, "bulkOrigin_markers.csv"))

markers$globalpct.1 <- bulk_comparison$pct.1
markers$globalpct.2 <- bulk_comparison$pct.2

markers_laidlaw <- subset(markers, avg_log2FC > 0)

markers_briggs <- subset(markers, avg_log2FC < 0)


#See what the mean pct is for the two conditions for the genes
mean(markers_briggs$globalpct.1)
mean(markers_briggs$globalpct.2)
mean(markers_briggs$pct.1)
mean(markers_briggs$pct.2)

mean(markers_laidlaw$globalpct.1)
mean(markers_laidlaw$globalpct.2)
mean(markers_laidlaw$pct.1)
mean(markers_laidlaw$pct.2)


