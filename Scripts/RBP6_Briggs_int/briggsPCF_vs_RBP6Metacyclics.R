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

merge_obj <- merge(tryp_obj, pcf_obj)

norm_scale_factor <- median(c(tryp_obj$nCount_RNA, pcf_obj$nCount_RNA))

merge_obj <- NormalizeData(merge_obj, scale.factor = norm_scale_factor)

merge_obj[["RNA"]] <- JoinLayers(merge_obj[["RNA"]])

#Order the cell type names
merge_obj$cellType <- factor(merge_obj$cellType, levels = c("Fresh",
                                                            "early_Procyclics", "EP_hi_GPEET_mid", "histone_EP_hi_GPEET_mid", "RHS_EP_hi_GPEET_mid",
                                                            "ZC3H36_EP_hi_GPEET_mid","late_RBP6_PAG", "late_possibleMetacyclics"))

deg_path <- "/data/Ross/lucy_rbp6_analysis/paper/DEG/RBP6Meta_vs_briggsPCF/"
output_table <- list()

#Cluster 4 (metacyclic) differences
meta_DE <- FindMarkers(merge_obj, 
                       ident.1 = "late_possibleMetacyclics", ident.2 = "Fresh",
                       logfc.threshold = 0.5, min.pct = 0.1)
meta_DE <- subset(meta_DE, p_val_adj < 0.05 & pct.1 > 0 & pct.2 > 0)
write.csv(meta_DE, paste0(deg_path, "metacyclics_FreshDE.csv"))


library(EnhancedVolcano)

subset(meta_DE, avg_log2FC > 0) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> meta_DE_laidlaw_10

subset(meta_DE, avg_log2FC < 0) %>%
  slice_min(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> meta_DE_briggs_10

pdf("plots/RBP6Meta_vs_briggsPCF/metacyclic_Fresh_volcanoPlot.pdf")
EnhancedVolcano(meta_DE,
                lab = rownames(meta_DE),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c( row.names(meta_DE_laidlaw_10), row.names(meta_DE_briggs_10) ))
dev.off()

Idents(merge_obj) <- "origin"
bulk_comparison <- FindMarkers(merge_obj, ident.1 = "Laidlaw", ident.2 = "Briggs", features = row.names(meta_DE), 
                               logfc.threshold = 0, min.pct = 0)

bulk_comparison <- bulk_comparison[order(match(row.names(bulk_comparison), row.names(meta_DE))),]
write.csv(bulk_comparison, paste0(deg_path, "bulkOrigin_meta_Fresh_DE.csv"))

meta_DE$globalpct.1 <- bulk_comparison$pct.1
meta_DE$globalpct.2 <- bulk_comparison$pct.2

library(mdthemes)
pdf("plots/RBP6Meta_vs_briggsPCF/metacyclic_Fresh_pctGlobal.pdf")
ggplot() + geom_point(aes(meta_DE$globalpct.1 * 100, meta_DE$globalpct.2 * 100)) + 
  theme_minimal() +
  xlab("Percentage of all RBP6 overexpression dataset cells that express the DE gene") +
  ylab("Percentage of all *in vivo* insect stage dataset cells that express the DE gene") +
  theme(axis.title.y = ggtext::element_markdown())
dev.off()

meta_DE_laidlaw <- subset(meta_DE, avg_log2FC > 0)
less_than_10percent_laidlaw <- subset(meta_DE_laidlaw, pct.2 < 0.1)
less_than_10percent_laidlaw_global <- subset(meta_DE_laidlaw, globalpct.2 < 0.1)

meta_DE_briggs <- subset(meta_DE, avg_log2FC < 0)
less_than_10percent_briggs <- subset(meta_DE_briggs, pct.1 < 0.1)
less_than_10percent_briggs_global <- subset(meta_DE_briggs, globalpct.1 < 0.1)

#See what the mean pct is for the two conditions for the genes
mean(meta_DE_briggs$globalpct.1)
mean(meta_DE_briggs$globalpct.2)

mean(meta_DE_laidlaw$globalpct.1)
mean(meta_DE_laidlaw$globalpct.2)
median(meta_DE_laidlaw$globalpct.1)
median(meta_DE_laidlaw$globalpct.2)

Idents(merge_obj) <- "cellType"

#Expression of Christiano 2017 genes
#ISG65, ISG64, Alternative oxidase
VlnPlot(merge_obj, features = c("Tb927.2.3270", "Tb927.5.1410", "Tb927.10.7090"), group.by = "cellType")

#Threonine catabolism - TDH, KBL
pdf("plots/RBP6Meta_vs_briggsPCF/metacyclic_Fresh_threonine.pdf")
VlnPlot(merge_obj, features = c("Tb927.6.2790", "Tb927.8.6060"))
dev.off()

#Proline catabolism - PDH , delta
pdf("plots/RBP6Meta_vs_briggsPCF/metacyclic_Fresh_proline.pdf")
VlnPlot(merge_obj, features = c("Tb927.7.210", "Tb927.10.3210"))
dev.off()

#VSG found signifciantly higher in RBP6 metacyclics
pdf("plots/RBP6Meta_vs_briggsPCF/metacyclic_Fresh_VSG.pdf")
VlnPlot(merge_obj, features = c("Tb927.11.17870", "Tb11.57.0047b", "Tb927.1.5060", "Tb927.2.2060",
                                "Tb927.9.7380", "Tb927.9.7390", "Tb11.v5.0213"))
dev.off()


