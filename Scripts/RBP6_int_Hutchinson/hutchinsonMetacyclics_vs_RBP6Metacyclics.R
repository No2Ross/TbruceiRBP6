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
library(batchelor)
library(SeuratWrappers)
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

merge_obj <- merge(tryp_obj, hutchinson_obj)

norm_scale_factor <- median(c(tryp_obj$nCount_RNA, hutchinson_obj$nCount_RNA))

#Check PAG expression
pdf("plots/Hutchinson_PAG1_2_4_5_Vlnplot.pdf")
VlnPlot(hutchinson_obj, features = c("Tb927.10.10240", "Tb927.10.10220", "Tb927.10.10210", "Tb927.10.10230") )
dev.off()

merge_obj <- NormalizeData(merge_obj, scale.factor = norm_scale_factor)

merge_obj[["RNA"]] <- JoinLayers(merge_obj[["RNA"]])

barp_genes_table <- read.csv("/data/Ross/lucy_rbp6_analysis/paper/BARP_genes.csv")

#Notes
#Tb927.5.4020 unlikely to be BARP. See: A receptor for the complement regulator factor H increases transmission of trypanosomes to tsetse flies
barp_genes_table <- barp_genes_table[which(barp_genes_table$Product.Description == "BARP protein"),]

barp_genes <- barp_genes_table[which(barp_genes_table$Gene.ID %in% row.names(hutchinson_obj)), "Gene.ID"]

VlnPlot(hutchinson_obj, features = barp_genes[1:10] )
VlnPlot(tryp_obj, features = barp_genes[11:20] )
VlnPlot(tryp_obj, features = barp_genes[21:24] )

#BARP protein expression
VlnPlot(merge_obj, features = c("Tb927.9.15510", "Tb927.9.15520", "Tb927.9.15530",
                                "Tb927.9.15540", "Tb927.9.15550", "Tb927.9.15560",
                                "Tb927.9.15570", "Tb927.9.15580", "Tb927.9.15590",
                                "Tb927.9.15600", "Tb927.9.15610", "Tb927.9.15620",
                                "Tb927.9.15630", "Tb927.9.15640"))

#Order the cell type names
merge_obj$cellType <- factor(merge_obj$cellType, levels = c("Midgut forms", "Attached Epimastigote", "Gametes", "Pre-metacyclic", "Metacyclics",
                                                            "early_Procyclics", "EP_hi_GPEET_mid", "histone_EP_hi_GPEET_mid", "RHS_EP_hi_GPEET_mid",
                                                            "ZC3H36_EP_hi_GPEET_mid","late_RBP6_PAG", "late_possibleMetacyclics"))

deg_path <- "/data/Ross/lucy_rbp6_analysis/paper/DEG/RBP6Meta_vs_HutchinsonMeta/"
output_table <- list()

#Cluster 4 (metacyclic) differences
meta_DE <- FindMarkers(merge_obj, 
                        ident.1 = "late_possibleMetacyclics", ident.2 = "Metacyclics",
                        logfc.threshold = 0.5, min.pct = 0.1)
meta_DE <- subset(meta_DE, p_val_adj < 0.05 & pct.1 > 0 & pct.2 > 0)
write.csv(meta_DE, paste0(deg_path, "metacyclicsDE.csv"))


library(EnhancedVolcano)

subset(meta_DE, avg_log2FC > 0) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> meta_DE_laidlaw_10

subset(meta_DE, avg_log2FC < 0) %>%
  slice_min(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> meta_DE_hutchinson_10

pdf("plots/RBP6Meta_vs_HutchinsonMeta/metacyclic_volcanoPlot.pdf")
EnhancedVolcano(meta_DE,
                lab = rownames(meta_DE),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c( row.names(meta_DE_laidlaw_10), row.names(meta_DE_hutchinson_10) ))
dev.off()

Idents(merge_obj) <- "origin"
bulk_comparison <- FindMarkers(merge_obj, ident.1 = "Laidlaw", ident.2 = "Hutchinson", features = row.names(meta_DE), 
                               logfc.threshold = 0, min.pct = 0)

bulk_comparison <- bulk_comparison[order(match(row.names(bulk_comparison), row.names(meta_DE))),]
write.csv(bulk_comparison, paste0(deg_path, "bulkOrigin_meta_DE_metacyclics.csv"))

meta_DE$globalpct.1 <- bulk_comparison$pct.1
meta_DE$globalpct.2 <- bulk_comparison$pct.2

library(mdthemes)
pdf("plots/RBP6Meta_vs_HutchinsonMeta/metacyclic_pctGlobal.pdf")
ggplot() + geom_point(aes(meta_DE$globalpct.1 * 100, meta_DE$globalpct.2 * 100)) + 
  theme_minimal() +
  xlab("Percentage of all RBP6 overexpression dataset cells that express the DE gene") +
  ylab("Percentage of all *in vivo* insect stage dataset cells that express the DE gene") +
  theme(axis.title.y = ggtext::element_markdown())
dev.off()

meta_DE_laidlaw <- subset(meta_DE, avg_log2FC > 0)
less_than_10percent_laidlaw <- subset(meta_DE_laidlaw, pct.2 < 0.1)
less_than_10percent_laidlaw_global <- subset(meta_DE_laidlaw, globalpct.2 < 0.1)

meta_DE_hutchinson <- subset(meta_DE, avg_log2FC < 0)
less_than_10percent_hutchinson <- subset(meta_DE_hutchinson, pct.1 < 0.1)
less_than_10percent_hutchinson_global <- subset(meta_DE_hutchinson, globalpct.1 < 0.1)

#See what the mean pct is for the two conditions for the genes
mean(meta_DE_hutchinson$globalpct.1)
mean(meta_DE_hutchinson$globalpct.2)

mean(meta_DE_laidlaw$globalpct.1)
mean(meta_DE_laidlaw$globalpct.2)
median(meta_DE_laidlaw$globalpct.1)
median(meta_DE_laidlaw$globalpct.2)

Idents(merge_obj) <- "cellType"

#Expression of Christiano 2017 genes
#ISG65, ISG64, Alternative oxidase
VlnPlot(merge_obj, features = c("Tb927.2.3270", "Tb927.5.1410", "Tb927.10.7090"), group.by = "cellType")

#PAG expression
#PAG1, PAG2, PAG4, PAG5
pdf("plots/RBP6Meta_vs_HutchinsonMeta/PAG1_2_4_5_Vlnplot.pdf")
VlnPlot(merge_obj, features = c("Tb927.10.10240", "Tb927.10.10220", "Tb927.10.10210", "Tb927.10.10230"))
dev.off()




