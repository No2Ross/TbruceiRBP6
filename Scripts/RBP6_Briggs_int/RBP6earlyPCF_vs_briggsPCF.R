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
#glycerol-3-phosphate dehydrogenase, alternative oxidase
#phosphofructokinase, phosphoglycerate kinase
#pyruvate kinase, phosphoglycerate mutase
lamour_metacyclic_genes <- c("Tb927.11.7380", "Tb927.10.7090",
                       "Tb927.3.3270", "Tb927.1.700",
                       "Tb927.10.14140", "Tb927.10.7930")

#Cytochrome C, cytochrome oxidase subunit IV
#complex V, ATP synthase alpha chain, complex III, rieske iron-sulfur protein
#Proline dehydrogenase, pyrroline-5 carboxylate dehydrogenase
#Glutamate dehydrogenase, alanine aminotransferase
#aconitase (ACO), citrate synthase (CS)
#mitochondrial malate dehydrogenase (mMDH), pyruvate dehydrogenase (PYDH)
#mitochondrial fumarate hydratase (mFH), mitochondrial fumarate reductase (mFRD)
#phosphoenolpyruvate carboxykinase (PEPCK), glycosomal fumarate hydratase (gFH)
#glycosomal malate dehydrogenase (gMDH), glycosomal fumarate hydratase (gFH)
#glycosomal fumarate reductase (gFRD), glyceraldehyde 3-phosphate dehydrogenase (G3PDH)
#triosephosphate isomerase (TPI), hexokinase (HK)
lamour_pcf_genes <- c("Tb927.8.5120", "Tb927.1.4100",
                      "Tb927.7.7420", "Tb927.9.14160",
                      "Tb927.7.210", "Tb927.10.3210",
                      "Tb927.9.5900", "Tb927.1.3950",
                      "Tb927.10.14000", "Tb927.10.13430",
                      "Tb927.10.2560", "Tb927.3.1790",
                      "Tb927.11.5050", "Tb927.5.940",
                      "Tb927.2.4210", "Tb927.3.4500",
                      "Tb927.10.15410", "Tb927.3.4500",
                      "Tb927.5.930", "Tb927.6.4280",
                      "Tb927.11.5520", "Tb927.10.2010")


deg_path <- "/data/Ross/lucy_rbp6_analysis/paper/DEG/RBP6earlyPCF_vs_briggsPCF/"
output_table <- list()

#Cluster 4 (pcf meta) differences
pcf_DE <- FindMarkers(merge_obj, 
                       ident.1 = "early_Procyclics", ident.2 = "Fresh",
                       logfc.threshold = 0.5, min.pct = 0.1)
pcf_DE <- subset(pcf_DE, p_val_adj < 0.05 & pct.1 > 0 & pct.2 > 0)
write.csv(pcf_DE, paste0(deg_path, "earlyPCF_FreshDE.csv"))

library(EnhancedVolcano)

subset(pcf_DE, avg_log2FC > 0) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> pcf_DE_laidlaw_10

subset(pcf_DE, avg_log2FC < 0) %>%
  slice_min(n = 10, order_by = avg_log2FC) %>%
  ungroup() -> pcf_DE_briggs_10

pdf("plots/RBP6earlyPCF_vs_briggsPCF/earlyPCF_Fresh_volcanoPlot.pdf")
EnhancedVolcano(pcf_DE,
                lab = rownames(pcf_DE),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c( row.names(pcf_DE_laidlaw_10), row.names(pcf_DE_briggs_10) ))
dev.off()

Idents(merge_obj) <- "origin"
bulk_comparison <- FindMarkers(merge_obj, ident.1 = "Laidlaw", ident.2 = "Briggs", features = row.names(pcf_DE), 
                               logfc.threshold = 0, min.pct = 0)

bulk_comparison <- bulk_comparison[order(match(row.names(bulk_comparison), row.names(pcf_DE))),]
write.csv(bulk_comparison, paste0(deg_path, "bulkOrigin_earlyPCF_Fresh_DE.csv"))

pcf_DE$globalpct.1 <- bulk_comparison$pct.1
pcf_DE$globalpct.2 <- bulk_comparison$pct.2

library(mdthemes)
pdf("plots/RBP6earlyPCF_vs_briggsPCF/earlyPCF_Fresh_pctGlobal.pdf")
ggplot() + geom_point(aes(pcf_DE$globalpct.1 * 100, pcf_DE$globalpct.2 * 100)) + 
  theme_minimal() +
  xlab("Percentage of all RBP6 overexpression dataset cells that express the DE gene") +
  ylab("Percentage of all *in vivo* insect stage dataset cells that express the DE gene") +
  theme(axis.title.y = ggtext::element_markdown())
dev.off()

pcf_DE_briggs <- subset(pcf_DE, avg_log2FC < 0)
pcf_DE_laidlaw <- subset(pcf_DE, avg_log2FC > 0)

intersect(lamour_pcf_genes, row.names(pcf_DE_briggs))
intersect(lamour_metacyclic_genes, row.names(pcf_DE_briggs))

intersect(lamour_pcf_genes, row.names(pcf_DE_laidlaw))
intersect(lamour_metacyclic_genes, row.names(pcf_DE_laidlaw))

less_than_10percent_laidlaw <- subset(pcf_DE_laidlaw, pct.2 < 0.1)
less_than_10percent_laidlaw_global <- subset(pcf_DE_laidlaw, globalpct.2 < 0.1)

less_than_10percent_briggs <- subset(pcf_DE_briggs, pct.1 < 0.1)
less_than_10percent_briggs_global <- subset(pcf_DE_briggs, globalpct.1 < 0.1)

#See what the mean pct is for the two conditions for the genes
mean(pcf_DE_briggs$globalpct.1)
mean(pcf_DE_briggs$globalpct.2)

mean(pcf_DE_laidlaw$globalpct.1)
mean(pcf_DE_laidlaw$globalpct.2)
median(pcf_DE_laidlaw$globalpct.1)
median(pcf_DE_laidlaw$globalpct.2)

pcf_briggs_mtx <- subset(merge_obj, cellType == "Fresh")[c("Tb927.8.1000", "Tb927.9.13570"),]@assays$RNA@layers$data
pcf_laidlaw_mtx <- subset(merge_obj, cellType == "early_Procyclics")[c("Tb927.8.1000", "Tb927.9.13570"),]@assays$RNA@layers$data

mean.fxn(pcf_laidlaw_mtx) - mean.fxn(pcf_briggs_mtx)
mean.fxn_old(pcf_laidlaw_mtx) - mean.fxn_old(pcf_briggs_mtx)

Idents(merge_obj) <- "cellType"

#Expression of Christiano 2017 genes
#ISG65, ISG64, Alternative oxidase
VlnPlot(merge_obj, features = c("Tb927.2.3270", "Tb927.5.1410", "Tb927.10.7090"), group.by = "cellType")

#Threonine catabolism - TDH, KBL
pdf("plots/RBP6earlyPCF_vs_briggsPCF/earlyPCF_Fresh_threonine.pdf")
VlnPlot(merge_obj, features = c("Tb927.6.2790", "Tb927.8.6060"))
dev.off()

#Proline catabolism - PDH , delta
pdf("plots/RBP6earlyPCF_vs_briggsPCF/earlyPCF_Fresh_proline.pdf")
VlnPlot(merge_obj, features = c("Tb927.7.210", "Tb927.10.3210"))
dev.off()

#VSG found signifciantly higher in RBP6 pcfcyclics
pdf("plots/RBP6earlyPCF_vs_briggsPCF/earlyPCF_Fresh_VSG.pdf")
VlnPlot(merge_obj, features = c("Tb927.11.17870", "Tb11.57.0047b", "Tb927.1.5060", "Tb927.2.2060",
                                "Tb927.9.7380", "Tb927.9.7390", "Tb11.v5.0213"))
dev.off()

