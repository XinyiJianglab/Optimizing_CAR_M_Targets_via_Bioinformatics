
library(Seurat)
library(dplyr)
library(qs2)

# downloadmac.atlas.zenodo.200524.rds from 
# https://zenodo.org/records/11222158

zenodo <- readRDS("mac.atlas.zenodo.200524.rds")

# match meta.data
identical(rownames(zenodo@meta.data),colnames(zenodo[["RNA"]]))
identical(rownames(zenodo@meta.data), colnames(zenodo[["SCT"]]))
identical(rownames(zenodo@meta.data), colnames(zenodo[["integrated"]]))

identical(zenodo@meta.data[["cellid"]], colnames(zenodo[["RNA"]]))
rownames(zenodo@meta.data) <- zenodo@meta.data$cellid

identical(rownames(zenodo@meta.data),colnames(zenodo[["RNA"]]))
identical(rownames(zenodo@meta.data), colnames(zenodo[["SCT"]]))
identical(rownames(zenodo@meta.data), colnames(zenodo[["integrated"]]))

# macrophage
zenodo_sub <- subset(zenodo,subset = 
                       (tissue %in% c("Tumor")) & 
                       !(short.label %in% c("20_TDoub","23_NA"))
                     )

# save by qs2

dir.create("input", recursive = TRUE, showWarnings = FALSE)

qs_save(zenodo_sub,"input/zenodo.200524.sub.qs", nthreads = 6)

# deg

DefaultAssay(zenodo_sub) <- "RNA"

Idents(zenodo_sub)<- "short.label"

macro_deg <- FindMarkers(
  zenodo_sub,
  ident.1 = c("8_IFNGMac", "17_IFNMac3", "22_IFNMac4", "10_InflamMac", "2_C3Mac"),
  ident.2 = c("1_MetM2Mac", "6_SPP1AREGMac", "18_ECMMac", "3_ICIMac1", "4_ICIMac2", "21_HemeMac", "9_AngioMac"),
  assay = "RNA",
  layer = "data",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0,
  only.pos = FALSE
) %>% filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.25)

write.csv(macro_deg, "input/macro_deg.csv", row.names = T)
