
library(Seurat)
library(harmony)
library(qs2)
library(data.table)
library(Matrix)
library(tidyverse)

# GSE232240 ------

# download
# GSE232240_Count_data_IMCISION.txt
# GSE232240_Meta_data_IMCISION
# from Supplementary file of GSE232240

count_file <- "GSE232240/GSE232240_Count_data_IMCISION.txt"
meta_file  <- "GSE232240/GSE232240_Meta_data_IMCISION.txt"

first_line <- readLines(count_file, n = 1)
cell_names <- strsplit(first_line,"\t", fixed = TRUE)[[1]]
count_df <- data.table::fread(count_file, sep = "\t", header = FALSE, skip = 1, 
                              fill = TRUE,data.table = FALSE, check.names = FALSE)
gene_names <- as.character(count_df[[1]])
count_df[[1]] <- NULL
stopifnot(length(cell_names) == ncol(count_df))
colnames(count_df) <- cell_names
count_matrix <- as.matrix(count_df)
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- make.unique(gene_names)
count_matrix <- Matrix::Matrix( count_matrix, sparse = TRUE)

obj1 <- CreateSeuratObject(counts = count_matrix,project = "GSE232240",min.cells = 0,min.features = 0)
meta <- data.table::fread(meta_file,sep = "\t",header = TRUE,data.table = FALSE,check.names = FALSE)
identical(meta$cell_id,colnames(obj1))
obj1 <- AddMetaData(object = obj1,metadata = meta)
obj1_sub <- subset(obj1,subset = mc_group == "Mono-macro")
obj1_sub@meta.data <- obj1_sub@meta.data %>%
  mutate(
    PATCS1 = case_when(
      timepoint == "pre" ~ "basic",
      timepoint == "post" & response == "RE" ~ "active",
      TRUE ~ NA_character_
    ),
    PATCS2 = as.character(patient),
    PATCS3 = "GSE232240"
  )
obj1_sub <- subset(obj1_sub, subset = !is.na(PATCS1))

PATCS_cell_count <- obj1_sub@meta.data %>%
  count(PATCS2, PATCS1, name = "n_cells") %>%
  arrange(PATCS1, PATCS2)
PATCS_cell_count
response_patient_list <- obj1_sub@meta.data %>%
  filter(response %in% c("NR", "RE")) %>%
  select(response, patient) %>%
  distinct() %>%
  arrange(response, patient)
response_patient_list
table(obj1_sub$PATCS1)

qs_save(obj1_sub,"GSE232240/GSE232240_sub.qs",nthreads = 4)

# GSE212217 ------

# download
# GSE212217_seurat_scRNAseq.v2.rds
# from Supplementary file of GSE212217

obj2 <- readRDS("GSE212217/GSE212217_seurat_scRNAseq.v2.rds")

obj2_sub <- subset(obj2,cells = rownames(obj2@meta.data)[
  grepl("TAM|mono", obj2@meta.data$finalIdent, ignore.case = TRUE)])
table(obj2_sub$finalIdent)
obj2_sub@meta.data <- obj2_sub@meta.data %>%
  dplyr::mutate(
    response_simple = dplyr::case_when(
      clinical %in% c("epiR", "mutR") ~ "RE",
      clinical == "NR"                ~ "NR",
      clinical == "Healthy"           ~ "Healthy",
      TRUE                            ~ NA_character_ ),
    PATCS1 = dplyr::case_when(
      timepointBinary == "pre" &
        clinical != "Healthy" ~ "basic",
      
      timepointBinary == "post" &
        clinical %in% c("epiR", "mutR") ~ "active",
      TRUE ~ NA_character_),
    PATCS2 = ifelse(
      grepl("^[0-9]+$", as.character(patient)),
      paste0("P", as.character(patient)),
      as.character(patient)),
    PATCS3 = "GSE212217"
    )
obj2_sub <- subset(obj2_sub,subset = !is.na(PATCS1))

PATCS_cell_count <- obj2_sub@meta.data %>%
  dplyr::count(PATCS3,PATCS2,PATCS1,response_simple,clinical,name = "n_cells" ) %>%
  dplyr::arrange(PATCS1, PATCS2)

response_patient_list <- obj2_sub@meta.data %>%
  dplyr::filter(response_simple %in% c("NR", "RE")) %>%
  dplyr::select(response = response_simple, clinical, patient = PATCS2) %>%
  dplyr::distinct() %>%
  dplyr::arrange(response, clinical, patient)

qs_save(obj2_sub,"GSE212217/GSE212217_sub.qs",nthreads = 4)

# GSE246613 ------

# download
# GSE246613_PembroRT_immune_R100_final.h5ad.gz
# from Supplementary file of GSE246613
# Format conversion needs to be performed using main.py

obj3_immune <- anndataR::read_h5ad(path = "GSE246613/GSE246613_PembroRT_immune_R100_final_updated.h5ad",as = "Seurat")

obj3_immune@meta.data <- obj3_immune@meta.data %>%
  dplyr::mutate(
    timepoint = dplyr::case_when(
      stringr::str_detect(patient_treatment, "_Base$")   ~ "pre",
      stringr::str_detect(patient_treatment, "_PD1$")    ~ "post_PD1",
      stringr::str_detect(patient_treatment, "_RTPD1$")  ~ "post_RTPD1",
      TRUE ~ NA_character_
    ),
    response = dplyr::case_when(
      as.character(pCR) == "R"  ~ "RE",
      as.character(pCR) == "NR" ~ "NR",
      TRUE ~ NA_character_
    ),
    patient = stringr::str_remove(
      as.character(cohort),
      "T[0-9]+$"
    ),
    tumor_id = stringr::str_extract(
      as.character(cohort),
      "T[0-9]+$"
    ),
    PATCS1 = dplyr::case_when(
      timepoint == "pre" ~ "basic",
      timepoint == "post_PD1" & response == "RE" ~ "active",
      timepoint == "post_RTPD1" & response == "RE" ~ "active",
      TRUE ~ NA_character_),
    PATCS2 = as.character(patient),
    PATCS3 = "GSE246613",
    sample_name = as.character(batch)
  )
obj3_immune <- subset(obj3_immune,subset = !is.na(PATCS1))
obj3_immune <- subset(obj3_immune,subset = celltype == "myeloid")
# delete myeloid_07 myeloid_08 myeloid_09 DC-like group
obj3_immune <- subset(obj3_immune,subset = !subcluster %in% c("myeloid_07","myeloid_08","myeloid_09"))
obj3_immune$subcluster <- droplevels(obj3_immune$subcluster)
obj3_immune$celltype <- "Macrophage"
qs_save(obj3_immune,"GSE246613/GSE246613_sub.qs",nthreads = 4)

# GSE207422 part1 ------

# download Supplementary file of GSE246613

count_file <- "GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
meta_file  <- "GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"

count_df <- data.table::fread(count_file,data.table = FALSE,check.names = FALSE)
gene_name <- make.unique(as.character(count_df[[1]]))
count_df[[1]] <- NULL
rownames(count_df) <- gene_name
count_matrix <- Matrix::Matrix(as.matrix(count_df),sparse = TRUE)
rm(count_df) ; gc()

cell_sample <- sub( "^(BD_immune[0-9]{2}).*$","\\1",colnames(count_matrix))
cell_metadata <- data.frame(sample_id = cell_sample,row.names = colnames(count_matrix),stringsAsFactors = FALSE)

GSE207422_obj <- CreateSeuratObject(
  counts = count_matrix,
  meta.data = cell_metadata,
  project = "GSE207422",
  min.cells = 3,
  min.features = 200
)

sample_info <- tribble(
  ~sample_id,      ~patient, ~timepoint, ~PATCS1,
  "BD_immune01",   "P01",    "pre",      "basic",
  "BD_immune02",   "P02",    "post",     NA_character_,
  "BD_immune03",   "P03",    "post",     "active",
  "BD_immune04",   "P04",    "post",     NA_character_,
  "BD_immune05",   "P05",    "pre",      "basic",
  "BD_immune06",   "P06",    "post",     "active",
  "BD_immune07",   "P07",    "post",     NA_character_,
  "BD_immune08",   "P08",    "pre",      "basic",
  "BD_immune09",   "P09",    "post",     NA_character_,
  "BD_immune10",   "P10",    "post",     NA_character_,
  "BD_immune11",   "P11",    "post",     "active",
  "BD_immune12",   "P12",    "post",     NA_character_,
  "BD_immune13",   "P13",    "post",     NA_character_,
  "BD_immune14",   "P14",    "post",     "active",
  "BD_immune15",   "P15",    "post",     NA_character_
)

match_index <- match(GSE207422_obj$sample_id,sample_info$sample_id)
stopifnot(!anyNA(match_index))

GSE207422_obj$patient <- sample_info$patient[match_index]
GSE207422_obj$timepoint <- sample_info$timepoint[match_index]
GSE207422_obj$PATCS1 <- sample_info$PATCS1[match_index]

qs_save(GSE207422_obj,"GSE207422/GSE207422_all.qs",nthreads = 4)

# GSE207422 part2 ------

GSE207422_obj[["percent.MT"]] <- PercentageFeatureSet(GSE207422_obj, pattern = "^MT-")
GSE207422_obj[["percent.HB"]] <- PercentageFeatureSet(GSE207422_obj, pattern = "^HBA|^HBB")

VlnPlot(GSE207422_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.MT", "percent.HB"),ncol = 2)

GSE207422_obj <- subset(GSE207422_obj, subset = nFeature_RNA > 200 & 
                          nFeature_RNA < 8000 & nCount_RNA > 500 & 
                          nCount_RNA < 60000 & percent.MT < 20 & percent.HB < 1) 

GSE207422_obj <- JoinLayers(object = GSE207422_obj,assay = "RNA")
GSE207422_obj <- NormalizeData(GSE207422_obj)
GSE207422_obj <- FindVariableFeatures(GSE207422_obj,selection.method = "vst",nfeatures = 2000)
GSE207422_obj <- ScaleData(GSE207422_obj,features = VariableFeatures(GSE207422_obj))
GSE207422_obj <- RunPCA(GSE207422_obj,features = VariableFeatures(object = GSE207422_obj))

ElbowPlot(GSE207422_obj, ndims = 50)
pc_sd <- GSE207422_obj[["pca"]]@stdev
pc_variance_percent <- pc_sd^2 / sum(pc_sd^2) * 100
pc_cumulative_percent <- cumsum(pc_variance_percent)
pc_cumulative_percent[30]

GSE207422_obj <- FindNeighbors(GSE207422_obj, dims = 1:30, reduction = "pca")
GSE207422_obj <- FindClusters(GSE207422_obj, resolution = 0.8, cluster.name = "unint_clusters")
GSE207422_obj <- RunUMAP(GSE207422_obj, dims = 1:30, reduction = "pca",reduction.name = "unint_UMAP")
GSE207422_obj <- GSE207422_obj %>% RunHarmony("sample_id", plot_convergence = TRUE)
GSE207422_obj <- FindNeighbors(GSE207422_obj, reduction = "harmony", dims = 1:30)
GSE207422_obj <- FindClusters(GSE207422_obj, resolution = 1.1, cluster.name = "harmony_clusters")
GSE207422_obj <- RunUMAP(GSE207422_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

MAEKER <- FindAllMarkers(GSE207422_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster_to_celltype <- setNames(
  ifelse(0:31 %in% c(3, 8, 11, 16, 17, 18, 27, 25),
         "Macrophage", "Other"),as.character(0:31))
annotated_celltype <- sapply(as.character(GSE207422_obj$harmony_clusters), function(cluster) {
  ifelse(cluster %in% names(cluster_to_celltype), 
         cluster_to_celltype[cluster], NA)  
})
names(annotated_celltype) <- colnames(GSE207422_obj)
celltype_levels <- unique(unname(cluster_to_celltype))
annotated_celltype <- factor(annotated_celltype, levels = celltype_levels)
GSE207422_obj$annotated_celltype <- annotated_celltype
GSE207422_obj <- subset(GSE207422_obj, subset = annotated_celltype == "Macrophage")
GSE207422_obj$annotated_celltype <- droplevels(GSE207422_obj$annotated_celltype)

GSE207422_obj@meta.data <- GSE207422_obj@meta.data %>%
  mutate(PATCS2 = as.character(patient),PATCS3 = "GSE207422")
GSE207422_obj <- subset(GSE207422_obj, subset = !is.na(PATCS1))
PATCS_cell_count <- GSE207422_obj@meta.data %>%
  dplyr::count(PATCS2,PATCS1,name = "n_cells") %>%
  dplyr::arrange(PATCS1, PATCS2)
qs_save(GSE207422_obj,"GSE207422/GSE207422_macro_sub.qs",nthreads = 4)

