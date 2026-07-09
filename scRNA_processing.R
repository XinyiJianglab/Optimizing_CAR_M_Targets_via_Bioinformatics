###### Load packages ------
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(tibble)
library(dplyr)
library(data.table)
library(Matrix)
library(qs)
library(harmony)
library(clustree)
library(readxl)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(ragg)
library(stringr)
library(anndataR)
library(future)
library(future.apply)

###### 1_GSE207422 data loading ------

count_file <- "./input/1_GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
meta_file  <- "./input/1_GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"

count_df <- data.table::fread(
  count_file,
  data.table = FALSE,
  check.names = FALSE
)

gene_name <- make.unique(as.character(count_df[[1]]))
count_df[[1]] <- NULL
rownames(count_df) <- gene_name

count_matrix <- Matrix::Matrix(
  as.matrix(count_df),
  sparse = TRUE
)

rm(count_df) ; gc()

cell_sample <- sub(
  "^(BD_immune[0-9]{2}).*$",
  "\\1",
  colnames(count_matrix)
)

table(cell_sample)

cell_metadata <- data.frame(
  sample_id = cell_sample,
  row.names = colnames(count_matrix),
  stringsAsFactors = FALSE
)

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
  "BD_immune12",   "P3",    "post",     NA_character_,
  "BD_immune13",   "P13",    "post",     NA_character_,
  "BD_immune14",   "P4",    "post",     "active",
  "BD_immune15",   "P15",    "post",     NA_character_
)

match_index <- match(
  GSE207422_obj$sample_id,
  sample_info$sample_id
)

stopifnot(!anyNA(match_index))

GSE207422_obj$patient <- sample_info$patient[match_index]
GSE207422_obj$timepoint <- sample_info$timepoint[match_index]
GSE207422_obj$PATCS1 <- sample_info$PATCS1[match_index]

qs_save(
  GSE207422_obj,
  "./input/1_GSE207422/1_GSE207422_all.qs",
  nthreads = 4
)


###### 1_GSE207422 preprocessing ------

# QC and doublet removal were performed in the original study.
# 92,330 cells were retained after initial filtering.

SCRNA <- qs_read(
  "./input/1_GSE207422/1_GSE207422_all.qs",
  nthreads = 4
)

SCRNA[["percent.MT"]] <- PercentageFeatureSet(
  SCRNA,
  pattern = "^MT-"
)

SCRNA[["percent.HB"]] <- PercentageFeatureSet(
  SCRNA,
  pattern = "^HBA|^HBB"
)

VlnPlot(
  SCRNA,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.MT",
    "percent.HB"
  ),
  ncol = 2
)

SCRNA <- subset(
  SCRNA,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 8000 &
    nCount_RNA > 500 &
    nCount_RNA < 60000 &
    percent.MT < 20 &
    percent.HB < 1
)


#######################################################################-

options(future.globals.maxSize = 80 * 1024^3)

SCRNA <- JoinLayers(
  object = SCRNA,
  assay = "RNA"
)

SCRNA <- NormalizeData(SCRNA)

SCRNA <- FindVariableFeatures(
  SCRNA,
  selection.method = "vst",
  nfeatures = 2000
)

SCRNA <- ScaleData(
  SCRNA,
  features = VariableFeatures(SCRNA)
)

SCRNA <- RunPCA(
  SCRNA,
  features = VariableFeatures(object = SCRNA)
)

ElbowPlot(SCRNA, ndims = 50)

pc_sd <- SCRNA[["pca"]]@stdev

pc_variance_percent <- 
  pc_sd^2 / sum(pc_sd^2) * 100

pc_cumulative_percent <- 
  cumsum(pc_variance_percent)

pc_cumulative_percent[30]


SCRNA <- FindNeighbors(
  SCRNA,
  dims = 1:30,
  reduction = "pca"
)

SCRNA <- FindClusters(
  SCRNA,
  resolution = 0.8,
  cluster.name = "unint_clusters"
)

SCRNA <- RunUMAP(
  SCRNA,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "unint_UMAP"
)


SCRNA <- SCRNA %>%
  RunHarmony(
    "sample_id",
    plot_convergence = TRUE
  )

SCRNA <- FindNeighbors(
  SCRNA,
  reduction = "harmony",
  dims = 1:30
)

SCRNA <- FindClusters(
  SCRNA,
  resolution = seq(1.3,1.5,by=0.1)
)

clustree(SCRNA)


SCRNA <- FindClusters(
  SCRNA,
  resolution = 1.1,
  cluster.name = "harmony_clusters"
)

SCRNA <- RunUMAP(
  SCRNA,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap.harmony"
)


p1 <- DimPlot(
  SCRNA,
  reduction = "unint_UMAP",
  group.by = c("sample_id", "unint_clusters"),
  combine = TRUE,
  label.size = 2,
  label = TRUE
)

p2 <- DimPlot(
  SCRNA,
  reduction = "umap.harmony",
  group.by = c("sample_id", "harmony_clusters"),
  combine = TRUE,
  label.size = 2,
  label = TRUE
)


agg_png(
  "./input/1_GSE207422/output/2.umap1.png",
  width = 12,
  height = 5,
  units = "in",
  res = 1200
)

print(p1)

dev.off()


agg_png(
  "./input/1_GSE207422/output/2.umap2.png",
  width = 12,
  height = 5,
  units = "in",
  res = 1200
)

print(p2)

dev.off()


qs_save(
  SCRNA,
  './input/1_GSE207422/output/2.SCRNA_filtered_harmony.qs',
  nthreads = 4
)


#######################################################################-

SCRNA <- qs_read(
  './input/1_GSE207422/output/2.SCRNA_filtered_harmony.qs',
  nthreads = 4
)

marker_file <- "./input/1.marker.xlsx"

marker_data <- read_excel(
  marker_file,
  col_names = TRUE
)


cell_types <- colnames(marker_data)

marker_list <- lapply(
  cell_types,
  function(cell_type) {
    markers <- marker_data[[cell_type]]
    markers <- markers[!is.na(markers)]
    return(markers)
  }
)

names(marker_list) <- cell_types

# Define marker visualization function
plot_cell_markers <- function(cell_type, markers, SCRNA) {
  
  DefaultAssay(SCRNA) <- "RNA"
  
  plots <- list()
  
  harmony_cluster_plot <- DimPlot(
    SCRNA,
    reduction = "umap.harmony",
    group.by = "harmony_clusters",
    label = TRUE
  ) +
    ggtitle("harmony_clusters")
  
  plots <- c(plots, list(harmony_cluster_plot))
  
  
  for (marker in markers) {
    
    if (!marker %in% rownames(SCRNA[["RNA"]])) {
      warning(
        paste(
          "Gene",
          marker,
          "not found in the dataset. Skipping."
        )
      )
      next
    }
    
    harmony_plot <- FeaturePlot(
      SCRNA,
      features = marker,
      reduction = "umap.harmony",
      raster = FALSE,
      order = FALSE
    ) +
      ggtitle(
        paste(
          "harmony",
          marker,
          sep = "_"
        )
      )
    
    plots <- c(plots, list(harmony_plot))
  }
  
  return(plots)
}



# Generate UMAP marker plots for all cell types
for (i in seq_along(cell_types)) {
  
  cell_type <- cell_types[i]
  markers <- marker_list[[cell_type]]
  
  plots <- plot_cell_markers(
    cell_type,
    markers,
    SCRNA
  )
  
  
  if (length(plots) < 2) {
    warning(
      paste(
        "No valid plots generated for cell type",
        cell_type,
        ". Skipping."
      )
    )
    next
  }
  
  
  num_plots <- length(plots)
  plot_height <- 5 * ceiling(num_plots / 2)
  
  
  png_filename <- paste0(
    "./input/1_GSE207422/output/3.harmony/1.umap_harmony_",
    cell_type,
    ".png"
  )
  
  
  agg_png(
    png_filename,
    width = 12,
    height = plot_height,
    units = "in",
    res = 600
  )
  
  grid.arrange(
    grobs = plots,
    ncol = 2
  )
  
  dev.off()
}



RES <- FindMarkers(
  SCRNA,
  ident.1 = 12,
  ident.2 = NULL,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


#######################################################################-

# Define macrophage clusters based on harmony clustering
cluster_to_celltype <- setNames(
  ifelse(
    0:31 %in% c(3,8,11,16,17,18,27,25),
    "Macrophage",
    "Other"
  ),
  as.character(0:31)
)


cluster_to_celltype


# Check whether all clusters are annotated
unmapped_clusters <- setdiff(
  unique(SCRNA$harmony_clusters),
  names(cluster_to_celltype)
)

if (length(unmapped_clusters) > 0) {
  warning(
    paste(
      "Unmapped clusters:",
      paste(unmapped_clusters, collapse = ", ")
    )
  )
}


# Assign cell type annotation
annotated_celltype <- sapply(
  as.character(SCRNA$harmony_clusters),
  function(cluster) {
    
    ifelse(
      cluster %in% names(cluster_to_celltype),
      cluster_to_celltype[cluster],
      NA
    )
  }
)


names(annotated_celltype) <- colnames(SCRNA)


celltype_levels <- unique(
  unname(cluster_to_celltype)
)


annotated_celltype <- factor(
  annotated_celltype,
  levels = celltype_levels
)


SCRNA$annotated_celltype <- annotated_celltype



# Cell type colors
celltype_palette <- c(
  "Macrophage" = "#D55E00",
  "Other" = "#999999"
)



# Plot annotated UMAP
p_umap_annotation <- DimPlot(
  SCRNA,
  reduction = "umap.harmony",
  group.by = "annotated_celltype",
  cols = celltype_palette,
  label = TRUE,
  pt.size = 0.8,
  raster = FALSE
) +
  ggtitle("Annotated UMAP (Harmony)")


agg_png(
  "./input/1_GSE207422/output/4.Annotated_UMAP.png",
  width = 8,
  height = 6,
  units = "in",
  res = 600
)

print(p_umap_annotation)

dev.off()


qs_save(
  SCRNA,
  './input/1_GSE207422/output/3.SCRNA_filtered_harmony_annotated.qs'
)


SCRNA <- subset(
  SCRNA,
  subset = annotated_celltype == "Macrophage"
)

SCRNA$annotated_celltype <- droplevels(
  SCRNA$annotated_celltype
)


qs_save(
  SCRNA,
  './input/1_GSE207422/output/1.GSE207422_macro_sub.qs'
)



###### Extract GSE207422 macrophage dataset ------

SCRNA@meta.data <- SCRNA@meta.data %>%
  mutate(
    PATCS2 = as.character(patient),
    PATCS3 = "GSE207422"
  )


SCRNA <- subset(
  SCRNA,
  subset = !is.na(PATCS1)
)



PATCS_cell_count <- SCRNA@meta.data %>%
  dplyr::count(
    PATCS2,
    PATCS1,
    name = "n_cells"
  ) %>%
  dplyr::arrange(
    PATCS1,
    PATCS2
  )


PATCS_cell_count



qs_save(
  SCRNA,
  "./input/1_GSE207422/output/1.GSE207422_macro_sub_2.qs",
  nthreads = 4
)



###### 2_GSE232240 ------


count_file <- "./input/2_GSE232240/GSE232240_Count_data_IMCISION.txt"
meta_file  <- "./input/2_GSE232240/GSE232240_Meta_data_IMCISION.txt"


first_line <- readLines(
  count_file,
  n = 1
)


cell_names <- strsplit(
  first_line,
  "\t",
  fixed = TRUE
)[[1]]


count_df <- data.table::fread(
  count_file,
  sep = "\t",
  header = FALSE,
  skip = 1,
  fill = TRUE,
  data.table = FALSE,
  check.names = FALSE
)


gene_names <- as.character(count_df[[1]])

count_df[[1]] <- NULL


length(cell_names)

ncol(count_df)


colnames(count_df) <- cell_names


count_matrix <- as.matrix(count_df)

storage.mode(count_matrix) <- "numeric"


rownames(count_matrix) <- make.unique(
  gene_names
)


count_matrix <- Matrix::Matrix(
  count_matrix,
  sparse = TRUE
)


P2 <- CreateSeuratObject(
  counts = count_matrix,
  project = "GSE232240",
  min.cells = 0,
  min.features = 0
)



meta <- data.table::fread(
  meta_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)


identical(
  meta$cell_id,
  colnames(P2)
)


P2 <- AddMetaData(
  object = P2,
  metadata = meta
)


P2_sub <- subset(
  P2,
  subset = mc_group == "Mono-macro"
)



P2_sub@meta.data <- P2_sub@meta.data %>%
  mutate(
    PATCS1 = case_when(
      timepoint == "pre" ~ "basic",
      timepoint == "post" & response == "RE" ~ "active",
      TRUE ~ NA_character_
    ),
    
    PATCS2 = as.character(patient),
    
    PATCS3 = "GSE232240"
  )



P2_sub <- subset(
  P2_sub,
  subset = !is.na(PATCS1)
)


PATCS_cell_count <- P2_sub@meta.data %>%
  count(
    PATCS2,
    PATCS1,
    name = "n_cells"
  ) %>%
  arrange(
    PATCS1,
    PATCS2
  )


PATCS_cell_count



response_patient_list <- P2_sub@meta.data %>%
  filter(
    response %in% c("NR","RE")
  ) %>%
  select(
    response,
    patient
  ) %>%
  distinct() %>%
  arrange(
    response,
    patient
  )


response_patient_list


table(P2_sub$PATCS1)



qs_save(
  P2_sub,
  "./input/2_GSE232240/2_GSE232240_sub.qs",
  nthreads = 4
)

###### 3_GSE212217 ------

# Extract macrophage/monocyte populations
P3_sub <- subset(
  P3,
  cells = rownames(P3@meta.data)[
    grepl(
      "TAM|mono",
      P3@meta.data$finalIdent,
      ignore.case = TRUE
    )
  ]
)


table(P3_sub$finalIdent)



P3_sub@meta.data <- P3_sub@meta.data %>%
  dplyr::mutate(
    
    response_simple = dplyr::case_when(
      clinical %in% c("epiR", "mutR") ~ "RE",
      clinical == "NR"                ~ "NR",
      clinical == "Healthy"           ~ "Healthy",
      TRUE                            ~ NA_character_
    ),
    
    
    PATCS1 = dplyr::case_when(
      
      timepointBinary == "pre" &
        clinical != "Healthy" ~ "basic",
      
      timepointBinary == "post" &
        clinical %in% c("epiR", "mutR") ~ "active",
      
      TRUE ~ NA_character_
    ),
    
    
    PATCS2 = ifelse(
      grepl(
        "^[0-9]+$",
        as.character(patient)
      ),
      paste0(
        "P",
        as.character(patient)
      ),
      as.character(patient)
    ),
    
    
    PATCS3 = "GSE212217"
  )



P3_sub <- subset(
  P3_sub,
  subset = !is.na(PATCS1)
)



PATCS_cell_count <- P3_sub@meta.data %>%
  dplyr::count(
    PATCS3,
    PATCS2,
    PATCS1,
    response_simple,
    clinical,
    name = "n_cells"
  ) %>%
  dplyr::arrange(
    PATCS1,
    PATCS2
  )


PATCS_cell_count



response_patient_list <- P3_sub@meta.data %>%
  dplyr::filter(
    response_simple %in% c("NR","RE")
  ) %>%
  dplyr::select(
    response = response_simple,
    clinical,
    patient = PATCS2
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(
    response,
    clinical,
    patient
  )


response_patient_list



qs_save(
  P3_sub,
  "./input/3_GSE212217/3_GSE212217_sub.qs",
  nthreads = 4
)




###### 4_GSE246613 ------


data_dir <- "./input/4_GSE246613"


file_all <- file.path(
  data_dir,
  "GSE246613_combined_RTPDv4_scvi_celltypist.h5ad"
)


file_immune <- file.path(
  data_dir,
  "GSE246613_PembroRT_immune_R100_final_updated.h5ad"
)



P4_all <- anndataR::read_h5ad(
  path = file_all,
  as = "Seurat"
)


P4_immune <- anndataR::read_h5ad(
  path = file_immune,
  as = "Seurat"
)



P4_immune@meta.data <- P4_immune@meta.data %>%
  dplyr::mutate(
    
    
    timepoint = dplyr::case_when(
      
      stringr::str_detect(
        patient_treatment,
        "_Base$"
      ) ~ "pre",
      
      stringr::str_detect(
        patient_treatment,
        "_PD1$"
      ) ~ "post_PD1",
      
      stringr::str_detect(
        patient_treatment,
        "_RTPD1$"
      ) ~ "post_RTPD1",
      
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
      
      timepoint %in% c(
        "post_PD1",
        "post_RTPD1"
      ) &
        response == "RE" ~ "active",
      
      TRUE ~ NA_character_
    ),
    
    
    
    PATCS2 = as.character(patient),
    
    
    PATCS3 = "GSE246613",
    
    
    sample_name = as.character(batch)
    
  )



P4_immune <- subset(
  P4_immune,
  subset = !is.na(PATCS1)
)



# Remove dendritic cell clusters according to the original annotation
P4_immune$celltype <- "Macrophage"


P4_immune <- subset(
  P4_immune,
  subset = !subcluster %in%
    c(
      "myeloid_07",
      "myeloid_08",
      "myeloid_09"
    )
)


P4_immune$subcluster <- droplevels(
  P4_immune$subcluster
)



qs_save(
  P4_immune,
  "./input/4_GSE246613/4_GSE246613_sub.qs",
  nthreads = 4
)





##### Load processed macrophage datasets ------


P2 <- qs_read(
  "./input/2_GSE232240/2_GSE232240_sub.qs",
  nthreads = 4
)


P3 <- qs_read(
  "./input/3_GSE212217/3_GSE212217_sub.qs",
  nthreads = 4
)


P4 <- qs_read(
  "./input/4_GSE246613/4_GSE246613_sub.qs",
  nthreads = 4
)


P1 <- qs_read(
  './input/1_GSE207422/output/1.GSE207422_macro_sub_2.qs',
  nthreads = 4
)



# Check RNA layers
Layers(P2[["RNA"]])

Layers(P3[["RNA"]])

Layers(P4[["RNA"]])

Layers(P1[["RNA"]])




##### Function for extracting counts and metadata ------


strip_to_counts_meta <- function(
    object,
    source_name,
    assay = "RNA",
    allow_X_as_counts = TRUE
) {
  
  
  if (!inherits(object, "Seurat")) {
    stop(
      "Object ",
      source_name,
      " is not a Seurat object."
    )
  }
  
  
  if (!assay %in% Assays(object)) {
    stop(
      paste0(
        "Object ",
        source_name,
        " does not contain assay: ",
        assay
      )
    )
  }
  
  
  available_layers <- Layers(object[[assay]])
  
  
  message(
    source_name,
    " RNA layers: ",
    paste(
      available_layers,
      collapse = ", "
    )
  )
  
  
  count_layers <- grep(
    pattern = "^counts($|\\.)",
    x = available_layers,
    value = TRUE
  )
  
  
  if (length(count_layers) > 1) {
    
    
    message(
      source_name,
      " contains multiple count layers. Joining layers."
    )
    
    
    object <- JoinLayers(
      object = object,
      assay = assay,
      layers = count_layers,
      new = "counts_joined"
    )
    
    
    counts_matrix <- LayerData(
      object = object,
      assay = assay,
      layer = "counts_joined"
    )
    
    
    used_layer <- "counts_joined"
    
    
  } else if (length(count_layers) == 1) {
    
    
    used_layer <- count_layers[1]
    
    
    counts_matrix <- LayerData(
      object = object,
      assay = assay,
      layer = used_layer
    )
    
    
  } else if (
    isTRUE(allow_X_as_counts) &&
    "X" %in% available_layers
  ) {
    
    
    message(
      source_name,
      " no counts layer detected. Using X layer."
    )
    
    
    used_layer <- "X"
    
    
    counts_matrix <- LayerData(
      object = object,
      assay = assay,
      layer = "X"
    )
    
    
  } else {
    
    stop(
      paste0(
        "No count matrix found in object ",
        source_name
      )
    )
    
  }
  
  
  # Validate expression matrix
  
  if (is.null(rownames(counts_matrix))) {
    stop(
      "Object ",
      source_name,
      " has no gene names."
    )
  }
  
  
  if (is.null(colnames(counts_matrix))) {
    stop(
      "Object ",
      source_name,
      " has no cell names."
    )
  }
  
  
  if (anyDuplicated(colnames(counts_matrix)) > 0) {
    stop(
      "Object ",
      source_name,
      " contains duplicated cell names."
    )
  }
  
  
  missing_meta_cells <- setdiff(
    colnames(counts_matrix),
    rownames(object@meta.data)
  )
  
  
  if (length(missing_meta_cells) > 0) {
    stop(
      "Object ",
      source_name,
      " contains cells missing from metadata."
    )
  }
  
  
  
  # Extract metadata according to expression matrix order
  
  metadata <- object@meta.data[
    colnames(counts_matrix),
    ,
    drop = FALSE
  ]
  
  
  metadata$merge_source <- source_name
  
  metadata$original_layer <- used_layer
  
  
  
  # Create minimal Seurat object
  
  object_minimal <- CreateSeuratObject(
    counts = counts_matrix,
    assay = "RNA",
    meta.data = metadata,
    project = source_name,
    min.cells = 0,
    min.features = 0
  )
  
  
  
  # Check cell and metadata consistency
  
  if (!identical(
    colnames(object_minimal),
    rownames(object_minimal@meta.data)
  )) {
    
    stop(
      "Cell names and metadata order are inconsistent after reconstruction."
    )
    
  }
  
  
  
  message(
    source_name,
    " processed: ",
    nrow(object_minimal),
    " genes, ",
    ncol(object_minimal),
    " cells; layer used: ",
    used_layer
  )
  
  
  return(object_minimal)
}






##### Extract minimal Seurat objects ------


# Define target object names

target_object_names <- paste0(
  "P",
  1:100
)



# Identify existing objects

existing_object_names <- target_object_names[
  vapply(
    target_object_names,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]



if (length(existing_object_names) == 0) {
  stop(
    "No Seurat objects found."
  )
}



message(
  "Found ",
  length(existing_object_names),
  " objects: ",
  paste(
    existing_object_names,
    collapse = ", "
  )
)



# Store objects into list

P_list <- mget(
  existing_object_names,
  envir = .GlobalEnv
)



# Check Seurat objects

is_seurat <- vapply(
  P_list,
  inherits,
  logical(1),
  what = "Seurat"
)



if (any(!is_seurat)) {
  
  stop(
    "The following objects are not Seurat objects: ",
    paste(
      names(P_list)[!is_seurat],
      collapse = ", "
    )
  )
}




# Extract counts and metadata

P_min_list <- lapply(
  names(P_list),
  function(object_name) {
    
    
    message(
      "Processing: ",
      object_name
    )
    
    
    strip_to_counts_meta(
      object = P_list[[object_name]],
      source_name = object_name,
      assay = "RNA"
    )
    
  }
)



names(P_min_list) <- names(P_list)




# Summarize extracted objects

object_summary <- do.call(
  rbind,
  lapply(
    names(P_min_list),
    function(object_name) {
      
      data.frame(
        object = object_name,
        genes = nrow(P_min_list[[object_name]]),
        cells = ncol(P_min_list[[object_name]]),
        stringsAsFactors = FALSE
      )
      
    }
  )
)



print(object_summary)




# Remove original large objects

rm(list = existing_object_names)

rm(P_list)

gc()






##### Merge datasets ------


if (length(P_min_list) == 1) {
  
  
  P_merge_tmp <- P_min_list[[1]]
  
  
} else {
  
  
  P_merge_tmp <- merge(
    x = P_min_list[[1]],
    y = P_min_list[-1],
    add.cell.ids = names(P_min_list),
    project = "PATCS_cohorts",
    merge.data = FALSE,
    merge.dr = FALSE
  )
  
}





##### Basic processing and save ------


options(
  future.globals.maxSize = 80 * 1024^3
)



P_merge_tmp <- JoinLayers(
  object = P_merge_tmp,
  assay = "RNA"
)



P_merge_tmp <- NormalizeData(
  P_merge_tmp
)



P_merge_tmp <- FindVariableFeatures(
  P_merge_tmp,
  selection.method = "vst",
  nfeatures = 2000
)



P_merge_tmp <- ScaleData(
  P_merge_tmp,
  features = VariableFeatures(P_merge_tmp)
)



P_merge_tmp <- RunPCA(
  P_merge_tmp,
  features = VariableFeatures(object = P_merge_tmp)
)



qs_save(
  P_merge_tmp,
  "./input/PATCS/P_merge.qs",
  nthreads = 4
)





##### Differential expression analysis ------


P_merge1 <- qs_read(
  "./input/PATCS/P_merge.qs",
  nthreads = 4
)



Idents(P_merge1) <- "PATCS1"



DEG1 <- FindMarkers(
  P_merge1,
  ident.1 = "active",
  ident.2 = "basic",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0
)


