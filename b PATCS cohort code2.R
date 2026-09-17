
library(Seurat)
library(qs2)

GSE232240 <- qs_read("GSE232240/GSE232240_sub.qs",nthreads = 4)
GSE212217 <- qs_read("GSE212217/GSE212217_sub.qs",nthreads = 4)
GSE246613 <- qs_read("GSE246613/GSE246613_sub.qs",nthreads = 4)
GSE207422 <- qs_read("GSE207422/GSE207422_macro_sub.qs", nthreads = 4)

strip_to_counts_meta <- function(
    object,
    source_name,
    assay = "RNA",
    allow_X_as_counts = TRUE
) {
  #  assay
  if (!inherits(object, "Seurat")) {
    stop("obj ", source_name, " is not Seurat")
  }
  
  if (!assay %in% Assays(object)) {
    stop(
      paste0(
        "obj ", source_name,
        " assay not exist :", assay,
        ";current assays :",
        paste(Assays(object), collapse = ", ")
      )
    )
  }
  
  available_layers <- Layers(object[[assay]])
  
  message(
    source_name,
    " current RNA layers:",
    paste(available_layers, collapse = ", ")
  )
  
  #  counts 
  count_layers <- grep(
    pattern = "^counts($|\\.)",
    x = available_layers,
    value = TRUE
  )
  
  # 3. layer
  
  if (length(count_layers) > 1) {
    
    message(
      source_name,
      "  counts layers:",
      paste(count_layers, collapse = ", "),
      " "
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
      " with no counts layer, ",
      " convert X layer to counts "
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
        " obj ", source_name,
        "  ", assay,
        " assay with no counts or X layer",
        " current layers is:",
        paste(available_layers, collapse = ", ")
      )
    )
  }
  if (is.null(rownames(counts_matrix))) {
    stop("obj ", source_name, " counts_matrix with no rownames")
  }
  
  if (is.null(colnames(counts_matrix))) {
    stop("obj ", source_name, " counts_matrix with no colnames")
  }
  
  if (anyDuplicated(colnames(counts_matrix)) > 0) {
    stop("obj ", source_name, " counts_matrix with duplicated colnames")
  }
  
  missing_meta_cells <- setdiff(
    colnames(counts_matrix),
    rownames(object@meta.data)
  )
  
  if (length(missing_meta_cells) > 0) {
    stop(
      "obj ", source_name,
      "  ", length(missing_meta_cells),
      " with no metadata"
    )
  }
  
  metadata <- object@meta.data[
    colnames(counts_matrix),
    ,
    drop = FALSE
  ]
  
  metadata$merge_source <- source_name
  metadata$original_layer <- used_layer
  
  object_minimal <- CreateSeuratObject(
    counts = counts_matrix,
    assay = "RNA",
    meta.data = metadata,
    project = source_name,
    min.cells = 0,
    min.features = 0
  )
  
  if (!identical(
    colnames(object_minimal),
    rownames(object_minimal@meta.data)
  )) {
    stop("obj ", source_name, "")
  }
  return(object_minimal)
}

# ------

existing_object_names <- c("GSE232240","GSE212217","GSE246613","GSE207422")
P_list <- mget(existing_object_names, envir = .GlobalEnv)
P_min_list <- lapply(
  names(P_list),
  function(object_name) {
    message("Processing:", object_name)
    strip_to_counts_meta(
      object = P_list[[object_name]],
      source_name = object_name,
      assay = "RNA"
    )
  }
)
names(P_min_list) <- names(P_list)
rm(list = existing_object_names)
rm(P_list)
gc()



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

P_merge_tmp <- JoinLayers(object = P_merge_tmp,assay = "RNA")
P_merge_tmp <- NormalizeData(P_merge_tmp)
P_merge_tmp <- FindVariableFeatures(P_merge_tmp,selection.method = "vst",nfeatures = 2000)
P_merge_tmp <- ScaleData(P_merge_tmp,features = VariableFeatures(P_merge_tmp))
P_merge_tmp[["timepoint"]] <- P_merge_tmp[["PATCS1"]]
P_merge_tmp[["patient"]] <- paste(P_merge_tmp$PATCS3, P_merge_tmp$patient, sep = "_")
P_merge_tmp[["cohort"]] <- P_merge_tmp[["PATCS3"]]
P_merge_tmp@meta.data <- P_merge_tmp@meta.data[c("timepoint", "patient", "cohort")]

qs_save(P_merge_tmp,"input/PATCS.qs", nthreads = 4)

