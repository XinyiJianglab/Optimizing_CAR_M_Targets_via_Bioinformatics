library(KEGGREST)

get_kegg_genes <- function(pathway_id) {
  gsInfo <- KEGGREST::keggGet(pathway_id)[[1]]
  
  if (is.null(gsInfo$GENE) || length(gsInfo$GENE) == 0) {
    return(list(
      pathway_id = pathway_id,
      gene_field_length = 0,
      estimated_gene_number = 0,
      cleaned_gene_number = 0,
      genes = character(0)
    ))
  }
  
  gene_field_length <- length(gsInfo$GENE)
  estimated_gene_number <- gene_field_length / 2
  
  geneSetRaw <- vapply(
    strsplit(
      unname(gsInfo$GENE),
      split = ";",
      fixed = TRUE
    ),
    function(x) x[1],
    character(1)
  )
  
  geneSetRaw <- trimws(geneSetRaw)
  
  geneSetRaw_clean <- geneSetRaw[
    !is.na(geneSetRaw) &
      nzchar(geneSetRaw) &
      !grepl("\\[|\\]", geneSetRaw) &
      !grepl("^[0-9]+$", geneSetRaw)
  ]
  
  genes <- unique(geneSetRaw_clean)
  
  return(list(
    pathway_id = pathway_id,
    gene_field_length = gene_field_length,
    estimated_gene_number = estimated_gene_number,
    cleaned_gene_number = length(genes),
    genes = genes
  ))
}


pathway_ids <- c(
  C_type_lectin         = "hsa04625",
  RIG_I_like_receptor   = "hsa04622",
  Toll_like_receptor    = "hsa04620",
  Cytosolic_DNA_sensing = "hsa04623",
  TNF_signaling         = "hsa04668",
  Phagosome             = "hsa04145",
  NOD_like_receptor     = "hsa04621",
  Cytokine_interaction  = "hsa04060",
  Antigen_presentation  = "hsa04612"
)

kegg_results <- lapply(pathway_ids, get_kegg_genes)

kegg_gene_stats <- do.call(
  rbind,
  lapply(names(kegg_results), function(x) {
    res <- kegg_results[[x]]
    data.frame(
      pathway_name = x,
      pathway_id = res$pathway_id,
      gene_field_length = res$gene_field_length,
      estimated_gene_number = res$estimated_gene_number,
      cleaned_gene_number = res$cleaned_gene_number,
      stringsAsFactors = FALSE
    )
  })
)

kegg_gene_stats

gene_lists <- lapply( kegg_results, `[[`, "genes")
names(gene_lists) <- names(pathway_ids)
allgenes <- unique( unlist( gene_lists, use.names = FALSE ))
allgenes2 <- data.frame( gene = allgenes, stringsAsFactors = FALSE)
