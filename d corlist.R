
library(Seurat)
library(future)
library(future.apply)
library(dplyr)
library(qs2)

dir.create("corlist", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

# 1. load ------

allgenes <- read.csv(
  "KEGG/allgenes-KEGG.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

macro <- qs_read("input/zenodo.200524.sub.qs", nthreads = 6)

DefaultAssay(macro) <- "RNA"

gene_list <- unique(na.omit(as.character(allgenes$gene)))
gene_list <- gene_list[gene_list != ""]

available_genes <- intersect(gene_list, rownames(macro))
missing_genes <- setdiff(gene_list, rownames(macro))

# 2. fetch ------

exp.mat <- FetchData(
  object = macro,
  vars = available_genes
)

exp.matnew <- exp.mat[
  ,
  colSums(exp.mat, na.rm = TRUE) != 0,
  drop = FALSE
]

X <- as.matrix(exp.matnew)

# 3. gene pairs ------

pair_index <- combn(
  seq_len(ncol(X)),
  2
)

n_pairs <- ncol(pair_index)

# 4. multisession ------

workers <- 8
chunk_size <- 5000

pair_chunks <- split(
  seq_len(n_pairs),
  ceiling(seq_len(n_pairs) / chunk_size)
)

options(
  future.globals.maxSize = 60 * 1024^3
)

plan(
  multisession,
  workers = workers
)

# 5. Spearman rho , P value ------

cor_results_list <- future_lapply(
  pair_chunks,
  function(chunk_ids) {
    
    chunk_result <- vector(
      mode = "list",
      length = length(chunk_ids)
    )
    
    for (k in seq_along(chunk_ids)) {
      
      pair_id <- chunk_ids[k]
      
      gene1_index <- pair_index[1, pair_id]
      gene2_index <- pair_index[2, pair_id]
      
      x1 <- X[, gene1_index]
      x2 <- X[, gene2_index]
      
      keep_cells <- (
        is.finite(x1) &
          is.finite(x2) &
          x1 > 0 &
          x2 > 0
      )
      
      n_common <- sum(keep_cells)
      
      if (n_common <= 3) {
        next
      }
      
      x1_keep <- x1[keep_cells]
      x2_keep <- x2[keep_cells]
      
      if (
        length(unique(x1_keep)) < 2 ||
        length(unique(x2_keep)) < 2
      ) {
        next
      }
      
      cor_test_result <- tryCatch(
        suppressWarnings(
          cor.test(
            x = x1_keep,
            y = x2_keep,
            method = "spearman",
            exact = FALSE
          )
        ),
        error = function(e) NULL
      )
      
      if (is.null(cor_test_result)) {
        next
      }
      
      rho <- unname(cor_test_result$estimate)
      p_value <- cor_test_result$p.value
      
      if (
        is.na(rho) ||
        is.na(p_value)
      ) {
        next
      }
      
      chunk_result[[k]] <- data.frame(
        gene1 = colnames(X)[gene1_index],
        gene2 = colnames(X)[gene2_index],
        gene_pair = paste0(
          colnames(X)[gene1_index],
          ",",
          colnames(X)[gene2_index]
        ),
        spearman_rho = rho,
        p_value = p_value,
        n_coexpress = n_common,
        stringsAsFactors = FALSE
      )
    }
    
    dplyr::bind_rows(chunk_result)
  },
  future.seed = FALSE,
  future.globals = TRUE
)

plan(sequential)

# 6. save ------

corlist_long <- dplyr::bind_rows(
  cor_results_list
)

corlist_long$p_adj <- p.adjust(
  corlist_long$p_value,
  method = "BH"
)

corlist <- data.frame(
  rho = corlist_long$spearman_rho,
  p_value = corlist_long$p_value,
  p_adj = corlist_long$p_adj,
  n_coexpress = corlist_long$n_coexpress,
  row.names = make.unique(
    paste0(
      corlist_long$gene1,
      ".",
      corlist_long$gene2
    )
  ),
  check.names = FALSE
)

write.csv(corlist,file = "corlist/allgenes-corlist.csv")

