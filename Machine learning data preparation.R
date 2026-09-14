library(Seurat)
library(readxl)
library(Matrix)
library(qs2)

dir.create("3 MLlist and data", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

macro_basic <- readRDS("./pan_cancer_macrophage.rds")
allgenes <- read.csv("./allgenes.csv")

genelist0 <- intersect(
  allgenes$gene,
  rownames(macro_basic)
)

genelist0 <- genelist0[
  Matrix::rowSums(
    macro_basic@assays$RNA$counts[
      genelist0,
      ,
      drop = FALSE
    ]
  ) > 0
]

PATCS_genelist <- read_excel(
  "./PATCS_genelist.xlsx",
  sheet = 1
)

PATCS_genes <- unique(
  as.character(
    na.omit(
      unlist(
        PATCS_genelist,
        use.names = FALSE
      )
    )
  )
)

genelist0 <- genelist0[
  !genelist0 %in% PATCS_genes
]

deg <- read.csv("./deg.csv")
deg_filtered <- deg[deg$avg_log2FC > 1.2 & deg$p_val_adj<0.05, ]
deg_genes <- unique(deg_filtered$X)

genelist_filt <- intersect(genelist0 , deg_genes)

length(genelist_filt)

corlist2 <- read.csv(
  "./corlist.csv",
  stringsAsFactors = FALSE
)

corlist2 <- corlist2[
  corlist2$rho > 0.4 &
    corlist2$p_adj < 0.05 &
    corlist2$n_coexpress > ncol(macro_basic) * 0.01,
]

pair_split <- strsplit(
  as.character(corlist2$X),
  ".",
  fixed = TRUE
)

corlist2$gene1 <- sapply(
  pair_split,
  `[`,
  1
)

corlist2$gene2 <- sapply(
  pair_split,
  `[`,
  2
)

corlist_filt <- corlist2[
  corlist2$gene1 %in% genelist_filt &
    corlist2$gene2 %in% genelist_filt,
]

nrow(corlist_filt)

set.seed(123)
split_ratio <- 0.7

zenodo_beforeharmony <- macro_basic

n_cells <- ncol(zenodo_beforeharmony)
n_train <- floor(n_cells * split_ratio)

train_indices <- sample(
  n_cells,
  n_train,
  replace = FALSE
)

train_set <- zenodo_beforeharmony[, train_indices]
test_set <- zenodo_beforeharmony[, -train_indices]

df <- PATCS_genelist

antigen_presenting <- as.character(
  na.omit(
    as.character(
      df$`Antigen presenting`
    )
  )
)

survival <- as.character(
  na.omit(
    as.character(
      df$Survival
    )
  )
)

cytokine <- as.character(
  na.omit(
    as.character(
      df$Cytokine
    )
  )
)

phago <- as.character(
  na.omit(
    as.character(
      df$Phagocytosis
    )
  )
)

recruit <- as.character(
  na.omit(
    as.character(
      df$Recruitment
    )
  )
)

score_cols <- c(
  "antigen1",
  "cytokine1",
  "phagocytosis1",
  "recruitment1",
  "survival1"
)

weights_df <- read.csv(
  "./final_weights.csv"
)

weights_df <- weights_df[
  match(
    score_cols,
    weights_df$X
  ),
]

weights <- weights_df$x
names(weights) <- weights_df$X

stable_pair_seed <- function(gene1, gene2, base_seed = 123) {

  pair_genes <- sort(
    c(
      as.character(gene1),
      as.character(gene2)
    )
  )

  pair_key <- paste(
    pair_genes,
    collapse = "."
  )

  ints <- utf8ToInt(
    enc2utf8(pair_key)
  )

  mod <- 2147483646
  h <- as.numeric(base_seed) %% mod

  for (v in ints) {
    h <- (h * 131 + as.numeric(v)) %% mod
  }

  as.integer(h + 1)
}

run_scoring <- function(obj, filtered_elements_df,
                        antigen_presenting, cytokine,
                        phago, recruit, survival, weights) {

  obj <- NormalizeData(
    obj,
    verbose = FALSE
  )

  score_names <- as.character(
    filtered_elements_df$X
  )

  for (i in seq_len(nrow(filtered_elements_df))) {

    gene1_i <- as.character(
      filtered_elements_df$gene1[i]
    )

    gene2_i <- as.character(
      filtered_elements_df$gene2[i]
    )

    pair_genes_i <- sort(
      c(
        gene1_i,
        gene2_i
      )
    )

    pair_seed <- stable_pair_seed(
      gene1_i,
      gene2_i,
      base_seed = 123
    )

    tmp_obj <- AddModuleScore(
      obj,
      features = list(pair_genes_i),
      name = "TMPPAIR",
      seed = pair_seed
    )

    obj@meta.data[[score_names[i]]] <- tmp_obj@meta.data[["TMPPAIR1"]]

    rm(tmp_obj)
  }

  obj <- AddModuleScore(
    obj,
    features = list(antigen_presenting),
    name = "antigen",
    seed = 123
  )

  obj <- AddModuleScore(
    obj,
    features = list(cytokine),
    name = "cytokine",
    seed = 123
  )

  obj <- AddModuleScore(
    obj,
    features = list(phago),
    name = "phagocytosis",
    seed = 123
  )

  obj <- AddModuleScore(
    obj,
    features = list(recruit),
    name = "recruitment",
    seed = 123
  )

  obj <- AddModuleScore(
    obj,
    features = list(survival),
    name = "survival",
    seed = 123
  )

  obj$Integrated_Score <- as.numeric(
    as.matrix(
      obj@meta.data[
        ,
        score_cols,
        drop = FALSE
      ]
    ) %*%
      weights[score_cols]
  )

  keep_cols <- c(
    score_names,
    "Integrated_Score"
  )

  obj@meta.data <- obj@meta.data[
    ,
    keep_cols,
    drop = FALSE
  ]

  return(obj)
}


train_object <- run_scoring(
  train_set,
  corlist_filt,
  antigen_presenting,
  cytokine,
  phago,
  recruit,
  survival,
  weights
)

test_object <- run_scoring(
  test_set,
  corlist_filt,
  antigen_presenting,
  cytokine,
  phago,
  recruit,
  survival,
  weights
)

train_genematrix <- FetchData(
  train_object,
  c(
    genelist_filt,
    corlist_filt$X,
    "Integrated_Score"
  )
)

test_genematrix <- FetchData(
  test_object,
  c(
    genelist_filt,
    corlist_filt$X,
    "Integrated_Score"
  )
)

qs_save(
  train_genematrix,
  "./train_genematrix.qs",
  nthreads = 4
)

qs_save(
  test_genematrix,
  "./test_genematrix.qs",
  nthreads = 4
)
