
library(Seurat)
library(qs2)
library(readxl)

dir.create("1 PATCS", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

# 1. PATCS genelist ------

obj <- qs_read("input/PATCS.qs", nthreads = 6)

df <- read_excel("PATCS_genelist.xlsx", sheet = 1)

antigen_presenting <- as.character(na.omit(as.character(df$`Antigen presenting`)))
survival <- as.character(na.omit(as.character(df$Survival)))
cytokine <- as.character(na.omit(as.character(df$Cytokine)))
phago <- as.character(na.omit(as.character(df$Phagocytosis)))
recruit <- as.character(na.omit(as.character(df$Recruitment)))

length(phago)
length(antigen_presenting)
length(recruit)
length(cytokine)
length(survival)

obj <- NormalizeData(obj,normalization.method = "LogNormalize")

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

baseline_obj <- subset(obj, subset = timepoint == "basic")

activated_obj <- subset(obj, subset = timepoint == "active")


# 2. CRITIC ------
set.seed(123)
n_use <- 1000

generate_repeated_samples <- function(seu1, seu2, n_use = NULL, n_iter = 100, seed = 123) {
  cells1 <- colnames(seu1)
  cells2 <- colnames(seu2)
  
  if (is.null(n_use)) {
    n_use <- min(length(cells1), length(cells2))
  } else {
    n_use <- min(n_use, length(cells1), length(cells2))
  }
  
  res_list <- vector("list", n_iter)
  
  for (i in seq_len(n_iter)) {
    set.seed(seed + i)
    
    sampled1 <- sample(cells1, n_use, replace = FALSE)
    sampled2 <- sample(cells2, n_use, replace = FALSE)
    
    res_list[[i]] <- list(
      baseline_cells = sampled1,
      activated_cells = sampled2
    )
  }
  
  return(res_list)
}

sample_list <- generate_repeated_samples(
  baseline_obj,
  activated_obj,
  n_use = 1000,
  n_iter = 100,
  seed = 123
)

critic_weight <- function(X) {
  
  X <- as.matrix(X)
  
  X_norm <- apply(X, 2, function(x) {
    xmax <- max(x, na.rm = TRUE)
    xmin <- min(x, na.rm = TRUE)
    
    if (xmax == xmin) {
      rep(0, length(x))
    } else {
      (x - xmin) / (xmax - xmin)
    }
  })
  
  X_norm <- as.matrix(X_norm)
  colnames(X_norm) <- colnames(X)
  rownames(X_norm) <- rownames(X)
  
  sigma <- apply(X_norm, 2, sd, na.rm = TRUE)
  
  R <- cor(
    X_norm,
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  
  conflict <- sapply(seq_len(ncol(R)), function(j) {
    sum(1 - R[, j], na.rm = TRUE)
  })
  
  names(conflict) <- colnames(X)
  
  info <- sigma * conflict
  weights <- info / sum(info)
  
  return(list(
    weights = weights,
    sigma = sigma,
    correlation = R,
    conflict = conflict,
    info = info,
    X_norm = X_norm
  ))
}

run_critic_on_sampled_cells <- function(
    sample_item,
    baseline_obj,
    activated_obj,
    score_cols) {
  
  baseline_cells <- sample_item$baseline_cells
  activated_cells <- sample_item$activated_cells
  
  baseline_mat <- baseline_obj@meta.data[
    baseline_cells,
    score_cols,
    drop = FALSE
  ]
  
  activated_mat <- activated_obj@meta.data[
    activated_cells,
    score_cols,
    drop = FALSE
  ]
  
  stopifnot(
    nrow(baseline_mat) == nrow(activated_mat)
  )
  

  difference_mat <- activated_mat - baseline_mat
  
  rownames(difference_mat) <- paste0(
    "Active_",
    activated_cells,
    "__minus__Baseline_",
    baseline_cells
  )
  
  critic_res <- critic_weight(
    difference_mat
  )
  
  return(
    list(
      weights = critic_res$weights,
      difference_matrix = difference_mat
    )
  )
}

score_cols <- c(
  "antigen1",
  "cytokine1",
  "phagocytosis1",
  "recruitment1",
  "survival1"
)

critic_results <- lapply(sample_list, function(one_sample) {
  run_critic_on_sampled_cells(
    sample_item = one_sample,
    baseline_obj = baseline_obj,
    activated_obj = activated_obj,
    score_cols = score_cols
  )
})

weight_mat <- do.call(rbind, lapply(critic_results, function(x) x$weights))
weight_mat <- as.data.frame(weight_mat)

weight_summary <- data.frame(
  Function = colnames(weight_mat),
  Mean_Weight = colMeans(weight_mat, na.rm = TRUE),
  SD_Weight = apply(weight_mat, 2, sd, na.rm = TRUE),
  Median_Weight = apply(weight_mat, 2, median, na.rm = TRUE),
  Lower_25 = apply(weight_mat, 2, quantile, probs = 0.25, na.rm = TRUE),
  Upper_75 = apply(weight_mat, 2, quantile, probs = 0.75, na.rm = TRUE)
)

weight_summary <- weight_summary[
  order(weight_summary$Mean_Weight, decreasing = TRUE),
]

final_weights <- weight_summary$Mean_Weight
names(final_weights) <- weight_summary$Function
final_weights <- final_weights / sum(final_weights)
final_weights <- final_weights[score_cols]

write.csv(
  weight_mat,
  "1 PATCS/1.weight_mat.csv"
)

write.csv(
  weight_summary,
  "1 PATCS/2.weight_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    X = names(final_weights),
    x = as.numeric(final_weights)
  ),
  "1 PATCS/3.final_weights.csv",
  row.names = FALSE
)

