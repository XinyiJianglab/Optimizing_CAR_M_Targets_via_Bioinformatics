library(Seurat)
library(readxl)

dir.create("1 PATCS", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

P_merge <- readRDS("./P_merge.rds")

df <- read_excel("./PATCS_genelist.xlsx", sheet = 1)

antigen_presenting <- as.character(na.omit(as.character(df$`Antigen presenting`)))
survival <- as.character(na.omit(as.character(df$Survival)))
cytokine <- as.character(na.omit(as.character(df$Cytokine)))
phago <- as.character(na.omit(as.character(df$Phagocytosis)))
recruit <- as.character(na.omit(as.character(df$Recruitment)))

P_merge <- NormalizeData(P_merge,normalization.method = "LogNormalize")

P_merge <- AddModuleScore(
  P_merge,
  features = list(antigen_presenting),
  name = "antigen",
  seed = 123
)

P_merge <- AddModuleScore(
  P_merge,
  features = list(cytokine),
  name = "cytokine",
  seed = 123
)

P_merge <- AddModuleScore(
  P_merge,
  features = list(phago),
  name = "phagocytosis",
  seed = 123
)

P_merge <- AddModuleScore(
  P_merge,
  features = list(recruit),
  name = "recruitment",
  seed = 123
)

P_merge <- AddModuleScore(
  P_merge,
  features = list(survival),
  name = "survival",
  seed = 123
)

baseline_obj <- subset(P_merge, subset = PATCS1 == "basic")

activated_obj <- subset(P_merge, subset = PATCS1 == "active")

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
  data.frame(
    X = names(final_weights),
    x = as.numeric(final_weights)
  ),
  "./final_weights.csv",
  row.names = FALSE
)


#external validation
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

df <- read_excel("./PATCS_genelist.xlsx", sheet = 1)

antigen_presenting <- as.character(na.omit(as.character(df$`Antigen presenting`)))
survival <- as.character(na.omit(as.character(df$Survival)))
cytokine <- as.character(na.omit(as.character(df$Cytokine)))
phago <- as.character(na.omit(as.character(df$Phagocytosis)))
recruit <- as.character(na.omit(as.character(df$Recruitment)))

final_weights <- read.csv("./final_weights.csv")

patcs_ex <- AddModuleScore(patcs_ex, features = list(antigen_presenting), name = "antigen")
patcs_ex <- AddModuleScore(patcs_ex, features = list(cytokine), name = "cytokine")
patcs_ex <- AddModuleScore(patcs_ex, features = list(phago), name = "phagocytosis")
patcs_ex <- AddModuleScore(patcs_ex, features = list(recruit), name = "recruitment")
patcs_ex <- AddModuleScore(patcs_ex, features = list(survival), name = "survival")

score_cols <- c("antigen1", "cytokine1","phagocytosis1","recruitment1","survival1")
setdiff(
  score_cols,
  colnames(patcs_ex@meta.data)
)

score_matrix <- as.matrix(
  patcs_ex@meta.data[, score_cols, drop = FALSE]
)

storage.mode(score_matrix) <- "numeric"

score_matrix_z <- scale(score_matrix)

colnames(score_matrix_z) <- paste0(
  score_cols,
  "_z"
)

patcs_ex <- AddMetaData(
  object = patcs_ex,
  metadata = score_matrix_z
)

patcs_ex$PATCS <- as.numeric(
  score_matrix_z %*% weights
)

patcs_ex$PATCS <- scales::rescale(
  patcs_ex$PATCS,
  to = c(0, 1)
)

patient_response <- patcs_ex@meta.data %>%
  filter(
    !is.na(patient),
    !is.na(PATCS),
    Response %in% c("R", "NR")
  ) %>%
  group_by(patient) %>%
  summarise(
    n_macrophages = n(),
    PATCS = mean(PATCS, na.rm = TRUE),
    response = first(Response),
    .groups = "drop"
  ) %>%
  filter(n_macrophages >= 50) %>%
  mutate(
    response = factor(
      response,
      levels = c("NR", "R")
    )
  )
write.csv(patient_response,file = "./patcs_ex.csv")

mac <- readRDS("./macrophage_survival.rds")
mac <- AddModuleScore(mac, features = list(antigen_presenting), name = "antigen")
mac <- AddModuleScore(mac, features = list(cytokine), name = "cytokine")
mac <- AddModuleScore(mac, features = list(phago), name = "phagocytosis")
mac <- AddModuleScore(mac, features = list(recruit), name = "recruitment")
mac <- AddModuleScore(mac, features = list(survival), name = "survival")

score_cols <- c("antigen1", "cytokine1","phagocytosis1","recruitment1","survival1")
setdiff(
  score_cols,
  colnames(mac@meta.data)
)

score_matrix <- as.matrix(
  mac@meta.data[, score_cols, drop = FALSE]
)

storage.mode(score_matrix) <- "numeric"

score_matrix_z <- scale(score_matrix)

colnames(score_matrix_z) <- paste0(
  score_cols,
  "_z"
)

mac <- AddMetaData(
  object = mac,
  metadata = score_matrix_z
)

mac$PATCS <- as.numeric(
  score_matrix_z %*% weights
)

mac$PATCS <- scales::rescale(
  mac$PATCS,
  to = c(0, 1)
)

patient_patcs <- mac@meta.data %>%
  filter(
    !is.na(patient),
    !is.na(PATCS)
  ) %>%
  group_by(patient) %>%
  summarise(
    n_cells = n(),
    mean_PATCS = mean(PATCS, na.rm = TRUE),
    median_PATCS = median(PATCS, na.rm = TRUE),
    sd_PATCS = sd(PATCS, na.rm = TRUE),
    sem_PATCS = sd_PATCS / sqrt(n_cells),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_PATCS))

cutoff <- median(patient_patcs$mean_PATCS, na.rm = TRUE)

patient_patcs <- patient_patcs %>%
  mutate(
    PATCS_group = ifelse(
      mean_PATCS >= cutoff,
      "High",
      "Low"
    )
  )

write.csv(patient_patcs,file = "./patient_patcs.csv")

library(survival)
library(survminer)
library(readxl)

survival_data <- read_xlsx("./survival_data.xlsx")
survival_data$Overall_survival <- as.numeric(
  survival_data$Overall_survival
)

survival_data$Status <- as.numeric(
  survival_data$OS_Status
)

survival_data$group <- factor(
  survival_data$group,
  levels = c("Low", "High")
)

fit <- survfit(
  Surv(Overall_survival, OS_Status) ~ group,
  data = survival_data
)

p <- ggsurvplot(
  fit,
  data = survival_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3,
  palette = c(
    "Low" = "#3C8DBC",
    "High" = "#D9534F"
  ),
  legend.title = "PATCS",
  legend.labs = c("Low", "High"),
  xlab = "Overall survival (months)",
  ylab = "Overall survival probability",
  break.time.by = 10,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  surv.median.line = "hv",
  ggtheme = theme_classic(base_size = 14)
)

p
