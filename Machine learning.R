library(qs2)
library(ranger)
library(dplyr)

dir.create("4 ML", recursive = TRUE, showWarnings = FALSE)

train_genematrix <- qs_read(
  "./train_genematrix.qs",
  nthreads = 4
)

test_genematrix <- qs_read(
  "./test_genematrix.qs",
  nthreads = 4
)

set.seed(123)
split_ratio <- 0.7

X <- as.matrix(
  train_genematrix[
    ,
    setdiff(colnames(train_genematrix), "Integrated_Score")
  ]
)
Y <- as.matrix(train_genematrix$Integrated_Score)

dim(X)
dim(Y)

X_train <- as.matrix(
  train_genematrix[
    ,
    setdiff(colnames(train_genematrix), "Integrated_Score")
  ]
)

X_test <- as.matrix(
  test_genematrix[
    ,
    setdiff(colnames(test_genematrix), "Integrated_Score")
  ]
)


Y_train <- as.matrix(train_genematrix$Integrated_Score)
Y_test <- as.matrix(test_genematrix$Integrated_Score)

normalize_feature <- function(x) {
  x <- gsub("HLA[._-]([A-Za-z0-9]+)", "HLA-\\1", x)
  x <- gsub("_", ".", x)
  x <- gsub("-(?=[A-Za-z0-9]+$)", ".", x, perl = TRUE)
  x
}

#lightGBM
library(lightgbm)
library(dplyr)

train_keep <- complete.cases(X_train) & !is.na(Y_train)
test_keep  <- complete.cases(X_test)  & !is.na(Y_test)

X_train_use <- X_train[train_keep, , drop = FALSE]
y_train_use <- Y_train[train_keep]

X_test_use <- X_test[test_keep, , drop = FALSE]
y_test_use <- Y_test[test_keep]


clean_feature_names <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  make.unique(x, sep = "_")
}

original_features <- colnames(X_train_use)
clean_features <- clean_feature_names(original_features)

colnames(X_train_use) <- clean_features
colnames(X_test_use)  <- clean_features

dtrain_lgb <- lgb.Dataset(
  data = X_train_use,
  label = y_train_use
)

eval_regression <- function(y_true, y_pred) {
  data.frame(
    Pearson_r = cor(y_true, y_pred, method = "pearson"),
    R2 = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2),
    RMSE = sqrt(mean((y_true - y_pred)^2)),
    MAE = mean(abs(y_true - y_pred))
  )
}

final_params_lgb <- list(
  objective = "regression",
  metric = "rmse",
  learning_rate = 0.05,
  num_leaves = 255,
  max_depth = 9,
  min_data_in_leaf = 150,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 0,
  lambda_l1 = 1,
  lambda_l2 = 5,
  min_gain_to_split=0.1,
  verbosity = -1,
  seed = 123
)

set.seed(123)

final_lgb_fit <- lgb.train(
  params = final_params_lgb,
  data = dtrain_lgb,
  nrounds = 1000
)

lgb_importance <- lgb.importance(
  model = final_lgb_fit
)

write.csv(lgb_importance,file = "./lgb_importance.csv")

pred_train_lgb <- predict(final_lgb_fit, X_train_use)
pred_test_lgb  <- predict(final_lgb_fit, X_test_use)

train_metrics_lgb <- eval_regression(y_train_use, pred_train_lgb)
test_metrics_lgb  <- eval_regression(y_test_use, pred_test_lgb)

lgb_metrics <- bind_rows(
  data.frame(Set = "Train", train_metrics_lgb),
  data.frame(Set = "test", test_metrics_lgb)
)

print(lgb_metrics)

plot_fit_scatter <- function(y_true, y_pred, n_plot = 3000, title = "Fit plot") {
  df <- data.frame(True = y_true, Pred = y_pred)
  
  if (nrow(df) > n_plot) {
    idx <- sample(seq_len(nrow(df)), n_plot)
    df_plot <- df[idx, ]
  } else {
    df_plot <- df
  }
  
  r <- cor(y_true, y_pred)
  r2 <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  rmse <- sqrt(mean((y_true - y_pred)^2))
  
  library(ggplot2)
  ggplot(df_plot, aes(x = True, y = Pred)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    labs(
      title = title,
      x = "Observed score",
      y = "Predicted score"
    ) +
    annotate(
      "text",
      x = min(df_plot$True),
      y = max(df_plot$Pred),
      hjust = 0, vjust = 1,
      label = paste0(
        "All cells:\n",
        "r = ", round(r, 3),
        "\nR2 = ", round(r2, 3),
        "\nRMSE = ", round(rmse, 3)
      )
    )
}

plot_fit_scatter(y_train_use, pred_train_lgb, n_plot = 3000, title = "Training set")
plot_fit_scatter(y_test_use, pred_test_lgb, n_plot = 3000, title = "test set")

# randomforest ------
train_df <- data.frame(y = y_train_use, X_train_use)

set.seed(123)

final_rf_fit <- ranger(
  formula = y ~ .,
  data = train_df,
  num.trees = 1000,
  mtry = 75,
  min.node.size = 50,
  sample.fraction = 0.8,
  importance = "permutation",
  num.threads = 6,
  seed = 123
)

pred_train_rf <- predict(
  final_rf_fit,
  data = train_df
)$predictions

pred_test_rf <- predict(
  final_rf_fit,
  data = data.frame(X_test_use)
)$predictions

eval_regression <- function(y_true, y_pred) {
  data.frame(
    Pearson_r = cor(y_true, y_pred, method = "pearson"),
    R2 = 1 - sum((y_true - y_pred)^2) /
      sum((y_true - mean(y_true))^2),
    RMSE = sqrt(mean((y_true - y_pred)^2)),
    MAE = mean(abs(y_true - y_pred))
  )
}

rf_metrics <- bind_rows(
  data.frame(
    Set = "Train",
    eval_regression(Y_train, pred_train_rf)
  ),
  data.frame(
    Set = "test",
    eval_regression(Y_test, pred_test_rf)
  )
)

print(rf_metrics)


rf_importance <- data.frame(
  Feature = normalize_feature(
    names(final_rf_fit$variable.importance)
  ),
  Importance = as.numeric(
    final_rf_fit$variable.importance
  )
) %>%
  arrange(desc(Importance))

print(head(rf_importance, 20))

write.csv(rf_importance,file = "./rf_importance.csv")

plot_fit_scatter <- function(y_true, y_pred, n_plot = 3000, title = "Fit plot") {
  df <- data.frame(True = y_true, Pred = y_pred)
  
  if (nrow(df) > n_plot) {
    idx <- sample(seq_len(nrow(df)), n_plot)
    df_plot <- df[idx, ]
  } else {
    df_plot <- df
  }
  
  r <- cor(y_true, y_pred)
  r2 <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  rmse <- sqrt(mean((y_true - y_pred)^2))
  
  library(ggplot2)
  ggplot(df_plot, aes(x = True, y = Pred)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    labs(
      title = title,
      x = "Observed score",
      y = "Predicted score"
    ) +
    annotate(
      "text",
      x = min(df_plot$True),
      y = max(df_plot$Pred),
      hjust = 0, vjust = 1,
      label = paste0(
        "All cells:\n",
        "r = ", round(r, 3),
        "\nR2 = ", round(r2, 3),
        "\nRMSE = ", round(rmse, 3)
      )
    )
}

plot_fit_scatter(y_train_use, pred_train_rf, n_plot = 3000, title = "Training set")
plot_fit_scatter(y_test_use, pred_test_rf, n_plot = 3000, title = "test set")

#XGBoost
library(xgboost)
library(dplyr)
library(ggplot2)
eval_regression <- function(y_true, y_pred) {
  data.frame(
    Pearson_r = cor(y_true, y_pred, method = "pearson"),
    R2 = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2),
    RMSE = sqrt(mean((y_true - y_pred)^2)),
    MAE = mean(abs(y_true - y_pred))
  )
}

train_keep <- complete.cases(X_train) & !is.na(Y_train)
X_train_use <- as.matrix(X_train[train_keep, , drop = FALSE])
y_train_use <- Y_train[train_keep]

test_keep <- complete.cases(X_test) & !is.na(Y_test)
X_test_use <- as.matrix(X_test[test_keep, , drop = FALSE])
y_test_use <- Y_test[test_keep]

dtrain <- xgb.DMatrix(data = X_train_use, label = y_train_use)
dtest  <- xgb.DMatrix(data = X_test_use, label = y_test_use)

final_params <- list(
  booster = "gbtree",
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.06,
  max_depth = 7,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 1,
  gamma = 0.05,
  lambda = 5,
  alpha = 0,
  seed = 123
)
set.seed(123)

final_xgb_fit <- xgb.train(
  params = final_params,
  data = dtrain,
  nrounds = 100,
  watchlist = list(train = dtrain),
  verbose = 0
)

pred_train_xgb <- predict(final_xgb_fit, newdata = dtrain)
pred_test_xgb  <- predict(final_xgb_fit, newdata = dtest)

train_metrics_xgb <- eval_regression(y_train_use, pred_train_xgb)
test_metrics_xgb  <- eval_regression(y_test_use, pred_test_xgb)

xgb_metrics <- bind_rows(
  data.frame(Set = "Train", train_metrics_xgb),
  data.frame(Set = "test", test_metrics_xgb)
)

print(xgb_metrics)

xgb_importance <- xgb.importance(
  model = final_xgb_fit,
  feature_names = colnames(X_train_use)
)

xgb_importance$Feature <- normalize_feature(xgb_importance$Feature)
write.csv(xgb_importance,file = "./xgb_importance.csv")

plot_fit_scatter <- function(y_true, y_pred, n_plot = 3000, title = "Fit plot") {
  df <- data.frame(True = y_true, Pred = y_pred)
  
  if (nrow(df) > n_plot) {
    idx <- sample(seq_len(nrow(df)), n_plot)
    df_plot <- df[idx, ]
  } else {
    df_plot <- df
  }
  
  r <- cor(y_true, y_pred)
  r2 <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  rmse <- sqrt(mean((y_true - y_pred)^2))
  
  library(ggplot2)
  ggplot(df_plot, aes(x = True, y = Pred)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    labs(
      title = title,
      x = "Observed score",
      y = "Predicted score"
    ) +
    annotate(
      "text",
      x = min(df_plot$True),
      y = max(df_plot$Pred),
      hjust = 0, vjust = 1,
      label = paste0(
        "All cells:\n",
        "r = ", round(r, 3),
        "\nR2 = ", round(r2, 3),
        "\nRMSE = ", round(rmse, 3)
      )
    )
}

plot_fit_scatter(y_train_use, pred_train_xgb, n_plot = 3000, title = "Training set")
plot_fit_scatter(y_test_use, pred_test_xgb, n_plot = 3000, title = "test set")
