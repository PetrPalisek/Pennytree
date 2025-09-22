# ── Packages ────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(rpart)
  library(rpart.plot)
})

df <- readxl::read_excel("sabolona_baserate (1).xlsx") %>%
  mutate(
    participant_id = as.character(id),
    response = as.factor(Response),
    type = as.factor(type),
    ratio = as.numeric(ratio),
    rt = as.numeric(rt),
    item_stem = as.character(item_stem)
    
  )

# ── Helper: k-fold accuracy for one participant ────────────────────────────
participant_tree_accuracy <- function(dat, seed = 123) {
  dat <- droplevels(dat)
  
  # Need at least 2 classes to train a classifier
  if (length(unique(dat$response)) < 2L) return(NA_real_)
  
  n <- nrow(dat)
  # Use 5-fold CV when possible, otherwise fall back to LOOCV
  k <- if (n >= 5) 5L else n
  
  set.seed(seed)
  folds <- sample(rep(seq_len(k), length.out = n))
  
  fold_acc <- numeric(k)
  
  for (f in seq_len(k)) {
    train <- dat[folds != f, , drop = FALSE]
    test  <- dat[folds == f, , drop = FALSE]
    
    # rpart prefers at least a few rows; guard tiny splits
    if (nrow(train) < 2 || length(unique(train$response)) < 2) {
      fold_acc[f] <- NA_real_
      next
    }
    
    fit <- rpart(
      response ~ type + ratio,
      data   = train,
      method = "class",
      control = rpart.control(minsplit = 2, cp = 0.001, xval = 0)
    )
    
    preds <- predict(fit, newdata = test, type = "class")
    fold_acc[f] <- mean(preds == test$response)
  }
  
  # Mean accuracy across folds (ignoring NA folds)
  if (all(is.na(fold_acc))) NA_real_ else mean(fold_acc, na.rm = TRUE)
}

# ── 1) Per-participant accuracies ───────────────────────────────────────────
acc_by_participant <- df %>%
  group_by(participant_id) %>%
  group_split() %>%
  map_dfr(function(d) {
    tibble(
      participant_id = d$participant_id[1],
      participant_accuracy = participant_tree_accuracy(d)
    )
  })

# ── 2) Append accuracy back to the long data ───────────────────────────────
df_with_acc <- df %>%
  left_join(acc_by_participant, by = "participant_id")


hist(df_with_acc$participant_accuracy)

# ── 3) Global (pooled) decision tree & plot ────────────────────────────────
# Predict response from type and ratio across all participants
tree_all <- rpart(
  rt ~  type + ratio ,
  data   = df,
  method = "anova",
  control = rpart.control(minsplit = 10, cp = 0.001, xval = 10)
)

fit <- fit_msdt(
  data    = df,
  N       = 3,
  formula = response ~ type + ratio,
  id_col  = "participant_id",
  minbucket = 50, cp = 0.001, max_iter = 500, conv_tol = 1e-3
  
  
)


state_hat <- max.col(fit$state_probs)           # 1..N
df$hmm_state <- factor(state_hat, levels = 1:fit$N)

fit$gamma %>% round(2)

# Inspect trees per state
fit$mod[[1]] %>% rpart.plot()
fit$mod[[2]] %>% rpart.plot()
fit$mod[[3]] %>% rpart.plot()
# ── Accuracy bins (10% width) ──────────────────────────────────────────────
acc_bins <- acc_by_participant %>%
  filter(!is.na(participant_accuracy)) %>%
  mutate(
    acc_pct = 100 * participant_accuracy,
    bin = cut(
      acc_pct,
      breaks = seq(0, 100, by = 10),
      include.lowest = TRUE,
      right = TRUE,
      labels = c("0–10%", "10–20%", "20–30%", "30–40%", "40–50%",
                 "50–60%", "60–70%", "70–80%", "80–90%", "90–100%")
    )
  )

# Sample 1 participant per (non-empty) bin
set.seed(123)
sampled_ids <- acc_bins %>%
  group_by(bin) %>%
  slice_sample(n = 1) %>%            # takes 1 if available
  ungroup() %>%
  select(bin, participant_id, participant_accuracy)

# ── Plot per-participant trees (up to 10) ───────────────────────────────────
# Layout: 2 rows × 5 columns (adjust if fewer bins are populated)
n_panels <- nrow(sampled_ids)
par(mfrow = c(2, 5), mar = c(2.5, 2, 3, 1))  # compact margins

invisible(purrr::pmap(
  sampled_ids,
  function(bin, participant_id, participant_accuracy) {
    d <- df %>% filter(participant_id == !!participant_id) %>% droplevels()
    # Guard tiny/degenerate cases
    if (nrow(d) < 2L || length(unique(d$response)) < 2L) {
      plot.new(); title(main = paste(bin, "\n", participant_id, "\n(no model)")); return(invisible())
    }
    fit <- rpart(
      response ~ type + ratio,
      data = d, method = "class",
      control = rpart.control(minsplit = 2, cp = 0.001, xval = 0)
    )
    rpart.plot(
      fit,
      main = sprintf("%s | %s\nAcc: %.1f%%", bin, participant_id, 100 * participant_accuracy),
      extra = 104, under = TRUE, faclen = 0, varlen = 0
    )
    invisible()
  }
))

# (Optional) Reset plotting layout
# par(mfrow = c(1, 1))

