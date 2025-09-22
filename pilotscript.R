library(readxl)
library(dplyr)
library(rstatix)
library(purrr)
library(rpart)
library(rpart.plot)
library(lme4)
library(restriktor)

df <- read_excel("sabolona_baserate (1).xlsx") %>%
  mutate(
    participant_id = as.character(id),
    Response = as.numeric(Response),
    response = as.factor(Response),
    type = as.factor(type),
    ratio = as.numeric(ratio),
    rt = as.numeric(rt),
    item_stem = as.character(item_stem)
  )  %>% group_by(participant_id) %>%
  arrange(participant_id) %>%
  mutate(trial_index = row_number()) %>%
  ungroup()


# Study 1 -----------------------------------------------------------------

## H1.1 -----------------------------------------------------------------

# Aggregate accuracy across type per participant
df_acc <- df %>%
  group_by(participant_id, type) %>%
  summarise(Response = mean(Response, na.rm = TRUE), .groups = "drop")

# Run ANOVA
anova_res <- anova_test(
  data = df_acc,
  dv = Response,
  wid = participant_id,
  within = type
)

anova_res

# Pairwise t-tests with Bonferroni correction
df_acc %>%
  pairwise_t_test(
    Response ~ type,
    paired = TRUE,         # within-subjects
    p.adjust.method = "bonferroni"
  )

df_acc %>%
  cohens_d(Response ~ type, paired = TRUE)

# item_level

glmm_acc <- glmer(Response ~ -1 + type + (1| participant_id) + (1 | item_stem), 
                  df %>%
                    mutate(type = recode(type,
                                         "konfliktní" = "incongruent",
                                         "nekonfliktní" = "congruent")) %>%
                    mutate(type = factor(type, levels = c("congruent", "incongruent"))),
      family = "binomial")

h1.1 <- "typeincongruent < -.85"

goric(glmm_acc, hypotheses = list(h1 = h1.1))
## H1.2+H1.3 -----------------------------------------------------------------

df_rt <- df %>%
  group_by(participant_id, type) %>%
  summarise(rt = log(mean(rt, na.rm = TRUE)), .groups = "drop")

# Run ANOVA
 anova_test(
  data = df_rt,
  dv = rt,
  wid = participant_id,
  within = type
)
 
 # Recode into 3-level within-subject factor
 df3 <- df %>%
   mutate(
     condition = case_when(
       type == "nekonfliktní" ~ "congruent",
       type == "konfliktní" & Response == 0 ~ "incon_resp0",
       type == "konfliktní" & Response == 1 ~ "incon_resp1",
       TRUE ~ NA_character_
     )
   ) %>%
   filter(!is.na(condition))

 df3_agg <- df3 %>%
   group_by(participant_id, condition = factor(
     condition,
     levels = c("congruent", "incon_resp0", "incon_resp1")
   )) %>%
   summarise(rt = log(mean(rt, na.rm = TRUE)), .groups = "drop")
 
 # congruent vs incongruent_resp0
 p10 <- df3_agg %>%
   tidyr::pivot_wider(names_from = condition, values_from = rt)
 
 t.test(p10$congruent, p10$incon_resp0, paired = TRUE)
 
 # congruent vs incongruent_resp1
 p11 <- df3_agg %>%
   tidyr::pivot_wider(names_from = condition, values_from = rt)
 
 t.test(p11$congruent, p11$incon_resp1, paired = TRUE)
 
 # item_level
 
 glmm_rt <- lmer(log(rt) ~ -1 + condition + (1| participant_id) + (1 | item_stem), 
                 df3 %>% filter(rt != 0) )
 h1.2 <- "conditionincon_resp0 > conditioncongruent; conditionincon_resp1 > conditioncongruent"
 
 goric(glmm_rt, hypotheses = list(h1 = h1.2))
 
# Study 2 --------------------------------------------------------------------

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
       control = rpart.control(minsplit = 2, cp = 0.001, xval = 0, maxdepth = 3)
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
 

## H2.1 --------------------------------------------------------------------

 # Predicting accuracy with age (NfCC was not collected in the pilot)
lm(participant_accuracy ~ V_k, df_with_acc) 
 
## H2.2 --------------------------------------------------------------------
 
 # Predicting initial RT with age (NfCC was not collected in the pilot)

lm(rt ~ V_k, df_with_acc %>% filter(trial_index <= 10)) 
 
## H2.3 --------------------------------------------------------------------
 
 participant_trees <- df %>%
   split(.$participant_id) %>%
   map(~ rpart(response ~ type + ratio, data = ., method = "class"))
 
 # 2) Extract all ratio cutpoints (+ split improvement, if available)
 ratio_rules_raw <- map_dfr(
   participant_trees,
   function(fit) {
     sp <- fit$splits
     if (is.null(sp)) return(tibble(ratio_cut = numeric(0), improve = numeric(0)))
     idx <- rownames(sp) == "ratio"
     if (!any(idx)) return(tibble(ratio_cut = numeric(0), improve = numeric(0)))
     tibble(
       ratio_cut = as.numeric(sp[idx, "index"]),
       improve   = suppressWarnings(as.numeric(sp[idx, "improve"]))
     )
   },
   .id = "participant_id"
 )
 
 # 3) Collapse to participant-level averages
 ratio_rules_avg <- ratio_rules_raw %>%
   group_by(participant_id) %>%
   summarise(
     n_ratio_splits   = n(),                                   # how many ratio splits they have
     ratio_cut_mean   = if (n_ratio_splits > 0) mean(unique(ratio_cut)) else NA_real_,
     ratio_cut_wmean  = if (n_ratio_splits > 0 && sum(pmax(improve, 0), na.rm = TRUE) > 0)
       stats::weighted.mean(ratio_cut, w = pmax(improve, 0), na.rm = TRUE)
     else ratio_cut_mean,                   # fallback to unweighted
     .groups = "drop"
   )
 
 # 4) Append to your participant-level table
 df_with_acc <- df_with_acc %>%
   left_join(ratio_rules_avg, by = "participant_id")

 # Predicting split rule with age (NfCC was not collected in the pilot)
 lm(ratio_cut ~ V_k, df_with_acc) 
 

  
 
 