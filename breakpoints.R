# --- Packages ---
library(dplyr)
library(purrr)
library(tibble)
library(segmented)

# --- Helper: one-break segmented fit with AIC grid fallback ---
fit_piecewise <- function(dat, min_trials = 8L, grid_range = 3:23) {
  dat <- dat %>% filter(is.finite(rt), is.finite(trial_index)) %>%
    mutate(trial_index = as.numeric(trial_index)) %>% arrange(trial_index)
  if (nrow(dat) < min_trials || dplyr::n_distinct(dat$trial_index) < 5)
    return(tibble(breakpoint = NA_real_, slope_pre = NA_real_, slope_post = NA_real_,
                  aic = NA_real_, method = "insufficient_data"))
  base_mod <- lm(rt ~ trial_index, data = dat)
  seg_try <- tryCatch(suppressWarnings(
    segmented(base_mod, seg.Z = ~ trial_index,
              psi = stats::median(dat$trial_index, na.rm = TRUE))
  ), error = function(e) NULL)
  if (!is.null(seg_try)) {
    bp <- as.numeric(seg_try$psi[1, "Est."])
    sl <- tryCatch(slope(seg_try)$trial_index, error = function(e) NULL)
    return(tibble(breakpoint = bp,
                  slope_pre  = if (is.null(sl)) NA_real_ else as.numeric(sl[1, "Est."]),
                  slope_post = if (is.null(sl)) NA_real_ else as.numeric(sl[2, "Est."]),
                  aic = AIC(seg_try), method = "segmented"))
  }
  grid <- map_dfr(grid_range, function(bp) {
    fit <- tryCatch(lm(rt ~ trial_index + I(pmax(trial_index - bp, 0)), data = dat),
                    error = function(e) NULL)
    if (is.null(fit)) return(tibble(bp = bp, AIC = Inf, slope_pre = NA_real_, slope_post = NA_real_))
    cf <- coef(fit)
    tibble(bp = bp, AIC = AIC(fit),
           slope_pre = unname(cf[["trial_index"]]),
           slope_post = unname(cf[["trial_index"]] + cf[["I(pmax(trial_index - bp, 0))"]]))
  })
  best <- dplyr::slice_min(grid, AIC, n = 1, with_ties = FALSE)
  tibble(breakpoint = best$bp, slope_pre = best$slope_pre, slope_post = best$slope_post,
         aic = best$AIC, method = "grid_search")
}


# --- Fit per participant on first 25 trials, add to data ---
bp_by_id <- df_with_acc %>%
  filter(trial_index <= 25, is.finite(rt)) %>%
  group_by(participant_id) %>%
  group_modify(~ fit_piecewise(.x)) %>%
  ungroup()

df_with_acc_bp <- df_with_acc %>% left_join(bp_by_id, by = "participant_id")

# --- Plots ---
library(ggplot2)

# Build predictions from the piecewise model using the chosen breakpoint per id
df_pred <- df_with_acc %>%
  filter(trial_index <= 25, is.finite(rt)) %>%
  left_join(bp_by_id, by = "participant_id") %>%
  group_by(participant_id) %>%
  group_modify(~{
    d <- .x
    bp <- d$breakpoint[1]
    if (!is.finite(bp)) return(d %>% mutate(fit = NA_real_))
    m <- lm(rt ~ trial_index + I(pmax(trial_index - bp, 0)), data = d)
    d$fit <- predict(m, newdata = d)
    d
  }) %>%
  ungroup()

# Faceted plot with fitted piecewise lines and breakpoint per participant
p <- ggplot(df_pred, aes(trial_index, rt)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fit), linewidth = 1) +
  geom_vline(data = bp_by_id, aes(xintercept = breakpoint),
             linetype = "dashed") +
  facet_wrap(~ participant_id, scales = "free_y") +
  labs(x = "Trial index", y = "RT", title = "RT by trial with estimated breakpoint") +
  theme_minimal()

print(p)                 # show in viewer
ggsave("rt_piecewise_by_participant.pdf", p, width = 12, height = 8)

# Distribution of breakpoints
p_hist <- ggplot(bp_by_id %>% filter(is.finite(breakpoint)),
                 aes(breakpoint)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
  labs(x = "Estimated breakpoint (trial)", y = "Count",
       title = "Distribution of estimated breakpoints") +
  theme_minimal()

print(p_hist)
ggsave("breakpoint_histogram.pdf", p_hist, width = 7, height = 4)

