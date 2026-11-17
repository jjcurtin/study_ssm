library(tidyverse)
library(tidymodels)

# Load prediction results
pred_out <- read_csv("S:/ssm_capstone/chtc/batch1/output/batch_results.csv")

# Assign each prediction to a 30-day window
pred_out <- pred_out |>
  mutate(
    window = (day - 8) %/% 30 + 1
  )

# Count predictions for each participant in each window
participant_window <- pred_out |>
  group_by(subid, window) |>
  summarise(
    n_days = n_distinct(day),
    .groups = "drop"
  )

# Summarize completeness by window
window_completeness <- participant_window |>
  mutate(
    expected_days = case_when(
      window %in% c(1, 2) ~ 30,
      window == 3 ~ 22
    )
  ) |>
  group_by(window) |>
  summarise(
    n_participants = n(),
    n_complete = sum(n_days == expected_days),
    n_incomplete = sum(n_days < expected_days),
    .groups = "drop"
  )

window_completeness

# Identify participants with complete first two windows
fixed_cohort <- participant_window |>
  filter(window %in% c(1, 2), n_days == 30) |>
  group_by(subid) |>
  summarise(
    n_complete_windows = n_distinct(window),
    .groups = "drop"
  ) |>
  filter(n_complete_windows == 2)

nrow(fixed_cohort)
# Check the fixed cohort
pred_fixed <- pred_out |>
  semi_join(fixed_cohort, by = "subid") |>
  filter(window %in% c(1, 2))

# Prepare predictions for evaluation
pred_fixed_eval <- pred_fixed |>
  filter(!is.na(actual), is.finite(pred)) |>
  mutate(
    truth = factor(
      actual,
      levels = c(0, 1),
      labels = c("no_lapse", "lapse")
    ),
    prob = pmax(0, pmin(1, pred)),
    prob_safe = pmax(1e-15, pmin(1 - 1e-15, prob))
  )
# Check observations and lapse distribution
pred_fixed_eval |>
  group_by(window) |>
  summarise(
    n_predictions = n(),
    n_lapse = sum(actual == 1),
    lapse_rate = mean(actual == 1),
    .groups = "drop"
  )
# AUROC by window
auc_fixed <- pred_fixed_eval |>
  group_by(window) |>
  roc_auc(
    truth = truth,
    prob,
    event_level = "second"
  )

auc_fixed
# Log Loss by window
logloss_fixed <- pred_fixed_eval |>
  group_by(window) |>
  mn_log_loss(
    truth = truth,
    prob_safe,
    event_level = "second"
  )

logloss_fixed

pred_fixed_eval |>
  group_by(window) |>
  summarise(
    mean_pred = mean(prob),
    mean_actual = mean(actual),
    min_pred = min(pred),
    max_pred = max(pred),
    .groups = "drop"
  )

pred_eval <- pred_out |>
  filter(!is.na(actual), is.finite(pred)) |>
  mutate(
    window = (day - 8) %/% 30 + 1,
    
    truth = factor(
      actual,
      levels = c(0, 1),
      labels = c("no_lapse", "lapse")
    ),
    
    # Clip predictions to [0, 1]
    prob = pmax(0, pmin(1, pred)),
    
    # Avoid log(0) when calculating Log Loss
    prob_safe = pmax(1e-15, pmin(1 - 1e-15, prob))
  )
# Overall AUROC
overall_auc <- pred_eval |>
  roc_auc(
    truth = truth,
    prob,
    event_level = "second"
  )

# Overall Log Loss
overall_logloss <- pred_eval |>
  mn_log_loss(
    truth = truth,
    prob_safe,
    event_level = "second"
  )

overall_auc
overall_logloss
# AUROC by window
auc_window <- pred_eval |>
  group_by(window) |>
  roc_auc(
    truth = truth,
    prob,
    event_level = "second"
  )

# Log Loss by window
logloss_window <- pred_eval |>
  group_by(window) |>
  mn_log_loss(
    truth = truth,
    prob_safe,
    event_level = "second"
  )
# Combine metrics
window_results <- pred_eval |>
  group_by(window) |>
  summarise(
    n_participants = n_distinct(subid),
    n_predictions = n(),
    lapse_rate = mean(actual == 1),
    .groups = "drop"
  ) |>
  left_join(
    auc_window |> select(window, AUROC = .estimate),
    by = "window"
  ) |>
  left_join(
    logloss_window |> select(window, LogLoss = .estimate),
    by = "window"
  )

window_results
