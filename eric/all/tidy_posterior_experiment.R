library(tidymodels)
library(tidyposterior)

data(two_class_dat)

ggplot(two_class_dat, aes(x = A, y = B, col = Class)) + 
  geom_point(alpha = 0.3, cex = 2) +
  coord_fixed()

set.seed(1)
cv_folds <- vfold_cv(two_class_dat)
cv_folds

logistic_spec <- logistic_reg() %>% set_engine("glm")

spline_rec <- 
  recipe(Class ~ ., data = two_class_dat) %>% 
  step_ns(A, B, deg_free = 3)

spline_wflow <- 
  workflow() %>% 
  add_recipe(spline_rec) %>% 
  add_model(logistic_spec)

spatial_sign_rec <- 
  recipe(Class ~ ., data = two_class_dat) %>% 
  step_normalize(A, B) %>% 
  step_spatialsign(A, B)

spatial_sign_rec %>% 
  prep() %>% 
  bake(new_data = NULL) %>% 
  ggplot(aes(x = A, y = B, col = Class)) + 
  geom_point(alpha = 0.3, cex = 2) +
  coord_fixed()

spatial_sign_wflow <- 
  workflow() %>% 
  add_recipe(spatial_sign_rec) %>% 
  add_model(logistic_spec)

compute_roc <- function(split, wflow) {
  # Fit the model to 90% of the data
  mod <- fit(wflow, data = analysis(split))
  # Predict the other 10%
  pred <- predict(mod, new_data = assessment(split), type = "prob")
  # Compute the area under the ROC curve
  pred %>% 
    bind_cols(assessment(split)) %>% 
    roc_auc(Class, .pred_Class1) %>% 
    pull(.estimate)
}

roc_values <- 
  cv_folds %>% 
  mutate(
    spatial_sign = map_dbl(splits, compute_roc, spatial_sign_wflow),
    splines      = map_dbl(splits, compute_roc, spline_wflow)
  )

roc_values

# Overall ROC statistics per workflow:
summarize(
  roc_values,
  splines = mean(splines),
  spatial_sign = mean(spatial_sign)
)

rset_mod <- perf_mod(roc_values, seed = 2, iter = 5000,hetero_var = TRUE, chains = 5, refresh = 1)

print(summary(rset_mod), digits = 3)

tidy(rset_mod, seed = 3)

autoplot(rset_mod)

preproc_diff <- contrast_models(rset_mod, seed = 4)
summary(preproc_diff, seed = 5)

autoplot(preproc_diff) + 
  xlab("Difference in ROC (spatial sign - splines)")

summary(preproc_diff, size = 0.02) %>% 
  select(contrast, starts_with("pract"))

autoplot(preproc_diff, size = 0.02)

data(precise_example)
