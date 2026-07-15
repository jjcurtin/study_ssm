library(tidymodels)
library(tidyposterior)
library(data.table)

auc_vals <- read.csv("synthetic_performance_dataframe.csv")

mod <- perf_mod(auc_vals, seed = 1,iter = 5000,
                hetero_var = FALSE, transform = tidyposterior::logit_trans,
                chains = 5, adapt_delta=0.99, refresh = 1)

print(summary(mod), digits = 3)
