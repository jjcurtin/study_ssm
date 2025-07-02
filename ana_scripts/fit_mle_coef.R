# fit full models on each subid to be used for priors

library("tidyverse")
library("foreach")
library("KFAS")
library("MASS",exclude=c('select'))
library("doParallel")
library("data.table")

source("_ssm/mle_coef_fit.R")
source("_ssm/helper_functions.R")
source("https://github.com/jjcurtin/lab_support/blob/main/format_path.R?raw=true")

path_data <- format_path("risk/data_processed/ssm")

data <- read_csv(here::here(path_data, "features_ema_ssm_1x_day_24h.csv"),
                 show_col_types = FALSE)
info <- read_csv(here::here(path_data, "subid_info.csv"),
                 show_col_types = FALSE)

subids <- info |> 
  pull(subid)

coef_out <- subids |> 
  map(\(subid) run_mle_fit(subid, data, info)) |> 
  list_rbind()

coef_out |> 
  write_csv(here::here(path_data, "mle_coef_fits.csv"))
