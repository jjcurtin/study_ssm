# Script to make subid info for model fitting scripts

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(foreach)
  library(tidyr) # for pivot_wider.
}))

devtools::source_url("https://github.com/jjcurtin/lab_support/blob/main/format_path.R?raw=true")
path_processed <- format_path("risk/data_processed/shared")
path_ssm <- format_path("risk/data_processed/ssm")


# read in EMA (morning only)
ema <- read_csv(here::here(path_processed, "ema_morning.csv"), 
                show_col_types = FALSE) |> 
  mutate(subid = as.numeric(subid),
         dttm_obs = with_tz(start_date, tz = "America/Chicago")) |> 
  select(subid, dttm_obs) 

# read in study dates
study_dates <- read_csv(here::here(path_processed, "study_dates_ema.csv"), 
                        show_col_types = FALSE) 


count_ema <- ema |> 
  group_by(subid) |> 
  count()

first_ema <- ema |> 
  group_by(subid) |> 
  arrange(dttm_obs) |> 
  slice_head(n = 1) |> 
  rename(first_ema = dttm_obs)

last_ema <- ema |> 
  # filter out two EMAs taken in 2020
  filter(dttm_obs < as_datetime("2020-01-01 00:00:00")) |> 
  group_by(subid) |> 
  arrange(dttm_obs) |> 
  slice_tail(n = 1) |> 
  rename(last_ema = dttm_obs)

subid_info <- count_ema |> 
  full_join(first_ema, by = "subid") |> 
  full_join(last_ema, by = "subid") |> 
  filter(subid %in% study_dates$subid)

# calculate adherence
subid_info <- subid_info |> 
  mutate(study_days = as.numeric(date(last_ema) - date(first_ema)),
         adherence = n/study_days) 



subid_info |> 
  write_csv(here::here(path_ssm, "subid_info.csv"))
