# This script is to process the raw data to avoid repeated preprocessing in chtc running

# Relevant libraries
library('tidyverse')

#--------------------------------------------------
# Read raw data
#--------------------------------------------------
raw_data <- read.csv("S:/ssm_capstone/data_raw/features_ema_ssm_1x_day_24h.csv")

#remove duplicates
raw_data <- raw_data[
  !duplicated(raw_data[, c("subid", "dttm_label")]),
]

#recode the lapse so that it would be numeric other that chr
raw_data$lapse <- ifelse(raw_data$lapse == "lapse", 1, 0)

# convert timestamp
raw_data$dttm_label <- as.POSIXct(
  raw_data$dttm_label,
  format = "%Y-%m-%dT%H:%M:%SZ",
  tz = "UTC"
)
# create relative day within each participant
raw_data$day <- ave(
  raw_data$subid,
  raw_data$subid,
  FUN = seq_along
)

write.csv(
  raw_data,
  "S:/ssm_capstone/data_processed/processed_features_ema_ssm_1x_day_24h.csv",
  row.names = FALSE
)


