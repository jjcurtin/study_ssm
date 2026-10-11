# This script is to process the raw data to avoid repeated preprocessing in chtc running

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
# remove participants with fewer than 8 model days
participant_n <- ave(
  raw_data$subid,
  raw_data$subid,
  FUN = length
)

raw_data <- raw_data[participant_n >= 8, ]

write.csv(
  raw_data,
  "S:/ssm_capstone/data_processed/processed_features_ema_ssm_1x_day_24h.csv",
  row.names = FALSE
)

subid_list <- unique(raw_data$subid)

jobs <- data.frame(
  job_num = seq_along(subid_list),
  subid = subid_list
)

write.table(
  jobs,
  "sherry_repro/jobs.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

