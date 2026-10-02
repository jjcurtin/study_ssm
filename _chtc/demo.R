# Relevant libraries
#library("dplyr")
#library("readr")
#library("stringr")

# Source helper files
source('helper_functions.R') 
source('mle_coef_fit.R') 


# get args
args <- commandArgs(trailingOnly = TRUE) 
job_num_arg <- args[1]
subid_arg <- args[2]

# ALL COMPUTATIONS HERE

tibble::tibble(subid = subid_arg) |>
  readr::write_csv(stringr::str_c("results_", job_num_arg, ".csv"))