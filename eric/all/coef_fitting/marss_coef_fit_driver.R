# Complete data
library('tidyverse')
library('tictoc')
library('foreach')
library('doParallel')
source('marss_coef_fit.R')
#datapath <- '/Users/eric/repos/aud/data/day_labels.csv'
datapath <- '/Users/eric/repos/aud/data/day_labels.csv'
infopath <- '/Users/eric/repos/aud/data/subject_info.csv'

data <- read.csv(datapath)
info <- read.csv(infopath)

subid_list <- unique(data$subid)
# Remove certain subids?
subid_list <- subid_list[subid_list != 193 & subid_list != 232]

# Parallel setup
n.cores <- detectCores() - 1
my.cluster <- makeCluster(
  n.cores,
  type="FORK"
)

doParallel::registerDoParallel(cl=my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

# Uncomment for parallel runs
tic('runtime')
#parout<-foreach(idx=1:length(subid_list)) %dopar% {
#parout<-foreach(idx=1:3) %dopar% {
for(test_subid in c(1)){
  print(idx)
  test_subid<-subid_list[idx]
  subid_data <- data[data$subid==test_subid,]
  subid_info <- info[info$subid==test_subid,]
  fit_info <- run_2ni_fit(subid_data,subid_info,FALSE)
  #print(subid_info$mema_count)
  #print(subid_info$last_morning_ema_day)
  #print(fit_info$fitting_error)
  return(fit_info)
}
#coef_out <- rbindlist(lapply(parout,function(x){x[[2]]}))
#writepath <- 'data/marss_long_coef_fits_lapse_not_fitted.csv'
#write.csv(coef_out,file=writepath)
toc()

# # Uncomment for single subject debugging
# source('marss_coef_fit.R')
# test_subid <- 23
# subid_data <- data[data$subid==test_subid,]
# subid_info <- info[info$subid==test_subid,]
# fit_info <- run_2ni_fit(subid_data,subid_info)
# #print(test_subid)
# #print(fit_info$fitting_error)