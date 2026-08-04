# this script came from eric in eric/all/coef_fitting folder

# Complete data
library('tidyverse')
library('tictoc')
library('foreach')
library('KFAS')
library("MASS",exclude=c('select'))
library('doParallel')
library("data.table")
source('mle_coef_fit.R')
source('/Users/eric/repos/aud/aud_helper_functions.R')
datapath <- '/Users/eric/repos/aud/data/day_labels.csv'
infopath <- '/Users/eric/repos/aud/data/subject_info.csv'
options(error=function()traceback(2))
data <- read.csv(datapath)
info <- read.csv(infopath)

# subid_list <- unique(data$subid)
# Remove certain subids
subid_list <- info$subid[info$responsiveness>=0.25]

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
parout<-foreach(idx=1:length(subid_list)) %dopar% {
#parout<-foreach(idx=1:3) %dopar% {
#for(idx in 105:115){
  test_subid<-subid_list[idx]
  cat(paste("Subject:",test_subid,"\n"))
  subid_data <- data[data$subid==test_subid,]
  subid_info <- info[info$subid==test_subid,]
  fit_info <- run_mle_fit(subid_data,subid_info,FALSE)
  return(fit_info)
}
toc()
coef_out <- rbindlist(lapply(parout,function(x){x[[2]]}))
writepath <- '/Users/eric/repos/aud/data/mle_coef_fits_lapse_not_fitted.csv'
write.csv(coef_out,file=writepath)

