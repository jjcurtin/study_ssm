library('tidyverse')

infopath <- '/Users/eric/repos/aud/data/subject_info.csv'
#datapath <- 'data/noisy_marss_coef_fits_lapse_not_fitted.csv'
outpath <- '/Users/eric/repos/aud/chtc/subid_list.txt'

info <- read.csv(infopath)
#data <- read.csv(datapath)

subid_count <- dim(info)[1]
print(subid_count)
first_write <- TRUE

for(i in 1:subid_count){
  info_entry <- info[i,]
  data_entry <- data[i,]
  id <- info_entry$subid
  write_bool <- FALSE
  # Conditions for inclusion
  if(info_entry$last_morning_ema_day>= 89 & info_entry$mema_count>= 50){
    write_bool <- TRUE
  }
  if(write_bool){
    if(first_write){
      cat(paste(id,sep=''),file=outpath,append=FALSE)
      first_write <- FALSE
    } else {
      cat(paste('\n',id,sep=''),file=outpath,append=TRUE)
    }
  }
}