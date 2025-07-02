#!/usr/bin/env Rscript
# Set to TRUE for CHTC runs
run_CHTC <- FALSE
options(error=function()traceback(2))
if(run_CHTC){
  args <- commandArgs(trailingOnly = TRUE)
}

# Relevant libraries
library("KFAS")
library('tidyverse')
library("MASS",exclude=c('select'))
library("comprehenr")
library("invgamma")
library("data.table")
library("truncnorm")
library("fitdistrplus")
library("tictoc")
library("mvtnorm")
library('MARSS')
library('jsonlite')

# Source helper files
source('helper_functions.R')
out<-data.frame()

# Driver function reference
if(run_CHTC){
  process_id <- args[1]
  cluster_id <- args[2]
  ref_key <- args[3]
  type <- args[4]
  file_str <- paste('outputs/',cluster_id,'_',process_id,'.csv',sep='')
} else{
  ref_key <- 'split5_rs10_fold1_id3'
  type <- 'MAP_mle'
}

# Some fixed values (for now, can make them arguments in future)
# UNCOMMENT FOR CHTC
if(run_CHTC){
  datapath <- 'day_labels.csv'
  infopath <- 'subject_info.csv'
  refpath <- 'ssm_fold_reference.json'
} else {
  datapath <- '../data/day_labels.csv'
  infopath <- '../data/subject_info.csv'
  refpath <- '../data/ssm_fold_reference.json'
}

# Process relevant json information
ref_fold_info <- read_json(refpath)
keyval <- ref_key
iter_fold <- ref_fold_info[[keyval]]
subid <- iter_fold[['test']]
prior_subids <- iter_fold[['train']]
cat(paste(keyval,subid))

# Process the data file using the id
raw_data <- read.csv(datapath)
data <- raw_data[raw_data$subid==subid,]

# Process the info file using the id
raw_info <- read.csv(infopath)
info <- raw_info[raw_info$subid==subid,]
last_ema_day <- info$last_morning_ema_day
first_ema_day <- info$first_morning_ema_day
mema_count <- info$mema_count


# # # # # # # # # # # # # # # # # # # # #
#             MODEL SETUP               #
# # # # # # # # # # # # # # # # # # # # #
# Initialize list for model
mod <- list()

# # # # # # # # # # # #
#     Dimensions      #
# # # # # # # # # # # #
# Set up dimensions for model (stored in list, stored in mod)
dims <- list()
# Fix observation dimension to 9 EMA + lapse (10 total)
# n:= dimension of observation vector
# g:= first dimension of observation noise coefficient matrix
n <- g <- 10
dims['n'] <- n
dims['g'] <- g

# Fix hidden state dimension at 2 (explore different dimensions in future work)
# m:= dimension of latent state vector
# h:= first dimension of latent state noise coefficient matrix
m <- h <- 2
dims['m'] <- m
dims['h'] <- h 

# Store in mod
mod[['dims']] <- dims

# # # # # # # # # # # #
#     Settings        #
# # # # # # # # # # # #
# Miscellaneous settings
settings<-list()
settings[['tinit']]=1

mod[['settings']]<-settings


# # # # # # # # # # # #
# Initial Parameters  #
# # # # # # # # # # # #
mod[['par']] <- init_par(mod[['dims']])

# # # # # # # # # # # #
#        Priors       #
# # # # # # # # # # # #
# Select the fitting distribution using the id
prior_path<-'mle_coef_fits_lapse_not_fitted.csv'

if(run_CHTC){
  raw_prior <- read.csv(prior_path)
} else {
  raw_prior <- read.csv(paste('../data/',prior_path,sep=''))
}

priors <- make_priors(mod, raw_prior, subid, TRUE,prior_subids)
mod[['priors']] <- priors

# # # # # # # # # # # #
#     Data Horizon    #
# # # # # # # # # # # #
# Check for the length of the time series, add to dims
# Tfinal := final time series index
Tfinal <- min(last_ema_day,89) + 1
# First EMA day using 1-indexing
first_ema_day_1i <- first_ema_day + 1
# Step amounts for data availability windows
width_increment <- 15
# Define the training sequence to be the set of 15-day increments that fit inside the length of this participant's study period
training_seq <- seq(from=0,to=75,by=width_increment)
training_seq <- training_seq[training_seq>=first_ema_day_1i & training_seq<=Tfinal]
# Store final study day for participant
mod$dims['T'] <- Tfinal

tic("Runtime")
# Loop over the training widths
for (train_width in training_seq){
  # Define each possible start/end day based on the training width (that fits into the study period)
  for(start_day in seq(from=first_ema_day_1i,to=(Tfinal-train_width+1),by=1)){
    end_day <- start_day+train_width-1
    TT <- train_width
    relative_Tfinal <- Tfinal-start_day+1
    print(paste("start:",start_day,", end:",end_day))
    mod$dims['TT']<- TT
    
    # # # # # # # # # # # #
    #     Attach Data     #
    # # # # # # # # # # # #
    # Put in the data
    feat_rows = c(3,4,5,6,7,8,9,10,11,14)
    full_raw_datamat <- t(as.matrix(data))[feat_rows,start_day:Tfinal]
    raw_datamat <- masked_datamat <- t(as.matrix(data))[feat_rows,start_day:end_day]
    # Mask the final lapse state
    masked_datamat[10,TT]<-NA
    mod[['data']]<-masked_datamat
    act_0 <- raw_datamat[10,TT]
    window_length <- 7
    
    # # # # # # # # # # # # # # # # # # # # #
    #             MODEL FITTING             #
    # # # # # # # # # # # # # # # # # # # # #
    if(type=='MAP_mle'){
      zero_priors_bool <- FALSE
    }else {
      zero_priors_bool <- TRUE
    }
    
    # Zero out the priors for the MLE fit
    if(zero_priors_bool){
      zero_list <- c("A","B","c","d","R")
      local_mod <- zero_priors(mod,zero_list)
    } else {
      # No zeroing of priors, just leave as is
      local_mod <- mod
    }
    
    # # # # # # # # # # # #
    #     Fit/Predict     #
    # # # # # # # # # # # #
    # Fit model
    local_mod <- run_kf(local_mod)
    iters <- 2500
    conv_tol <- 0.0001
    fit_obj <- run_em(local_mod,iters,conv_tol,zero_priors_bool)
    local_mod <- fit_obj[['model']]
    
    # Prepare coefficient dataframe
    need_flip <- flip_check(local_mod)
    coef_out <- make_coef_frame(local_mod,need_flip,FALSE)
    
    # Make 0-7 preds
    # Same day prediction
    pred_0 <- local_mod$kf_proc$ytT[10,TT]
    # # # # # # # # # # # #
    #     Iter Output     #
    # # # # # # # # # # # #
    # File outputs
    out_frame_mod_obj <- fit_obj
    mod_pars <- local_mod$par
    meanTT <- local_mod$kf_proc$xtT[,TT]
    varTT <- local_mod$kf_proc$VtT[,,TT]

    # Label names
    prob_col_names <- c("w0_pred","w1_pred","w2_pred","w3_pred","w4_pred","w5_pred","w6_pred","w7_pred")
    act_col_names <- c("w0_act","w1_act","w2_act","w3_act","w4_act","w5_act","w6_act","w7_act")
    
    # Get the actual lapse outcomes for the different windows
    actuals_mat <- matrix(NA, nrow = 1, ncol = (window_length+1))
    actuals_df <- data.frame(actuals_mat)
    for(j in 0:(window_length)){
      if((TT+j)<=relative_Tfinal){
        relevant_lapse_data <- full_raw_datamat[10,TT:(TT+j)]
        lapse_bool <- as.numeric(1 %in% relevant_lapse_data)
        val <- lapse_bool
      } else {
        val <- NA
      }
      actuals_df[(j+1)]<- val
    }
    # Make window predictions
    raw_window_preds <- get_window_preds(mod_pars,masked_datamat,meanTT,varTT,type,TT,relative_Tfinal,window_length)
    proc_window_preds <- vapply(raw_window_preds,mean,numeric(1))
    probs_df <- data.frame(t(proc_window_preds))
    # Adjust column names
    colnames(probs_df)<-prob_col_names
    colnames(actuals_df)<-act_col_names
   
    pred_act_df <- cbind(probs_df,actuals_df)
    fit_out <- get_rolling_out_frame(type,subid,TT,end_day,out_frame_mod_obj,pred_0,act_0)
    full_out <- cbind(fit_out,pred_act_df)
    full_out$key <- keyval
    out <- rbind(out, full_out)
  }
}
toc()
# # # # # # # # # # # # # # # # # # # # #
#             CHTC OUTPUT               #
# # # # # # # # # # # # # # # # # # # # #
# UNCOMMENT FOR CHTC
if(run_CHTC){
  write.csv(out,file_str,row.names=FALSE)
} else {
  print(out)
  write.csv(out,'../outputs/local_test.csv',row.names=FALSE)
}
