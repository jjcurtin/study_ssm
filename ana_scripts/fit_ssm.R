# Setup------

library("tidyverse")

# confirm all libraries below are needed
library("KFAS")
library("MASS",exclude=c('select'))
library("comprehenr")
library("invgamma")
library("data.table")
library("truncnorm")
library("fitdistrplus")
library("mvtnorm")
library('MARSS')
library('jsonlite')

source("scripts/helper_functions.R")

devtools::source_url("https://github.com/jjcurtin/lab_support/blob/main/format_path.R?raw=true")
path_processed <- format_path("risk/data_processed/shared")

out <- data.frame() # for binding results to

# read in data

data <- read_csv(here::here(path_processed, "features_ema_ssm_1x_day_24h.csv"),
                 show_col_types = FALSE)


# Model Setup------

# Initialize list for model
mod <- list()

# Set up dimensions for model (stored in list, stored in mod)
dims <- list()

# Fix observation dimension to 9 EMA + lapse (10 total)
# n:= dimension of observation vector
# g:= first dimension of observation noise coefficient matrix
n <- g <- 10
dims["n"] <- n
dims["g"] <- g

# Fix hidden state dimension at 2 (explore different dimensions in future work)
# m:= dimension of latent state vector
# h:= first dimension of latent state noise coefficient matrix
m <- h <- 2
dims["m"] <- m
dims["h"] <- h 

# Store in mod
mod[["dims"]] <- dims

# Settings------

# Miscellaneous settings
settings <- list()
settings[["tinit"]] = 1

mod[["settings"]] <- settings


# Initial parameters------
mod[["par"]] <- init_par(mod[["dims"]])

# Priors------

# Select the fitting distribution using the id
# not sure where this data file is made or format it is in
# KW: ask Eric about raw prior distribution
prior_path<-'mle_coef_fits.csv'

raw_prior <- read_csv(here::here("../data/", prior_path))

priors <- make_priors(mod, raw_prior, subid, TRUE, prior_subids)
mod[["priors"]] <- priors

# Time series information------
# Check for the length of the time series, add to dims
# Tfinal := final time series index
t_final <- 60 # 60 days on study in simulation data
# First EMA day
first_ema_day <- 1


# Will use all available data up until model is fit
# Model will be fit at 1 week, 2 weeks, 1 month, 2 months, 3 months, 4 months...12 moths
# Make priors on 9/10 held in subids and fit up to 14 individual models for 1/10 held out subids
# temporal test/train splits for individual models (e.g., first model trained on 1 week of data and predictions made for week 2). 

# Update code below to be growing windows of data instead of sliding windows
# training  = all data before model fit/refit timepoint
# test = all data after model fit/refit timepoint and before next model fit/refit timepoint




# Step amounts for data availability windows
width_increment <- 15
# Define the training sequence to be the set of 15-day increments that fit inside the length of this participant's study period
training_seq <- seq(from = 0, to = 45, by = width_increment)
training_seq <- training_seq[training_seq >= first_ema_day & training_seq <= t_final]
# Store final study day for participant
mod$dims['T'] <- t_final


# Test one training width at a time
train_width = training_seq[1]

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



# Print output------
print(out)
write.csv(out,'../outputs/local_test.csv',row.names=FALSE)

