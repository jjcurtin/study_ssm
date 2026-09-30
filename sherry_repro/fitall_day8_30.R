# This script is to fit SSM on all patient's data
# and predict day8 - end
# Relevant libraries
library("KFAS")
library('tidyverse')
library("MASS",exclude=c('select'))
library("data.table")
library("tictoc")
library('MARSS')

# Source helper files
source('sherry_repro/helper_functions.R') # from eric's 
source('sherry_repro/mle_fit/mle_coef_fit.R') # from eric's

#--------------------------------------------------
# Read participants
#--------------------------------------------------
raw_data <- read.csv("S:/ssm_capstone/data_processed/processed_features_ema_ssm_1x_day_24h.csv")
subid_info <- read.csv("S:/ssm_capstone/data_raw/subid_info.csv")

subid_list <- subid_info$subid

#--------------------------------------------------
# Define Matrix
#--------------------------------------------------

# # # # # # # # # # # # # # # # # # # # #
#             MODEL SETUP               #
# # # # # # # # # # # # # # # # # # # # #
# Initialize list for model
mod <- list()

# # # # # # # # # # # #
#     Dimensions      #
# # # # # # # # # # # #

# Set up dimensions for model (stored in mod-list)
dims <- list()

#Y_t = aX_t + c + v_t (observation equation)
# Fix observation dimension to 9 EMA + lapse (10 total)
# n:= dimension of observation vector
# g:= dimension of observation noise coefficient matrix v
n <- g <- 10
dims['n'] <- n
dims['g'] <- g

#X_t+1 = bX_t + d + w_t (transition equation)
# Fix hidden state dimension at 2
# m:= dimension of latent state vector
# h:= dimension of latent state noise coefficient matrix w
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

# create Observation Matrix
obs_cols <- c(
  "ema_2.p0.l0.rrecent_response",
  "ema_3.p0.l0.rrecent_response",
  "ema_4.p0.l0.rrecent_response",
  "ema_5.p0.l0.rrecent_response",
  "ema_6.p0.l0.rrecent_response",
  "ema_7.p0.l0.rrecent_response",
  "ema_8.p0.l0.rrecent_response",
  "ema_9.p0.l0.rrecent_response",
  "ema_10.p0.l0.rrecent_response",
  "lapse"
)

# # # # # # # # # # # # # # # # # # # # #
#             MODEL FITTING             #
# # # # # # # # # # # # # # # # # # # # #

# ZERO PRIORS (placeholder so run_em() works)
priors <- make_zero_priors(mod)
mod[['priors']] <- priors
zero_priors_bool <- TRUE

# Store predictions
pred_out <- data.frame(
  subid = integer(),
  day = integer(),
  pred = numeric(),
  actual = numeric()
)

base_mod <- mod

for (subid in subid_list) {
  
  # Current participant
  data_i <- raw_data[raw_data$subid == subid, ]
  data_i <- data_i[order(data_i$dttm_label), ]
  
  # Observation matrix: 10 x T
  Y <- t(as.matrix(data_i[, obs_cols]))
  Tfinal <- ncol(Y)
  
  # Fewer than 8 model days -> skip
  if (Tfinal < 8) {
    next
  }
  
  # Final length of this participant's study period
  # current processed timeline starts at Day 1
  participant_mod <- base_mod
  participant_mod$dims[["T"]] <- Tfinal
  
  iters <- 15000
  conv_tol <- 0.0001
  
  # Predict Day 8 through last day
  for (pred_day in 8:Tfinal) {
    
    # For expanding window starting at Day 1,
    # window length = prediction day
    TT <- pred_day
    
    # Start every window from the initial model
    local_mod <- participant_mod
    local_mod$dims[["TT"]] <- TT
    
    # Day 1 through prediction day
    raw_datamat <- Y[, 1:TT, drop = FALSE]
    
    # Save true lapse
    act_0 <- raw_datamat[10, TT]
    
    # Mask lapse on prediction day
    masked_datamat <- raw_datamat
    masked_datamat[10, TT] <- NA
    local_mod[["data"]] <- masked_datamat
    
    # Fit model
    local_mod <- run_kf(local_mod)
    fit_obj <- run_em(local_mod, iters, conv_tol, zero_priors_bool)
    local_mod <- fit_obj[["model"]]
    
    # Prediction for masked lapse
    pred_0 <- local_mod$kf_proc$ytT[10, TT]
    
    pred_out <- rbind(
      pred_out,
      data.frame(
        subid = subid,
        day = pred_day,
        pred = pred_0,
        actual = act_0
      )
    )
  }
}

# # #
# Above model fitting part basicaly came from 
# function run_mle_fit() in mle_coef_fit.R
# # #

# make prediction for all participants, day 8 - ,
# first ignore time and get overall AUROC(for all participants) and LOG LOSS(all observations)(make negative prediction zero, above 1-1)
# include day 30, second day 30, should be better.(probability on each 30 days and overall)

write.csv(
  pred_out,
  "sherry_repro/pred_out.csv",
  row.names = FALSE
)

