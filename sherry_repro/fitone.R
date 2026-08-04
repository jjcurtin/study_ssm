# This script is to fit SSM on one patient's data(subid=2)
# and predict day83 through day81-82
# Relevant libraries
library("KFAS")
library('tidyverse')
library("MASS",exclude=c('select'))
library("data.table")
library("tictoc")
library('MARSS')

# Source helper files
source('sherry_repro/helper_functions.R')
source('sherry_repro/mle_fit/mle_coef_fit.R')

#--------------------------------------------------
# Read one participant
#--------------------------------------------------
data_all <- read.csv("S:/ssm_capstone/data_raw/features_ema_ssm_1x_day_24h.csv")
subject_data <- subset(data_all, subid == 2)
unique(subject_data$lapse)
subject_data$lapse <- ifelse(subject_data$lapse == "lapse", 1, 0)

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

# Observation Matrix (g x TT)
Y <- t(as.matrix(subject_data[, obs_cols]))
dim(Y)

# # # # # # # # # # # #
#     Data Horizon    #
# # # # # # # # # # # #

# Length of this participant's time series
TT <- ncol(Y)
mod$dims[['TT']] <- TT

# # # # # # # # # # # #
#     Attach Data     #
# # # # # # # # # # # #

# Keep a copy of the complete observation matrix
raw_datamat <- Y

# Create masked data matrix for prediction
masked_datamat <- raw_datamat

# Mask the final lapse state
masked_datamat[10, TT] <- NA

# Attach to model object
mod[['data']] <- masked_datamat

# Save the true final lapse for later comparison
act_0 <- raw_datamat[10, TT]

# # # # # # # # # # # # # # # # # # # # #
#             MODEL FITTING             #
# # # # # # # # # # # # # # # # # # # # #

# ZERO PRIORS (placeholder so run_em() works)
priors <- make_zero_priors(mod)
mod[['priors']] <- priors
zero_priors_bool <- TRUE

local_mod <- mod
# Fit model
local_mod <- run_kf(local_mod)
iters <- 15000
conv_tol <- 0.0001
fit_obj <- run_em(local_mod,iters,conv_tol,zero_priors_bool)
local_mod <- fit_obj[['model']]

# # #
# Above model fitting part basicaly came from 
# function run_mle_fit() in mle_coef_fit.R
# # #

fit_obj$error

mod$par$B

fit_obj$model$par$B

pred_0 <- local_mod$kf_proc$ytT[10, TT]
pred_0
act_0
local_mod$kf_proc$xtT[, TT]
A <- local_mod$par$A
c <- local_mod$par$c
x <- local_mod$kf_proc$xtT[, TT]

(A %*% matrix(x, 2, 1) + c)
