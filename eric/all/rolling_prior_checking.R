#!/usr/bin/env Rscript
options(error=function()traceback(2))

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
source('aud_helper_functions.R')
out<-data.frame()

# Pull in the text file
test_list_path <- "/Users/eric/repos/aud/chtc/subid_list_data_eff.txt"
test_list <- read.table(test_list_path)

datapath <- 'data/day_labels.csv'
infopath <- 'data/subject_info.csv'
refpath <- 'data/ssm_fold_reference.json'

# Process relevant json information
ref_fold_info <- read_json(refpath)

for(i in 1:nrow(test_list)){
  iter_row <- test_list[i,]
  keyval <- iter_row$V1
  type <- iter_row$V2
  
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
  #print(last_ema_day)
  #print(mema_count)
  
  
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
  # TEMPORARY, fix observation dimension to 9 EMA + lapse (10)
  # n:= dimension of observation vector
  # TODO g:= first dimension of observation noise coefficient matrix
  n <- g <- 10
  dims['n'] <- n
  dims['g'] <- g
  
  # TEMPORARY - Fix hidden state dimension at 2
  # m:= dimension of latent state vector
  # TODO h:= first dimension of latent state noise coefficient matrix
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
  if (type=='MAP_mle'){
    prior_path<-'mle_coef_fits_lapse_not_fitted.csv'
  } else if (type=='MAP_marss'){
    prior_path<-'marss_long_coef_fits_lapse_not_fitted.csv'
  } else {
    prior_path<-'mle_coef_fits_lapse_not_fitted.csv'
  }
  
  raw_prior <- read.csv(paste('data/',prior_path,sep=''))
  
  priors <- make_priors(mod, raw_prior, subid, TRUE,prior_subids)
  mod[['priors']] <- priors
}