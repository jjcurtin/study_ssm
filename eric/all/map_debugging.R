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

# Source helper files
source('aud_helper_functions.R')

infopath <- '/Users/eric/repos/aud/data/subject_info.csv'
raw_info <- read.csv(infopath)

subjects = raw_info$subid

for(idx in 1:length(subjects)){
  subid <- subjects[idx]
  print(subid)
  
  mod = list()
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
  prior_path
  raw_prior <- read.csv('/Users/eric/repos/aud/data/mle_long_coef_fits_lapse_not_fitted.csv')
  
  priors <- make_priors(mod, raw_prior, subid)
}