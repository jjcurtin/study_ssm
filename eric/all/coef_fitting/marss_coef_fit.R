library('tidyverse')
library('MARSS')
library('data.table')
run_2ni_fit <- function(input_data,info,fit_lapse=FALSE){
  # set the horizon based on info (note conversion from 0 to 1 indexing)
  last_ema_day <- min(info$last_morning_ema_day+1,90)
  #print(last_ema_day)
  # Dataframe modifications
  data_select <- select(input_data,
                        ema_2,ema_3,ema_4,ema_5,ema_6,ema_7,ema_8,ema_9,ema_10,lapse)
  
  datamat <- t(as.matrix(data_select))

  # pull relevant columns
  datamat <- datamat[,1:last_ema_day]
  
  # Model definition, specific to 2NI
  B1 <- matrix(list(0),2,2)
  diag(B1)<-list("b11","b22")
  U1 <- matrix(list("u1","u2"),2,1)
  Q1 <- diag(1, 2)
  Z1 <- matrix(c(
    "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101",
    "z12", "z22", "z32", "z42", "z52", "z62", "z72", "z82", "z92", "z102"
  ), 10, 2)
  A1 <- matrix(list("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10"),10,1)
  R1 <- matrix(list(
    "r1", 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, "r2", 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, "r3", 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, "r4", 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, "r5", 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, "r6", 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, "r7", 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, "r8", 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, "r9", 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, "r10"
  ),10,10)
  if(!fit_lapse){
    R1[10,10] <- 0.25
  }
  pi1 <- matrix(0, 2, 1)
  V1 = diag(1, 2, 2)
  model<-model2ni <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)
  
  em_fitting_error <- FALSE
  # em_fit <- MARSS(datamat, model = model, control = list(maxit = 1000))
  # if(any(em_fit$errors %like% 'Stopped at iter')){
  #   bfgs_fit <- MARSS(datamat, model = model, method = "BFGS")
  #   em_fitting_error <- TRUE
  # } else {
  #   bfgs_fit <- MARSS(datamat, model = model, method = "BFGS", inits = em_fit)
  # }
  em_fit <- try(MARSS(datamat, model = model, control = list(maxit = 10000)),silent=TRUE)
  if(inherits(em_fit,"try-error")){
    print("encountered fatal error on EM fit, proceeding to BFGS with default init")
    bfgs_fit <- MARSS(datamat, model = model, method = "BFGS",silent=TRUE)
    em_fitting_error <- TRUE
  } else{
    print("no fatal error in EM fit, using EM as BFGS init if convergence indicates no fitting issues")
    if(em_fit[['convergence']]==2 || em_fit[['convergence']]==52){
      print("convergence indicated fitting issues")
      bfgs_fit <- MARSS(datamat, model = model, method = "BFGS",silent=TRUE)
      em_fitting_error <- TRUE
    } else {
      print("convergence indicated no fitting issues")
      bfgs_fit <- MARSS(datamat, model = model, method = "BFGS", inits = em_fit,silent=TRUE)
    }
  }
  coefmat<-coef(bfgs_fit,type='matrix')
  coeflist <-coef(bfgs_fit,type='list')
  need_flip <- check_flip(coefmat)
  subid_coef <- make_coef_frame(info,coeflist,need_flip,fit_lapse)
  return_list <- list()
  return_list[['fit']] <- bfgs_fit
  return_list[['coefs']] <- subid_coef
  return_list[['fitting_error']] <- em_fitting_error
  return(return_list)
}

check_flip <- function(coefmat){
  z1 <- coefmat$Z[,1]
  z2 <- coefmat$Z[,2]
  nz1 <-norm(z1,type='2')
  nz2 <-norm(z2,type='2')
  switch <- (nz2>nz1)
  return(switch)
}

make_coef_frame <- function(info,coefs,flip,fit_var){
  id <- info$subid
  mema_count <- info$mema_count
  gendf <- data.frame("id"=id,"mema_count"=mema_count,"fit_var"=fit_var)
  # Use convention that larger rank Z column determines feature order
  if(flip){
    # Pull in the original coefficient df
    rawcoefdf <-data.frame(t(do.call(rbind,coefs)))
    # Make a copy for modification and final export
    coefdf <- rawcoefdf
    # Switch the diagonal B values
    coefdf$b11<-rawcoefdf$b22
    coefdf$b22<-rawcoefdf$b11
    # Switch the Z values (swap columns)
    # One side of the swap
    coefdf$z11<-rawcoefdf$z12
    coefdf$z21<-rawcoefdf$z22
    coefdf$z31<-rawcoefdf$z32
    coefdf$z41<-rawcoefdf$z42
    coefdf$z51<-rawcoefdf$z52
    coefdf$z61<-rawcoefdf$z62
    coefdf$z71<-rawcoefdf$z72
    coefdf$z81<-rawcoefdf$z82
    coefdf$z91<-rawcoefdf$z92
    coefdf$z101<-rawcoefdf$z102
    # Other side of the swap
    coefdf$z12<-rawcoefdf$z11
    coefdf$z22<-rawcoefdf$z21
    coefdf$z32<-rawcoefdf$z31
    coefdf$z42<-rawcoefdf$z41
    coefdf$z52<-rawcoefdf$z51
    coefdf$z62<-rawcoefdf$z61
    coefdf$z72<-rawcoefdf$z71
    coefdf$z82<-rawcoefdf$z81
    coefdf$z92<-rawcoefdf$z91
    coefdf$z102<-rawcoefdf$z101
    # Switch the U values
    coefdf$u1 <-rawcoefdf$u2
    coefdf$u2 <-rawcoefdf$u1
  }
  else{
    coefdf <- data.frame(t(do.call(rbind,coefs)))
  }
  df <- cbind(gendf,coefdf)
  return(df)
}

# datapath <- '/Users/eric/repos/aud/data/day_labels.csv'
# infopath <- '/Users/eric/repos/aud/data/subject_info.csv'
# data <- read.csv(datapath)
# info <- read.csv(infopath)
# test_subid <- 98
# fit_lapse_var <- FALSE
# subid_data <- data[data$subid==test_subid,]
# subid_info <- info[info$subid==test_subid,]
# #debug(MARSSkem)
# fit_info <- run_2ni_fit(subid_data,subid_info,fit_lapse_var)
# # #print(test_subid)
# # #print(fit_info$fitting_error)