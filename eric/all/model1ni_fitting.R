rm(list=ls())
library('tidyverse')
library('MARSS')
library('foreach')
library('doParallel')
library('data.table')
library('tictoc')
library('mvtnorm')

# Import the data
path <- '/Users/eric/repos/aud/data/subject_marss_60_90_covar.csv'
data <- read.csv(path)
subid_list <- unique(data$subid)
use_covar <- 'no'
tic("runtime")


fit_model <- function(input_data,input_model){
  em_fit <- MARSS(input_data, model = input_model, control = list(maxit = 300),silent=TRUE)
  if(any(em_fit$errors %like% 'Stopped at iter')){
    bfgs_fit <- MARSS(input_data, model = input_model, method = "BFGS")
  } else {
    bfgs_fit <- MARSS(input_data, model = input_model, method = "BFGS", inits = em_fit)
  }
  #bfgsfit <- MARSS(input_data, model = input_model, method = "BFGS", inits = em_fit,silent=TRUE)
  return(bfgs_fit)
}

make_pred_frame <- function(model_name, subject, AIC_val,horizon,t,pred,actual,use_covar){
  df <- data.frame("model_name"=model_name,"id"=subject,"AIC"=AIC_val,"train_horizon"=horizon,"pred_t"=t,"pred"=pred,"actual"=actual,"covar"=use_covar)
  return(df)
}
make_coef_frame <- function(model_name, subject, horizon,coefs,flip,use_covar){
  gendf <- data.frame("model_name"=model_name,"id"=subject,"train_horizon"=horizon,"covar"=use_covar)
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

check_flip <- function(coefmat){
  z1 <- coefmat$Z[,1]
  z2 <- coefmat$Z[,2]
  nz1 <-norm(z1,type='2')
  #print(nz1)
  nz2 <-norm(z2,type='2')
  #print(nz2)
  switch <- (nz2>nz1)
  #print(switch)
  return(switch)
}

# Function for forecasting the next step
fr_step <-function(B,Urep,x,w){
  xt1 <- B %*% x + Urep + w
  return(xt1)
}
# Function for returning the clipped lapse probability
d_step <- function(Z10,A10rep,x){
  #d <- Z10 %*% x + A10rep + v
  d <- Z10 %*% x + A10rep
  return(d)
}
# Function for clipping predictions
clip <- function(d){
  eps <- 0.000000000001
  if(d < eps){
    return(eps)
  }
  if(d > 1-eps){
    return(1-eps)
  }
  return(d)
}
# Function for returning a vector of lapse probabilities for the given window
calc_lapse_window <- function(mat){
  prod_val <- apply(mat,2,prod)
  p_window_lapse <- 1-prod_val
  return(p_window_lapse)
}
####
get_window_preds <- function(modelfit,window_length){
  n<- 10000
  # Create list for storing results
  fr_x <- list()
  fr_d <- list()
  
  # Get the coefficients from latest fit
  coefs <- coef(modelfit,type='matrix')
  B <- coefs$B
  U <- coefs$U
  Z <- coefs$Z
  A <- coefs$A
  A10 <- A[10]
  Z10 <- Z[10,]
  
  # Tile the relevant vectors for matrix processing
  Urep<- t(matrix(U, nrow=n, ncol=length(U), byrow=TRUE))
  A10rep <- t(matrix(A10, nrow=n, ncol=length(A10), byrow=TRUE))
  
  # Pull out the mean and variance of the state estimate for the end of the trajectory
  # Get the index of the final entry
  data_length = dim(modelfit[['marss']]$data)[2]
  # Run the kalman filtering/smoothing (on the model with the extended dataset, both training + new data)
  kf <- MARSSkf(modelfit)
  # Pull out the mean and variance
  x90m<-c(kf$xtt[,data_length])
  print(x90m)
  x90v<-matrix(kf$Vtt[,,data_length])
  print(as.matrix(x90v))
  
  # Generate the set of initial state estimates and store in list
  xt0 <- t(rmvnorm(n,mean=x90m,sigma=x90v))
  
  # Store the initial states as the first entry of the fr_x list
  fr_x[[1]] <- xt0
  # Store same day prediction values (complement of probability of no lapse)
  fr_d[[1]] <- 1-vapply(d_step(Z10,A10rep,xt0),clip,numeric(1))
  
  # Loop over desired prediction horizon
  # Next day is denoted as 2 (since same-day is 1 in the lists)
  start_val <- 2
  # Window length > 1 denotes future predictions (beyond same-day)
  if(window_length>1){
    # Loop from 2 to the final day of the window
    end_val <- start_val+window_length-2
    for(i in start_val:end_val){
      # Use the previous day's x value (each trajectory is independent)
      x_iter <- fr_x[[i-1]]
      # Generate identity noise for every trajectory
      w_iter <- t(rmvnorm(n,rep(0,1),sigma=diag(1)))
      # Create current step's hidden state using the fr_step function (just a step forward in the transition equation)
      xt <- fr_step(B,Urep,x_iter,w_iter)
      # Store the current step's hidden state
      fr_x[[i]] <- xt
      # Noiseless evaluation of the unclipped lapse probability
      dt <- d_step(Z10,A10rep,xt)
      # Clip the lapse probability
      clipped_dt <- vapply(dt,clip,numeric(1))
      # Store the probability of no lapse
      comp_dt <- 1-clipped_dt
      fr_d[[i]] <- comp_dt
    }
  }
  
  final <- list()
  lapse_mat <- do.call("rbind",fr_d)
  for(i in 1:length(fr_d)){
    mat <- matrix(lapse_mat[1:i,],nrow=i)
    final[[i]] <- calc_lapse_window(mat)
  }
  return(final)
}

####

n.cores <- detectCores() - 1
my.cluster <- makeCluster(
  n.cores,
  type="FORK"
)

doParallel::registerDoParallel(cl=my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

# Use this model
model_name<- "1ni"
window_length <- 8
prob_col_names <- c("w0_pred","w1_pred","w2_pred","w3_pred","w4_pred","w5_pred","w6_pred","w7_pred")
act_col_names <- c("w0_act","w1_act","w2_act","w3_act","w4_act","w5_act","w6_act","w7_act")


# Solve using these training horizons
#horizon_list <- c(80)
horizon_list <- c(40,60,80)
# Loop over subjects (parallelized)
parout<-foreach(subject_idx=1:length(subid_list)) %dopar% {
#parout<-foreach(subject_idx=1:2) %dopar% {
#for(subject_idx in 1:2){
  # Get the current subject id
  current_id <- subid_list[[subject_idx]]
  # Pull the data associated with the current individual
  current_data <- data[data$subid==current_id,]
  # Convert to a matrix and trim extra rows
  #datamat <- t(as.matrix(current_data))
  #datamat <- datamat[-c(1,2,3,14),]
  datamat <- t(as.matrix(current_data))
  datamat <- datamat[-c(1,2,3,14,15,16,17,18,19,20),]
  covmat <- t(as.matrix(current_data))
  covmat <- covmat[14:19,]
  # Create a storage dataframe for this subject
  subject_pred_store <-data.frame()
  subject_coef_store <-data.frame()
  # Loop over each horizon and train/test the model for each one
  for(horizon_idx in 1:length(horizon_list)){
    # Get the specific horizon value
    horizon<-horizon_list[horizon_idx]
    
    # Set up the complete data from 1:horizon as model fitting training set
    training <- datamat[,1:horizon]
    # Specify the model using the covariates
    # Model definition
    init_var<-1
    # 1 latent states
    B1 <- matrix(list("b"),1,1)
    U1 <- matrix(list("u1"),1,1)
    Q1 <- diag(1, 1)
    Z1 <- matrix(c(
      "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101"
    ), 10, 1)
    A1 <- matrix(list("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10"),10,1)
    # R1 <- matrix(list(
    #   "r1", 0, 0, 0, 0, 0, 0, 0, 0, 0,
    #   0, "r2", 0, 0, 0, 0, 0, 0, 0, 0,
    #   0, 0, "r3", 0, 0, 0, 0, 0, 0, 0,
    #   0, 0, 0, "r4", 0, 0, 0, 0, 0, 0,
    #   0, 0, 0, 0, "r5", 0, 0, 0, 0, 0,
    #   0, 0, 0, 0, 0, "r6", 0, 0, 0, 0,
    #   0, 0, 0, 0, 0, 0, "r7", 0, 0, 0,
    #   0, 0, 0, 0, 0, 0, 0, "r8", 0, 0,
    #   0, 0, 0, 0, 0, 0, 0, 0, "r9", 0,
    #   0, 0, 0, 0, 0, 0, 0, 0, 0, "r10"
    # ),10,10)
    R1 <- "unconstrained"
    pi1 <- matrix(0, 1, 1)
    V1 = diag(init_var, 1, 1)
    if(use_covar=='x'){
      chosen_model<-model1ni <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1,C=C1,c=covmat[,1:horizon], x0 = pi1, V0 = V1, tinitx = 1)
    } else if (use_covar=='y') {
      chosen_model<-model1ni <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1,D=D1,d=covmat[,1:horizon], x0 = pi1, V0 = V1, tinitx = 1)
    }else if (use_covar=='lapse') {
      chosen_model<-model1ni <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1,D=Dlapse,d=covmat[,1:horizon], x0 = pi1, V0 = V1, tinitx = 1)
    }else {
      chosen_model<-model1ni <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)
    }
    
    
    # Fit the model
    #print("prefit")
    mlefit <- fit_model(training,chosen_model)
    #print("postfit")
    coefs<-coef(mlefit,type='list')
    coefmat<-coef(mlefit,type='matrix')
    aic_val<-AIC(mlefit)
    
    # Maximum horizon length for prediction steps
    mhl <- 90-horizon
    horizon_pred_store <- data.frame()
    # Add 1 day of new data at a time and make same-day lapse prediction
    for(hlength in 1:mhl){
      if(mhl>1){
        # Same day prediction
        # Start index for new data
        start <- horizon+1
        # End index for new data
        end <- horizon+hlength
        # Use start and end to create a new data matrix
        new_dat <-as.matrix(datamat[,start:end])
        new_cov <-as.matrix(covmat[,start:end])
        # Pull the actual lapse value out of the new data matrix
        actual_lapse <- new_dat[[10,ncol(new_dat)]]
        # Put in NA for use in forecasting
        new_dat[10,ncol(new_dat)]<-NA
        # Make the prediction (uses the new data to smooth xt-1 and get filtered xt?)
        if(use_covar=='x'){
          fr<-predict(mlefit,n.ahead=hlength,newdata=list(y=new_dat,c=new_cov),type='ytT')
        } else if(use_covar=='y'|use_covar=='lapse'){
          fr<-predict(mlefit,n.ahead=hlength,newdata=list(y=new_dat,d=new_cov),type='ytT')
        } else {
          fr<-predict(mlefit,n.ahead=hlength,newdata=list(y=new_dat),type='ytT')
        }
        output <- fr$pred
        pred_val <- fr$pred[nrow(output),ncol(output)]
        
        # Window prediction
        window_data <- cbind(training,new_dat)
        mlefit_w <- mlefit
        mlefit_w[['marss']]$data <- window_data
        mlefit_w[['model']]$data <- window_data
        attributes(mlefit_w[['marss']])$model.dims$data <- dim(window_data)
        
        window_probs <- vapply(get_window_preds(mlefit_w,window_length),mean,numeric(1))
        probs_df <- data.frame(t(window_probs))

        actuals_mat <- matrix(NA, nrow = 1, ncol = window_length)
        actuals_df <- data.frame(actuals_mat)
        for(j in 1:window_length){
          if((end+j-1)<=90){
            relevant_lapse_data <- datamat[10,end:(end+j-1)]
            subset <- relevant_lapse_data[1:j]
            lapse_bool <- as.numeric(1 %in% subset)
            actuals_df[j]<-lapse_bool
          } else {
            probs_df[j] <- NA
          }
          
        }
        # if((end+window_length-1) <= 90){
        #   window_probs <- vapply(get_window_preds(mlefit_w,window_length),mean,numeric(1))
        #   probs_df <- data.frame(t(window_probs))
        #   relevant_lapse_data <- datamat[10,start:(start+window_length-1)]
        #   actuals_mat <- matrix(NA, nrow = 1, ncol = window_length)
        #   for(j in 1:window_length){
        #     subset <- relevant_lapse_data[1:j]
        #     lapse_bool <- as.numeric(1 %in% subset)
        #     actuals_mat[1,j]<-lapse_bool
        #   }
        #   actuals_df<-data.frame(actuals_mat)
        # }else{
        #   probs_df<- data.frame(matrix(NA, nrow = 1, ncol = window_length))
        #   actuals_df <-data.frame(matrix(NA, nrow = 1, ncol = window_length))
        # }
        colnames(probs_df)<-prob_col_names
        colnames(actuals_df)<-act_col_names
      }
      else{
        end <- NA
        actual_lapse <-NA
        pred_val <- NA
        
        probs_df<- data.frame(matrix(NA, nrow = 1, ncol = window_length))
        actuals_df <-data.frame(matrix(NA, nrow = 1, ncol = window_length))
        colnames(probs_df)<-prob_col_names
        colnames(actuals_df)<-act_col_names
      }
      iter_pred_df <- make_pred_frame(model_name,current_id, aic_val,horizon,end,pred_val,actual_lapse,use_covar)
      iter_pred_df <- cbind(iter_pred_df,probs_df,actuals_df)
      horizon_pred_store<-rbind(horizon_pred_store,iter_pred_df)
    }
    subject_pred_store<-rbind(subject_pred_store,horizon_pred_store)
    #check_flip(coefmat)
    #subject_coef_store<-rbind(subject_coef_store,make_coef_frame(model_name,current_id, horizon,coefs,check_flip(coefmat),use_covar))
  }
  return(list(subject_pred_store,subject_coef_store))
}
pred_out <- as.data.frame(rbindlist(lapply(parout,function(x){x[[1]]})))
#coef_out <- rbindlist(lapply(parout,function(x){x[[2]]}))
write.csv(pred_out, "/Users/eric/repos/aud/model_eval/1ni/pred_nocovar_window_fullvarest.csv", row.names=FALSE)
#write.csv(coef_out, "/Users/eric/repos/aud/model_eval/1ni/coef_nocovar_window.csv", row.names=FALSE)
stopCluster(cl = my.cluster)
toc()