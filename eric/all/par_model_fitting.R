library('tidyverse')
library('MARSS')
library('foreach')
library('doParallel')
library('tictoc')

# Import the data
path <- '/Users/eric/repos/aud/subject_marss_60_90.csv'
data <- read.csv(path)
subid_list <- unique(data$subid)
tic("runtime")
# Model definition
init_var<-10

# 2 latent states with state interactions
# B1 <- matrix(list("b11", "b21", "b12", "b22"), 2, 2)
B1 <- "unequal"
U1 <- "unequal"
Q1 <- diag(1, 2)
Z1 <- matrix(c(
  "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101",
  "z12", "z22", "z32", "z42", "z52", "z62", "z72", "z82", "z92", "z102"
), 10, 2)
A1 <- "unequal"
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
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25
),10,10)
pi1 <- matrix(0, 2, 1)
V1 = diag(init_var, 2, 2)
model2i <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

# 3 latent states with state interactions
B1 <- matrix(list("b11", "b21","b31", "b12", "b22","b32","b13","b23","b33"), 3, 3)
U1 <- "unequal"
Q1 <- diag(1, 3)
Z1 <- matrix(c(
  "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101",
  "z12", "z22", "z32", "z42", "z52", "z62", "z72", "z82", "z92", "z102",
  "z13", "z23", "z33", "z43", "z53", "z63", "z73", "z83", "z93", "z103"
), 10, 3)
A1 <- "unequal"
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
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25
),10,10)
pi1 <- matrix(0,3,1)
V1 = diag(init_var, 3, 3)
model3i <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

# 4 latent states with state interactions
B1 <- "unequal"
U1 <- "unequal"
Q1 <- diag(1, 4)
Z1 <- matrix(c(
  "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101",
  "z12", "z22", "z32", "z42", "z52", "z62", "z72", "z82", "z92", "z102",
  "z13", "z23", "z33", "z43", "z53", "z63", "z73", "z83", "z93", "z103",
  "z14", "z24", "z34", "z44", "z54", "z64", "z74", "z84", "z94", "z104"
), 10, 4)
A1 <- "unequal"
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
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25
),10,10)
pi1 <- matrix(0,4,1)
V1 = diag(init_var, 4, 4)
model4i <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

# 1 latent states
B1 <- matrix(list("b"),1,1)
U1 <- matrix(list("u1"),1,1)
Q1 <- diag(1, 1)
Z1 <- matrix(c(
  "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101"
), 10, 1)
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
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25
),10,10)
pi1 <- matrix(0, 1, 1)
V1 = diag(init_var, 1, 1)
model1ni <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

# 2 latent states, no state interactions
B1 <- matrix(list(0),2,2)
diag(B1)<-list("b1","b2")
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
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25
),10,10)
pi1 <- matrix(0, 2, 1)
V1 = diag(init_var, 2, 2)
model2ni <- list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

# 3 latent states, no state interactions
B1 <- matrix(list(0),3,3)
diag(B1)<-list("b1","b2","b3")
U1 <- matrix(list("u1","u2","u3"),3,1)
Q1 <- diag(1, 3)
R1 <- matrix(list(0),10,10)
diag(R1)<-list("r1","r2","r3","r4","r5","r6","r7","r8","r9",0.25)
Z1 <- matrix(c(
  "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101",
  "z12", "z22", "z32", "z42", "z52", "z62", "z72", "z82", "z92", "z102",
  "z13", "z23", "z33", "z43", "z53", "z63", "z73", "z83", "z93", "z103"
), 10, 3)
A1 <- "unconstrained"
pi1 <- matrix(0, 3, 1)
V1 = diag(init_var, 3, 3)
model3ni <-list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

# 4 latent states, no state interactions
B1 <- matrix(list(0),4,4)
diag(B1)<-list("b1","b2","b3","b4")
U1 <- "unequal"
Q1 <- diag(1, 4)
R1 <- matrix(list(0),10,10)
diag(R1)<-list("r1","r2","r3","r4","r5","r6","r7","r8","r9",0.25)
Z1 <- matrix(c(
  "z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z91", "z101",
  "z12", "z22", "z32", "z42", "z52", "z62", "z72", "z82", "z92", "z102",
  "z13", "z23", "z33", "z43", "z53", "z63", "z73", "z83", "z93", "z103",
  "z14", "z24", "z34", "z44", "z54", "z64", "z74", "z84", "z94", "z104"
), 10, 4)
A1 <- "unconstrained"
pi1 <- matrix(0, 4, 1)
V1 = diag(init_var, 4, 4)
model4ni <-list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

# 9 latent states, no state interactions
B1 <- "diagonal and unequal"
U1 <- "unequal"
Q1 <- diag(1, 9)
R1 <- matrix(list(0),10,10)
diag(R1)<-list("r1","r2","r3","r4","r5","r6","r7","r8","r9",0.25)
Z1 <- diag(1,10,9)
Z1[10,] <- c("d1","d2","d3","d4","d5","d6","d7","d8","d9")
A1 <- "unequal"
pi1 <- matrix(0, 9, 1)
V1 = diag(init_var, 9, 9)
model9ni <-list(B = B1, U = U1, Q = Q1, Z = Z1, A = A1, R = R1, x0 = pi1, V0 = V1, tinitx = 1)

fit_model <- function(input_data,input_model){
  em_fit <- MARSS(input_data, model = input_model, control = list(maxit = 200),silent=TRUE)
  bfgsfit <- MARSS(input_data, model = input_model, method = "BFGS", inits = em_fit,silent=TRUE)
  return(bfgsfit)
}
# fit_model <- function(input_data,input_model){
#   try1<-try(em_fit <- MARSS(input_data, model = input_model, control = list(maxit = 200)),silent=TRUE)
#   #try2<-try(bfgsfit <- MARSS(input_data, model = input_model, method = "BFGS", inits = em_fit),silent=TRUE)
#   #print(try1)
#   #print(try2)
#   return(bfgsfit)
# }
backup_fit_model <- function(input_data,input_model){
  bfgsfit <- MARSS(input_data, model = input_model, method = "BFGS")
  return(bfgsfit)
}

make_frame <- function(model_name, subject, AIC_val,horizon,t,pred,actual,coefs){
  gendf <- data.frame("model_name"=model_name,"id"=subject,"AIC"=AIC_val,"train_horizon"=horizon,"pred_t"=t,"pred"=pred,"actual"=actual)
  coefdf <- data.frame(t(do.call(rbind,coefs)))
  df <- cbind(gendf,coefdf)
  return(df)
}

n.cores <- detectCores() - 1
my.cluster <- makeCluster(
  n.cores,
  type="FORK"
)

doParallel::registerDoParallel(cl=my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

# Use this model
model_name<- "2ni"
chosen_model<- model2ni

# Solve using these training horizons
horizon_list <- c(40,60,80)
#horizon_list <- c(80)
# Loop over subjects (parallelized)
parout<-foreach(subject_idx=1:length(subid_list),.combine=rbind.data.frame) %dopar% {
#parout<-foreach(subject_idx=17,.combine=rbind.data.frame) %dopar% {
#for(subject_idx in c(17,18)){
  # Get the current subject id
  current_id <- subid_list[[subject_idx]]
  # Pull the data associated with the current individual
  current_data <- data[data$subid==current_id,]
  # Convert to a matrix and trim extra rows
  datamat <- t(as.matrix(current_data))
  datamat <- datamat[-c(1,2,3,14),]
  # Create a storage dataframe for this subject
  subject_store <-data.frame()
  
  # Loop over each horizon and train/test the model for each one
  for(horizon_idx in 1:length(horizon_list)){
    # Get the specific horizon value
    horizon<-horizon_list[horizon_idx]
    
    # Set up the complete data from 1:horizon as model fitting training set
    training <- datamat[,1:horizon]
    
    # Fit the model
    mlefit <- fit_model(training,chosen_model)
    coefs <- coefs<-coef(mlefit,type='list')
    # error_bool=FALSE
    # tryCatch({mlefit <- fit_model(training,chosen_model)},
    #          message=function(m){print("we see a message");next},
    #          warning=function(w){print("we see a warning");next},
    #          error=function(e){print("we see an error");error_bool<<-TRUE}
    # )
    # if(error_bool){
    #   next
    # }
    aic_val<-AIC(mlefit)
    
    # Maximum horizon length for prediction steps
    mhl <- 90-horizon
    #mhl <- 1
    horizon_store <- data.frame()
    # Add 1 day of new data at a time and make same-day lapse prediction
    for(hlength in 1:mhl){
      # Start index for new data
      start <- horizon+1
      # End index for new data
      end <- horizon+hlength
      # Use start and end to create a new data matrix
      new_dat <-as.matrix(datamat[,start:end])
      # Pull the actual lapse value out of the new data matrix
      actual_lapse <- new_dat[[10,ncol(new_dat)]]
      # Put in NA for use in forecasting
      new_dat[10,ncol(new_dat)]<-NA
      # Make the prediction (uses the new data to smooth xt-1 and get filtered xt?)
      fr<-predict(mlefit,n.ahead=hlength,newdata=list(y=new_dat),type='ytT')
      output <- fr$pred
      pred_val <- fr$pred[nrow(output),ncol(output)]
      iter_df <- make_frame(model_name,current_id, aic_val,horizon,end,pred_val,actual_lapse,coefs)
      horizon_store<-rbind(horizon_store,iter_df)
    }
    subject_store<-rbind(subject_store,horizon_store)
  }
  return(subject_store)
}

write.csv(parout, "/Users/eric/repos/aud/model_eval/2ni6090.csv", row.names=FALSE)
stopCluster(cl = my.cluster)
toc()