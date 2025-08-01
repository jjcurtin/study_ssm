run_mle_fit <- function(the_subid, data, info, fit_lapse=FALSE){

  subid_info <- info |> 
    filter(subid == the_subid)
  
  subid_data <- data |> 
    filter(subid == the_subid)
  
  first_ema_day <- 1 
  # KW: labels skip first day so days on study (last ema - first ema) has an extra observation
  # remove 1 from count?
  last_ema_day <- subid_info$study_days - 1
  ema_count <- subid_info$n
  
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
  # n: dimension of the observation vector - set to 10 (9 EMA variables + 1 lapse variable)
  # g: dimension of the observation noise (same as n, each observation has its own noise term)
  n <- g <- 10
  dims['n'] <- n
  dims['g'] <- g
  
  # Fix hidden state dimension at 2
  # m: dimension of the latent state vector, set to 2
  # h: dimension of the latent state noise, set to 2
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
  settings[['tinit']]<-1 # intialization time (first row of data matrix)
  settings[['subid']] <-subid_info$subid
  settings[['mema_count']] <- ema_count # number of emas
  settings[['last_ema']] <- last_ema_day
  settings[['first_ema']] <- first_ema_day
  
  mod[['settings']]<-settings
  
  # # # # # # # # # # # #
  # Initial Parameters  #
  # # # # # # # # # # # #
  
  # Initialize the parameter values (uses a fixed value obtained as a population average of early MARSS MLE fits)
  # returns:
  #  initial states (x0, v0)
  #  transition equation, offset, variance/covariance, coefficient matrix on noise
  #  observation equation, offset, variance/covariance, coefficient matrix on noise
  mod[['par']] <- init_par(mod[['dims']])
  
  # # # # # # # # # # # #
  #        Priors       #
  # # # # # # # # # # # #
  priors <- make_zero_priors(mod)
  mod[['priors']] <- priors
  zero_priors_bool <- TRUE
  
  # # # # # # # # # # # #
  #     Data Horizon    #
  # # # # # # # # # # # #
  # Check for the length of the time series, add to dims
  # Tfinal = final time series index( also equals length of data frame)
  # TT = The number of time steps (i.e., duration or length of the time series)

  Tfinal <- min(last_ema_day,89) 
  start_day <- first_ema_day 
  mod$dims['T'] <- Tfinal 
  mod$dims['TT'] <- Tfinal - start_day
  
  # # # # # # # # # # # #
  #     Attach Data     #
  # # # # # # # # # # # #
  # Put in the data
  datamat <- t(as.matrix(subid_data))[c(3, 4,5,6,7,8,9,10,11,12),start_day:Tfinal]
  
  mod[['data']]<-datamat
  
  # # # # # # # # # # # # # # # # # # # # #
  #             MODEL FITTING             #
  # # # # # # # # # # # # # # # # # # # # #
  local_mod <- mod
  
  # # # # # # # # # # # #
  #     Fit/Predict     #
  # # # # # # # # # # # #
  
  # Fit model
  # error 1: TT needs to be number of data rows
  
  # error 2 with SSMCustom: 
  # Error in `[[<-.data.frame`(`*tmp*`, i, value = c(25L, 25L, 25L, 25L, 25L,  : 
  # replacement has 830 rows, data has 83
  local_mod <- run_kf(local_mod) # MARSSkf() as another option?
  iters <- 15000
  conv_tol <- 0.0001
  fit_obj <- run_em(local_mod,iters,conv_tol,zero_priors_bool)
  local_mod <- fit_obj[['model']]
  return_list <- list()
  return_list[['fit']] <- local_mod
  need_flip <- check_mle_flip(local_mod$par)
  return_list[['coefs']] <- make_mle_coef_frame(local_mod,need_flip,fit_lapse)
  return_list[['fitting_error']] <- fit_obj[['error']]
  return(return_list)
}

make_zero_priors <- function(model){
  priors <- list()
  n <- model$dims$n
  m <- model$dims$m
  
  # (mu, sigma) parameters for priors on A
  mu_A1 <- rep(0,n)
  mu_A2 <- rep(0,n)
  mu_A <- cbind(mu_A1,mu_A2)
  
  sig_A1 <- rep(0,n)
  sig_A2 <- rep(0,n)
  sig_A <- matrix(0,n,m)
  
  comb_A <- matrix(0,n,m)
  
  priors[['mu_A']] <- mu_A
  priors[['sig_A']] <- sig_A
  priors[['comb_A']] <- comb_A
  
  # (mu, sigma) parameters for priors on B
  mu_B1 <- 0
  mu_B2 <- 0
  
  sig2_B1 <- 0
  sig2_B2 <- 0
  
  mu_B <- diag(c(mu_B1,mu_B2),m)
  sig_B <- matrix(0,m,m)
  
  priors[['mu_B']] <- mu_B
  priors[['sig_B']] <- sig_B
  
  # (alpha,beta) parameters for prior on R
  r_al1 <- -1
  r_al2 <- -1
  r_al3 <- -1
  r_al4 <- -1
  r_al5 <- -1
  r_al6 <- -1
  r_al7 <- -1
  r_al8 <- -1
  r_al9 <- -1
  r_al10 <- -1
  vec_r_al <- c(r_al1,r_al2,r_al3,r_al4,r_al5,r_al6,r_al7,r_al8,r_al9,r_al10)
  r_alpha_tilde <- diag(vec_r_al+1,n,n)
  
  r_bet1 <- 0
  r_bet2 <- 0
  r_bet3 <- 0
  r_bet4 <- 0
  r_bet5 <- 0
  r_bet6 <- 0
  r_bet7 <- 0
  r_bet8 <- 0
  r_bet9 <- 0
  r_bet10 <- 0
  vec_r_bet <- c(r_bet1,r_bet2,r_bet3,r_bet4,r_bet5,r_bet6,r_bet7,r_bet8,r_bet9,r_bet10)
  r_beta_tilde <- diag(vec_r_bet,n,n)
  
  priors[['r_alpha_tilde']]<-r_alpha_tilde
  priors[['r_beta_tilde']]<-r_beta_tilde
  
  # Priors on c
  # Normal distributions (mu_c,sig2_c)
  mu_c1 <- 0
  mu_c2 <- 0
  mu_c3 <- 0
  mu_c4 <- 0
  mu_c5 <- 0
  mu_c6 <- 0
  mu_c7 <- 0
  mu_c8 <- 0
  mu_c9 <- 0
  mu_c10 <- 0
  mu_c <- matrix(c(mu_c1,mu_c2,mu_c3,mu_c4,mu_c5,mu_c6,mu_c7,mu_c8,mu_c9,mu_c10),n,1)
  
  sig2_c1 <- 0
  sig2_c2 <- 0
  sig2_c3 <- 0
  sig2_c4 <- 0
  sig2_c5 <- 0
  sig2_c6 <- 0
  sig2_c7 <- 0
  sig2_c8 <- 0
  sig2_c9 <- 0
  sig2_c10 <- 0
  #sig_c <- diag(c(sig2_c1,sig2_c2,sig2_c3,sig2_c4,sig2_c5,sig2_c6,sig2_c7,sig2_c8,sig2_c9,sig2_c10),m)
  sig_c <- matrix(0,n,n)
  
  priors[['mu_c']] <- mu_c
  priors[['sig_c']] <- sig_c
  
  # Priors on d
  # Normal distributions (mu_d,sig2_d)
  mu_d1 <- 0
  mu_d2 <- 0
  mu_d <- matrix(c(mu_d1,mu_d2),m,1)
  
  sig2_d1 <- 0
  sig2_d2 <- 0
  #Sig_d <- diag(c(sig2_d1,sig2_d2),m)
  #Sig_inv_d <- solve(Sig_d)
  sig_d <- matrix(0,m,m)
  
  priors[['mu_d']] <- mu_d
  priors[['sig_d']] <- sig_d
  
  return(priors)
}

check_mle_flip <- function(coefmat){
  A1 <- coefmat$A[,1]
  A2 <- coefmat$A[,2]
  nA1 <-norm(A1,type='2')
  nA2 <-norm(A2,type='2')
  #print(paste(nA1,nA2))
  switch <- (nA2>nA1)
  return(switch)
}

make_mle_coef_frame <- function(model,need_flip,fit_lapse=FALSE){
  # Pull information from model
  pars <- model$par
  n <- model$dims$n
  m <- model$dims$m
  
  # Set up placeholder lists for storing information and later export
  A <- list()
  B <- list()
  c <- list()
  d <- list()
  R <- list()
  
  parA <- pars$A
  parB <- pars$B
  pard <- pars$d
  
  # Flip parameters as needed (A, B, and d)
  if(need_flip){
    parA <- matrix(0,n,m)
    parA[,1] <- pars$A[,2]
    parA[,2] <- pars$A[,1]
    
    parB <- matrix(0,m,m)
    parB[1,1] <- pars$B[2,2]
    parB[2,2] <- pars$B[1,1]
    
    pard <- matrix(NA,m,1)
    pard[1] <- pars$d[2]
    pard[2] <- pars$d[1]
  }
  
  # Write out values to lists using (potentially flipped) parameter matrices
  for (j in 1:m){
    for (i in 1:n){
      A[[paste("z",i,j,sep='')]] <- parA[i,j]
      if (j==1){
        c[[paste("a",i,sep='')]] <- pars$c[i]
        if (i < n){
          R[[paste("r",i,sep='')]] <- pars$R[i,i]
        }
      }
    }
    B[[paste("b",j,j,sep='')]] <- parB[j,j]
    d[[paste("u",j,sep='')]] <- pard[j]
    
  }
  # Assemble into dataframe, including subid info at the beginning
  info_frame <- data.frame(id=model$settings$subid,mema_count=model$settings$mema_count,fit_var=fit_lapse)
  out<-cbind(info_frame,data.frame(A),data.frame(c),data.frame(R),data.frame(B),data.frame(d))
  return(out)
}
