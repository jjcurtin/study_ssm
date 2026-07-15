# Clear workspace
#rm(list=ls())
# Libraries
library("KFAS")
library("comprehenr")
library('tictoc')
# Set up dummy model object
get_dummy_model <- function(data_in){
  # Set up the storage object for the model
  mod <- list()
  
  # Set up data-related parameters
  # List for storing relevant dimensions
  dims <- list()
  # Check for the length of the time series, add to dims
  # TT := final time series index
  TT <- dim(data_in)[2]
  dims['TT']<- TT
  
  # First trim out covariate rows
  dat <- data_in[-c(1,2,3,14,15,16,17,18,19,20),1:TT]
  mod[['data']]<-dat
  
  # Get number of observations left
  # n:= dimension of observation vector
  # TODO g:= first dimension of observation noise coefficient matrix
  n <- g <- dim(dat)[1] 
  dims['n'] <- n
  dims['g'] <- g
  
  # TEMPORARY - Fix hidden state dimension at 2
  # m:= dimension of latent state vector
  # TODO h:= first dimension of latent state noise coefficient matrix
  m <- h  <- 2
  dims['m'] <- m
  dims['h'] <- h 
  
  # Relevant parameters (start with observed means from an early run of parameter estimation)
  # Parent list
  par.1 <- list()
  
  # INITIAL STATES
  par.1[['x0']] <- matrix(c(0,0),2,1)
  par.1[['V0']] <- matrix(c(1,0,0,1),2,2)
  
  # TRANSITION EQUATION
  
  # Transition matrix B (assume diagonal)
  b11 <- 0.2
  b22 <- 0.85
  par.1[["B"]] <- matrix(c(b11,0,0,b22),2,2)
  
  # Transition offset 
  d1 <- -0.24
  d2 <- -0.54
  par.1[["d"]] <- matrix(c(d1,d2),2,1)
  
  # Transition variance covariance
  par.1[["Q"]] <- diag(1,m,m)
  
  # Coefficient matrix on noise
  par.1[["H"]] <- diag(1,h,m)
  
  # OBSERVATION EQUATION
  
  # Observation matrix A
  # Values for each EMA question
  a1 <- c(0.13,0.07)
  a2 <- c(0.14,0.07)
  a3 <- c(0.11,0.06)
  a4 <- c(0.11,0.08)
  a5 <- c(0.07,0.07)
  a6 <- c(0.10,0.07)
  a7 <- c(0.08,0.04)
  a8 <- c(0.09,0.06)
  a9 <- c(0.06,0.04)
  a10 <- c(0.03,0.02)
  # Concatenate all rows together to form observation matrix A
  par.1[["A"]] <-rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
  
  # Observation offset
  # Offset for each EMA question
  c1 <- 0.14
  c2 <- 0.14
  c3 <- 0.16
  c4 <- 0.29
  c5 <- 0.60
  c6 <- 0.40
  c7 <- 0.51
  c8 <- 0.61
  c9 <- 0.21
  c10 <- 0.06
  # Concatenate into vector
  par.1[["c"]] = rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
  
  # Observation noise variance covariance
  r1 <- 0.02
  r2 <- 0.01
  r3 <- 0.03
  r4 <- 0.04
  r5 <- 0.03
  r6 <- 0.05
  r7 <- 0.01
  r8 <- 0.01
  r9 <- 0.01
  r10 <- 0.25
  # Build into matrix
  par.1[["R"]] <- diag(c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10),n,n)
  
  # Coefficient matrix on noise
  par.1[["G"]] <- diag(1,n,n)
  
  # Priors
  priors <- list()
  
  # (mu, sigma) parameters for priors on A
  mu_A1 <- rep(0,n)
  mu_A2 <- rep(0,n)
  mu_A <- cbind(mu_A1,mu_A2)
  
  sig_A1 <- rep(0,n)
  sig_A2 <- rep(0,n)
  # Sig_A <- cbind(1/sig_A1,1/sig_A2)
  Sig_A <- matrix(0,n,m)
  
  priors[['mu_A']] <- mu_A
  priors[['Sig_A']] <- Sig_A
  
  # (mu, sigma) parameters for priors on B
  mu_B1 <- 0
  mu_B2 <- 0
  
  sig2_B1 <- 0
  sig2_B2 <- 0
  
  mu_B <- diag(c(mu_B1,mu_B2),m)
  #Sig_B <- diag(c(1/sig2_B1,1/sig2_B2),m)
  Sig_B <- matrix(0,m,m)
  
  priors[['mu_B']] <- mu_B
  priors[['Sig_B']] <- Sig_B
  
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
  #priors[['alpha_r']] <- vec_r_al
  #priors[['beta_r']] <- vec_r_bet
  
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
  Sig_c <- diag(c(sig2_c1,sig2_c2,sig2_c3,sig2_c4,sig2_c5,sig2_c6,sig2_c7,sig2_c8,sig2_c9,sig2_c10),m)
  #Sig_inv_c <- solve(Sig_c)
  Sig_inv_c <- matrix(0,n,n)
  priors[['mu_c']] <- mu_c
  priors[['Sig_inv_c']] <- Sig_inv_c
  
  # Priors on d
  # Normal distributions (mu_d,sig2_d)
  mu_d1 <- 0
  mu_d2 <- 0
  mu_d <- matrix(c(mu_d1,mu_d2),m,1)
  
  sig2_d1 <- 0
  sig2_d2 <- 0
  Sig_d <- diag(c(sig2_d1,sig2_d2),m)
  #Sig_inv_d <- solve(Sig_d)
  Sig_inv_d <- matrix(0,m,m)
  priors[['mu_d']] <- mu_d
  priors[['Sig_inv_d']] <- Sig_inv_d
  
  # Miscellaneous settings
  settings<-list()
  settings[['tinit']]=1
  
  # Storage of model dims and parameters
  mod[['dims']]<-dims
  mod[['par']]<-par.1
  mod[['settings']]<-settings
  mod[['priors']]<-priors
  return(mod)
}

# Kalman filter section
run_kf <- function(model,lag1=TRUE){
  # Extract data and relevant dimensions from model object
  datamat <- model$data
  n <- model$dims$n
  m <- model$dims$m
  TT <- model$dims$TT
  h <- model$dims$h
  g <- model$dims$g
  par.1 <- model$par
  
  # Check for NA values and replace in the yt matrix accordingly
  YM <- matrix(as.numeric(!is.na(datamat)),n,TT)
  yt <- datamat
  y <-datamat
  # For use in expectation calculations
  y[!YM] <- 0
  # For input to KFAS
  yt[!YM] <- as.numeric(NA)
  
  # Convert orientation to rows
  yt <- t(yt)
  
  # Construction of Zt matrix - assumes hidden state has an additional dimension with value 1 for offset
  # In AUD notation this is an extension of the A matrix to include the constant "c" vector
  Zt <- cbind(par.1$A,par.1$c)

  # Lag-1 version (needs to have 0's in matrix to account for state augmentation with previous state)
  stack.Zt <- matrix(0,n,2*(m+1))
  stack.Zt[1:n, 1:(m+1)] <- Zt
  
  # Construction of Tt matrix - similar to above, incorporation of offset term into matrix 
  # Equivalent of B in AUD notation with constant vector 'd' included
  Tt <- cbind(rbind(par.1$B,matrix(0,1,m)),matrix(c(par.1$d,1),m+1,1))
  
  # Lag-1 version (block identity matrix propogates the augmented hidden state correctly)
  stack.Tt <- matrix(0,2*(m+1),2*(m+1))
  stack.Tt[1:(m + 1), 1:(m + 1)] <- Tt
  stack.Tt[(m+2):(2*m+2),1:(m+1)] <- diag(1,m+1)
  
  # Build the Ht matrix
  # Equivalent of (G^T R G) in AUD notation -- observation process noise with noise coefficient matrix
  # Note that AUD models typically use identity G matrices, making this effectively the R matrix
  # No need to create a lag-1 version since the augmented state dimension doesn't impact dimension of additive terms in obs. eqn.
  Ht <- tcrossprod(par.1$G %*% par.1$R, par.1$G)
  
  # Build the Qt matrix
  # Same notation in KFAS as in AUD, Q is the transition equation noise variance-covariance
  Qt <- par.1$Q
  # Lag-1 version (need to add 0's to account for augmented hidden state)
  # Dimensions get sorted out when multiplying by KFAS Rt matrix (to account for constant state value)
  stack.Qt <- matrix(0,2*h,2*h)
  stack.Qt[1:h,1:h]<-Qt
  
  # Build the Rt matrix
  # Equivalent of H matrix in AUD notation (coefficient matrix on state noise)
  # m+1 needed to account for constant latent state
  Rt <- diag(1,m+1,h)
  Rt[1:m, 1:h] <- par.1$H
  # Lag-1 version (need to add 0's to account for augmented hidden state)
  stack.Rt <- matrix(0, 2 * (m + 1), 2 * h)
  stack.Rt[1:(m + 1), 1:h] <- Rt
  
  # Start with tinit=1
  x10 <- par.1$x0
  x00 <- matrix(10,m,1) # dummy?
  V10 <- par.1$V0
  V00 <- diag(0,m) # dummy? 
  
  # a1 vector is mean of the initial state
  a1 <- rbind(x10,1)
  stack.a1 <- rbind(x10,1,x00,1)
  
  # Set up the initial state covariance matrices (P1)
  # Rest of the variance covariance is determined by P1 and system equations
  # KFAS accepts diffuse and non-diffuse inputs, inf corresponds to part of diffuse handling
  # TODO read into diffuse handling
  # Current model structure considers no diffuse part
  P1inf <- matrix(0, m + 1, m + 1)
  stack.P1inf <- matrix(0, 2 * (m + 1), 2 * (m + 1))
  # Non-diffuse part of P1
  P1 <- matrix(0, m + 1, m + 1)
  P1[1:m, 1:m] <- V10
  # Following notation from MARSS for lag-1 variance covariance
  # P1 is the var-cov matrix for the stacked x1,x0
  # it is matrix(c(V10,B*V00,V00*t(B),V00),2,2)
  # x1=B*x0+U+w; in the var-cov mat looks like
  # E[x1*t(x1)] E[(Bx0+U+w1)*t(x0)] -   E[x1]*E[x1]     E[(Bx0+U+w1)]E[t(x0)]
  # E[x0*t(Bx0+U+w1)] E[x0*t(x0)]       E[(Bx0+U+w1)]E[t(x0)] - E[x0]E[t(x0)]
  # Note dependence on choice of tinit=0 vs. 1
  stack.P1 <- matrix(0, 2 * (m + 1), 2 * (m + 1))
  stack.P1[1:m, 1:m] <- V10
  stack.P1[(m + 2):(2 * m + 1), (m + 2):(2 * m + 1)] <- V00
  stack.P1[1:m, (m + 2):(2 * m + 1)] <- par.1$B %*% V00
  stack.P1[(m + 2):(2 * m + 1), 1:m] <- tcrossprod(V00, par.1$B)
  
  # Create the model
  # Pass different matrices depending on lag-1 or not
  if(lag1){
    kfas.model <- SSModel(yt ~ -1 + SSMcustom(Z = stack.Zt, T = stack.Tt, R = stack.Rt, Q = stack.Qt, a1 = stack.a1, P1 = stack.P1, P1inf = stack.P1inf), H = Ht)
  }else{
    kfas.model <- SSModel(yt ~ -1 + SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), H = Ht)
  }
  kfas.model$tol <- 0

  # For n>1
  diag.R <- unname(par.1$R)[1 + 0:(n - 1) * (n + 1)]

  # Actual call to the filter and smoother
  if (any(diag.R == 0)) { # because KFAS 1.0.4 added a warning message
    ks.out <- suppressWarnings(KFS(kfas.model, simplify = FALSE))
  } else {
    ks.out <- KFS(kfas.model, simplify = FALSE)
  }
  ks.out$a <- t(ks.out$a)
  ks.out$alphahat <- t(ks.out$alphahat)
  ks.out$att <- t(ks.out$att)
  model[["ks_raw"]]<-ks.out
  
  
  # Processed KF quantities directly available from KFS
  proc <-list()
  
  # VtT is the state variance covariance matrix based on the complete dataset y1:yT
  VtT <-ks.out$V[1:m, 1:m, , drop=FALSE]
  
  Vtt1T <- ks.out$V[1:m, (m + 2):(2 * m + 1), , drop = FALSE]
  if (model$settings$tinit == 1) Vtt1T[, , 1] <- matrix(NA, m, m)
  
  # Initial state related estimates
  x01T <- ks.out$alphahat[1:m, 1, drop = FALSE]
  V10T <- matrix(VtT[, , 1], m, m)
  # Note, currently only use the tinit=1 setting
  if (model$settings$tinit == 1) {
    x00 <- matrix(NA, m, 1)
    V00 <- matrix(NA, m, m)
    x0T <- x01T
    V0T <- V10T
  }
  # xtT is the predicted state values based on the complete dataset y1:yT
  xtT<-ks.out$alphahat[1:m, , drop = FALSE]
  
  # Expectation of observations given data (this accounts for missing data)
  hatyt <- matrix(0, n, TT)
  hatOt <- array(0, dim = c(n, n, TT))
  hatyxt <- hatyxttp <- array(0, dim = c(n, m, TT))
  hatxtp <- cbind(xtT[, 2:TT, drop = FALSE], NA)
  hatVtpt <- array(NA, dim = dim(Vtt1T))
  hatVtpt[, , 1:(TT - 1)] <- Vtt1T[, , 2:TT, drop = FALSE]
  I.n <- diag(1, n)
  for (t in 1:TT){
    # No missing data case (unlikely)
    if (all(YM[, t] == 1)) {
      hatyt[, t] <- y[, t, drop = FALSE]
      hatOt[, , t] <- tcrossprod(hatyt[, t, drop = FALSE])
      hatyxt[, , t] <- tcrossprod(hatyt[, t, drop = FALSE], xtT[, t, drop = FALSE])
      hatyxttp[, , t] <- tcrossprod(hatyt[, t, drop = FALSE], hatxtp[, t, drop = FALSE])
    } else { # Missing data case
      # IMPORTANT: We assume that R is diagonal (different calculations required if R is not diagonal)
      I.2 <- I.r <- I.n
      # Recal YM is boolean matrix where !na = TRUE, na = FALSE
      I.2[YM[, t] == 1, ] <- 0 # zero out rows with actual observations
      I.r[YM[, t] == 0, ] <- 0 # zero out rows with missing observations
      Delta.r <- I.n - I.r
      
      # See MARSS paper for derivations
      hatyt[, t] <- y[, t, drop = FALSE] - Delta.r %*% (y[, t, drop = FALSE] - par.1$A %*% xtT[, t, drop = FALSE] - par.1$c)
      # Grabs the entries of A associated with missing observations for this time step t(Delta.r * A)
      t.DA <- matrix(Delta.r %*% par.1$A, m, n, byrow = TRUE)
      hatOt[, , t] <- I.2 %*% (Delta.r %*% par.1$R + Delta.r %*% par.1$A %*% VtT[, , t] %*% t.DA) %*% I.2 + tcrossprod(hatyt[, t, drop = FALSE])
      hatyxt[, , t] <- tcrossprod(hatyt[, t, drop = FALSE], xtT[, t, drop = FALSE]) + Delta.r %*% par.1$A %*% VtT[, , t]
      hatyxttp[, , t] <- tcrossprod(hatyt[, t, drop = FALSE], hatxtp[, t, drop = FALSE]) + Delta.r %*% tcrossprod(par.1$A, hatVtpt[, , t])
    }
  }
  
  # Attach quantities
  proc[['VtT']] <- VtT
  proc[['Vtt1T']] <- Vtt1T
  proc[['xtT']] <- xtT
  proc[['ytT']] <- hatyt
  proc[['yytT']] <- hatOt
  proc[['yxtT']] <- hatyxt
  
  proc[['loglik']] <- ks.out$logLik
  
  model[['kf_proc']]<-proc
  # # # # # # # # # # # #
  # # Debugging block   #
  # # # # # # # # # # # #
  
  # # Params
  # cat(n,m,TT)
  # print(par.1)
  
  # # Z
  # print("Zt")
  # print(Zt)
  # print("stack.Zt")
  # print(stack.Zt)
  # 
  # # T
  # print("Tt")
  # print(Tt)
  # print("stack.Tt")
  # print(stack.Tt)
  # 
  # # H
  # print("H")
  # print(Ht)
  # 
  # # Q
  # print("Q")
  # print(Qt)
  # print("stack.Qt")
  # print(stack.Qt)
  # 
  # # R
  # print("Rt")
  # print(Rt)
  # print("stack.Rt")
  # print(stack.Rt)
  # 
  # # Effective state noise variance covariance matrix
  # print(stack.Rt%*% stack.Qt %*% t(stack.Rt))
  
  # # alpha
  # print("a1")
  # print(a1)
  # print("stack.a1")
  # print(stack.a1)
  
  # # P1
  # print("p1_inf")
  # print(P1inf)
  # print("stack.P1inf")
  # print(stack.P1inf)
  # print("P1")
  # print(P1)
  # print("stack.P1")
  # print(stack.P1)
  
  
  return(model)
}
run_em <- function(model){
  
  # Based on the tinit chosen, choose the relevant string
  kf.x0 <- ifelse(model$settings$tinit == 1, "x10", "x00")
  
  # Pull in the relevant model/data as local variables
  y <- model$data
  
  # Relevant initializations
  kf <- NULL
  loglik_old <- -Inf
  
  # Run the EM loop
  # EM parameters, can later be shifted as model settings
  maxit <- 5000
  tol <- 0.000001
  converged <- 0
  conv_check_min_iter <- 15
  error_status<-0
  tic('runtime')
  for(iter in 1:(maxit+1)){
    #print(iter)
    # Expectation/E Step
    # Store the last iteration's processed KF outputs
    kf.last <- kf
    # Pull in the current iteration's KF output
    model <- run_kf(model)
    kf <- model$kf_proc
    # Get specific expectations and likelihood
    xtT <- kf$xtT
    VtT <- kf$VtT
    Vtt1T <- kf$Vtt1T
    ytT <- kf$ytT
    yytT <- kf$yytT
    yxtT <- kf$yxtT
    
    loglik<-kf$loglik
    #print(loglik)
    if(converged){
      print(loglik)
      break
    }
    # Model parameters
    n <- model$dims$n
    m <- model$dims$m
    TT <- model$dims$TT
    
    # Convergence check
    cvg <- loglik-loglik_old 
    # Check for a drop in log likelihood
    if(cvg < -sqrt(.Machine$double.eps)){
      error_status<-1
    }
    
    if(error_status==1){
      cat("Error. Log-likelihood dropped on iter:",iter)
      break
    }
    if(iter>conv_check_min_iter){
      if(cvg>=0 && cvg<tol){
        cat("EM satisfied convergence criterion\n")
        converged<-TRUE
      }
    }
    if(iter>maxit){
      iter <- maxit
      break
    }
    loglik_old <- loglik
      
    # Maximization/M Step
    # Initial parameter values
    inits <- list()
    par.1 <-model$par
    A <- model$par$A
    B <- model$par$B
    R <- model$par$R
    c <- model$par$c
    d <- model$par$d
    # Assumes initial model structure
    # A matrix - fully estimated
    init_A <- par.1$A
    inits[['A']] <-init_A
    # B matrix - diagonal is estimated
    init_B <- to_vec(for(i in 1:dim(par.1$B)[1]) par.1$B[i,i])
    inits[['B']] <-init_B
    # c vector - fully estimated
    init_c <- par.1$c
    inits[['c']]
    # d vector - fully estimated
    init_d <-par.1$d
    inits[['d']]<init_d
    
    # A update (fully estimated), note assumption that R is diagonal
    Sig_A <- model$priors$Sig_A
    mu_A <- model$priors$mu_A
    # lhs_sum_term <- 0
    # rhs_sum_term <- 0
    Amat <- t(Sig_A) %*% R
    Bmat <- 0
    Cmat <- 0
    for(i in 1:TT){
      Bmat <- Bmat + VtT[,,i] + tcrossprod(xtT[,i,drop=FALSE])
      Cmat <- Cmat + t(yxtT[,,i])-(xtT[,i,drop=FALSE] %*% t(c))
    }
    A_new <- matrix(0,n,m)
    for(i in 1:n){
      a_i <- Amat[,i]
      diag_a_i <- diag(a_i,m)
      c_i <- Cmat[,i,drop=FALSE]
      a_new_i <- solve(diag_a_i + Bmat) %*% c_i
      A_new[i,]<- a_new_i
    }
    model$par$A <- A_new
    # for(i in 1:TT){
    #   # XtXt^T coefficient term
    #   lhs_sum_term <- lhs_sum_term + t(VtT[,,i]+tcrossprod(xtT[,i,drop=FALSE]))
    #   # Terms without A as a coefficient
    #   rhs_sum_term <- rhs_sum_term + t(yxtT[,,i]) - xtT[,i,drop=FALSE] %*% t(c)
    # }
    # A_new <- matrix(0,n,m)
    # for(i in 1:n){
    #   for(j in 1:m){
    #     A_new[j,i] <- (rhs_sum_term[j,i]+t(Sig_A)[j,i]*t(mu_A)[j,i]*R[i,i])/()
    #   }
    # }
    
    # B update (relies on diagonal structure)
    Sig_B <- model$priors$Sig_B
    mu_B <- model$priors$mu_B
    lhs_sum_term <- 0
    rhs_sum_term <- 0
    
    for(i in 2:TT){
      lhs_sum_term <- lhs_sum_term + (VtT[,,i-1]+tcrossprod(xtT[,i-1,drop=FALSE]))
      rhs_sum_term <- rhs_sum_term + (-d %*% t(xtT[,i-1,drop=FALSE])) + Vtt1T[, , i] + tcrossprod(xtT[,i,drop=FALSE], xtT[,i-1,drop=FALSE])
    }
    B_raw <- solve(Sig_B + lhs_sum_term) %*% (rhs_sum_term + Sig_B %*% mu_B)
    B_new <- diag(diag(B_raw),m,m)
    # Need to add a check here for values > 1
    #print(B_new)
    model$par$B <- B_new
    
    # R update (assumes diagonal structure)
    alpha_tilde <- model$priors$r_alpha_tilde
    beta_tilde <- model$priors$r_beta_tilde
    # alpha_r <- model$priors$alpha_r
    # beta_r <- model$priors$beta_r
    sum_term <- 0
    for(i in 1:TT){
      sum_term_it <- (yytT[,,i]
                  - (yxtT[,,i] %*% t(A))
                  - (A %*% t(yxtT[,,i]))
                  - (ytT[,i,drop=FALSE] %*% t(c))
                  - (c %*% t(ytT[,i,drop=FALSE]))
                  + (A %*% (VtT[,,i]+tcrossprod(xtT[,i,drop=FALSE])) %*% t(A))
                  + (A %*% xtT[,i,drop=FALSE] %*% t(c))
                  + (c %*% t(xtT[,i,drop=FALSE]) %*% t(A))
                  + (c %*% t(c))
                  )
      sum_term <- sum_term + sum_term_it
    }
    # R_new <- matrix(0,n,n)
    # for(i in 1:n){
    #   R_new[i,i] <- 1/(TT + 2*alpha_r[i]+2)*(sum_term[i,i]+2*beta_r[i])
    # }
    R_raw <- solve(2*alpha_tilde+TT*diag(1,n))%*%(sum_term+2*beta_tilde)
    R_new <- diag(diag(R_raw),n,n)
    # Variance assignment (fixed lapse variance case)
    R_new[n,n]<-0.25
    negativity_check <- any(R_new<=-1*10^-10)
    if(negativity_check){
      print("Negative variance")
      print(R_new)
      break
    }
    #print(R_new)
    #(diag(R_raw))
    model$par$R <- R_new
    
    # c update
    mu_c <- model$priors$mu_c
    Sig_inv_c <- model$priors$Sig_inv_c
    
    sum_term<-0
    for(i in 1:TT){
      sum_term_it <- ytT[,i,drop=FALSE]-A%*%xtT[,i,drop=FALSE]
      sum_term<-sum_term +sum_term_it
    }
    c_new <- solve(TT*diag(1,n)+R%*%Sig_inv_c) %*% (sum_term + R%*%Sig_inv_c%*%mu_c)
    model$par$c <- c_new
    
    # d update
    mu_d <- model$priors$mu_d
    Sig_inv_d <- model$priors$Sig_inv_d
    
    sum_term<-0
    x0 <- xtT[,1,drop=FALSE]
    # Pulled in negative sign & distributed 1/2 since Q=I
    t1_term <- x0-(B%*%x0)
    sum_term<- sum_term + t1_term
    for(i in 2:TT){
      sum_term_it <- xtT[,i,drop=FALSE]-B%*%xtT[,i-1,drop=FALSE]
      sum_term<-sum_term +sum_term_it
    }
    d_new <- solve(TT*diag(1,m)+Sig_inv_d) %*% (sum_term + Sig_inv_d%*%mu_d)
    model$par$d <- d_new
    #print(d_new)
  }
  toc()
  return(model)
}
# Data import
# Import the data
path <- '/Users/eric/repos/aud/data/subject_marss_60_90_covar.csv'
data <- read.csv(path)
subid_list <- unique(data$subid)

debug_id <- 54
current_data <- data[data$subid==debug_id,]
datamat <- t(as.matrix(current_data))
model <- get_dummy_model(datamat)
model <- run_kf(model)
model <- run_em(model)

