# Initialize the parameter values (uses a fixed value obtained as a population average of early MARSS MLE fits)
# Future work to incorporate multiple random initialization to evaluate if this leads to materially better fits
init_par <- function(dims){
  par.1 <- list()
  n <- dims$n
  m <- dims$m
  g <- dims$g
  h <- dims$h
  
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
  
  return(par.1)
}

# Take the MLE fits from other individuals and form prior distributions for the new subject
make_priors <- function(model,prior_dist,id,included_format=FALSE,included_list=c(),bound=TRUE){
  # Note: function was originally written for leave-one-out structure, so all other subject fits could be used
  # included_format structure was added to accommodate repeated cross validation, where a specific subset was used for fitting
  
  # Set up list for storing priors
  priors <- list()
  # Remove the current individual from the prior distribution
  if(included_format){
    prior_dist <- prior_dist[prior_dist$id %in% included_list,]
  } else {
    prior_dist <- prior_dist[prior_dist$id != id,]
  }
  print(dim(prior_dist))
  # Dimension variables
  m <- model$dims$m
  n <- model$dims$n
  
  # Preallocate matrices
  mu_A <- matrix(0,n,m)
  sig_A <- matrix(0,n,m)
  comb_A <- matrix(0,n,m)
  mu_B <- matrix(0,m,m)
  sig_B <- matrix(0,m,m)
  mu_c <- matrix(0,n,1)
  sig_c <- matrix(0,n,n)
  mu_d <- matrix(0,m,1)
  sig_d <- matrix(0,m,m)
  al_r <- matrix(0,n,n)
  bet_r <- matrix(0,n,n)
  
  # Loop over the model dimensions
  for(i in 1:n){
    for(j in 1:m){
      # Normal fitting for A values (Z in MARSS notation)
      iter_label <- paste('z',i,j,sep='')
      iter_dist <- prior_dist[[iter_label]]
      est <- fitdist(iter_dist,'norm')
      mu_A[i,j] <- est$estimate[1]
      sig_A[i,j] <- 1/est$estimate[2]^2
      comb_A[i,j] <- mu_A[i,j]*sig_A[i,j]
      
      # Run d and B fits during the j loop on the first iteration
      if(i==1){
        # Normal fitting for d values (u in MARSS notation)
        iter_label <- paste('u',j,sep='')
        iter_dist <- prior_dist[[iter_label]]
        # Used in the manuscript to address a few outliers that strongly affected the prior fitting
        if(bound){
          d_lb <- -5
          d_ub <- 5
          iter_dist <- iter_dist[iter_dist > d_lb]
          iter_dist <- iter_dist[iter_dist < d_ub]
        }
        est <- fitdist(iter_dist,'norm')
        mu_d[j] <- est$estimate[1]
        sig_d[j,j] <- 1/est$estimate[2]^2
        
        # Truncated normal fitting for B (same in MARSS notation)
        iter_label <- paste('b',j,j,sep='')
        iter_dist <- prior_dist[[iter_label]]
        #print(iter_dist)
        eps <- 1*10^(-10)
        lb <- -1+eps
        ub <- 1-eps
        iter_dist[iter_dist<lb]<-lb
        iter_dist[iter_dist>ub]<-ub
        
        if(j==1){
          est <- fitdist(iter_dist,dtruncnorm,method='mle',start=list(mean=mean(iter_dist),sd=0.5),fix.arg=list(a=lb-eps,b=ub+eps))
        } else{
          est <- fitdist(iter_dist,dtruncnorm,method='mle',start=list(mean=max(iter_dist),sd=10),fix.arg=list(a=lb-eps,b=ub+eps))
        }
        
        mu_B[j,j]=est$estimate[1]
        sig_B[j,j] <- 1/est$estimate[2]^2
      }
    }
    # Normal fitting for c values (a in MARSS notation)
    iter_label <- paste('a',i,sep='')
    iter_dist <- prior_dist[[iter_label]]
    est <- fitdist(iter_dist,'norm')
    mu_c[i] <- est$estimate[1]
    sig_c[i,i] <- 1/est$estimate[2]^2
    
    # Inverse gamma fitting for R (same in MARSS notation)
    # i<n needed for case where lapse variance is not fitted
    if(i<n){
      iter_label <- paste('r',i,sep='')
      #print(iter_label)
      iter_dist <- prior_dist[[iter_label]]
      #print(length(iter_dist))
      lb <- 0.001
      #iter_dist <- iter_dist[iter_dist > lb]
      iter_dist[iter_dist<lb]<-lb
      #print(length(iter_dist))
      est <- fitdist(iter_dist,dinvgamma,method='mle',start=list(shape=2,rate=.1),lower=0)
      al_r[i,i] <- est$estimate[1]
      bet_r[i,i] <- est$estimate[2]
    }
  }
  
  # Store priors
  priors[['mu_A']]<-mu_A
  priors[['sig_A']]<-sig_A
  priors[['comb_A']]<-comb_A
  
  priors[['mu_B']]<-mu_B
  priors[['sig_B']]<-sig_B
  
  priors[['mu_c']]<-mu_c
  priors[['sig_c']]<-sig_c
  
  priors[['mu_d']]<-mu_d
  priors[['sig_d']]<-sig_d
  
  priors[['r_alpha_tilde']]<- al_r + diag(c(rep(1,n-1),0),n,n)
  priors[['r_beta_tilde']] <- bet_r
  return(priors)
}

# Model fits can be degenerate based on which latent state coefficients are attached to
# This uses a convention that the first column of A should have the greatest norm
# Columns are switched as needed, this function determines if a switch is necessary
flip_check <- function(model){
  A<-model$par$A
  a1<-A[,1]
  a2<-A[,2]
  na1 <- norm(a1,type='2')
  na2 <- norm(a2,type='2')
  return(na2>na1)
}

# Zero out the priors (i.e., set priors such that they have no impact on parameter updates, leading to equivalent fit to MLE)
zero_priors <- function(model,param_list){
  # Priors
  priors <- model$priors
  
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
  
  if('A' %in% param_list){
    priors[['mu_A']] <- mu_A
    priors[['sig_A']] <- sig_A
    priors[['comb_A']] <- comb_A
  }
  
  # (mu, sigma) parameters for priors on B
  mu_B1 <- 0
  mu_B2 <- 0
  
  sig2_B1 <- 0
  sig2_B2 <- 0
  
  mu_B <- diag(c(mu_B1,mu_B2),m)
  sig_B <- matrix(0,m,m)
  
  if('B' %in% param_list){
    priors[['mu_B']] <- mu_B
    priors[['sig_B']] <- sig_B
  }
  
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
  if('R' %in% param_list){
    priors[['r_alpha_tilde']]<-r_alpha_tilde
    priors[['r_beta_tilde']]<-r_beta_tilde
  }
  
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
  if('c' %in% param_list){
    priors[['mu_c']] <- mu_c
    priors[['sig_c']] <- sig_c
  }

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
  if('d' %in% param_list){
    priors[['mu_d']] <- mu_d
    priors[['sig_d']] <- sig_d
  }
  
  model[['priors']] <- priors
  return(model)
}

# For MARSS models only
check_flip <- function(coefmat){
  z1 <- coefmat$Z[,1]
  z2 <- coefmat$Z[,2]
  nz1 <-norm(z1,type='2')
  nz2 <-norm(z2,type='2')
  switch <- (nz2>nz1)
  return(switch)
}

# Calculate the component of the total likelihood associated with the prior parameters
calc_prior_loglik<-function(mod){
  ploglik <- 0
  
  # Terms that don't include estimated parameters are not included
  # They are the same each iteration and have no impact on convergence
  # A terms
  for(i in 1:mod$dims$n){
    for(j in 1:mod$dims$m){
      param <- mod$par$A[i,j]
      prior_mu <- mod$priors$mu_A[i,j]
      prior_sig2 <- mod$priors$sig_A[i,j] # this is 1/var
      term_contribution <- -1/2*(param-prior_mu)^2*prior_sig2
      #print(term_contribution)
      ploglik <- ploglik + term_contribution
    }
  }
  
  # c terms
  for(i in 1:mod$dims$n){
    param <- mod$par$c[i]
    prior_mu <- mod$priors$mu_c[i]
    prior_sig2 <- mod$priors$sig_c[i]
    term_contribution <- -1/2*((param-prior_mu)^2)*prior_sig2
    ploglik <- ploglik + term_contribution
  }
  
  # B terms
  for(i in 1:mod$dims$m){
    param <- mod$par$B[i,i]
    prior_mu <- mod$priors$mu_B[i,i]
    prior_sig2 <- mod$priors$sig_B[i,i]
    term_contribution <- -1/2*(param-prior_mu)^2*prior_sig2
    ploglik <- ploglik + term_contribution
  }
  
  # d terms
  for(i in 1:mod$dims$m){
    param <- mod$par$d[i]
    prior_mu <- mod$priors$mu_d[i]
    prior_sig2 <- mod$priors$sig_d[i]
    term_contribution <- -1/2*(param-prior_mu)^2*prior_sig2
    ploglik <- ploglik + term_contribution
  }
  
  # R terms
  # NOTE, assumes lapse not fit
  for(i in 1:(n-1)){
    param<- mod$par$R[i,i]
    prior_a_tilde <- mod$priors$r_alpha_tilde[i,i]
    prior_b_tilde <- mod$priors$r_beta_tilde[i,i]
    term_contribution <- -prior_a_tilde*log(param) - prior_b_tilde/param
    ploglik <- ploglik + term_contribution
  }
  #print(loglik)
  return(ploglik)
}

# Helper function for formatting outputs
get_out_frame <- function(type,subid,horizon,day,mod_obj,pred_0,act_0){
  output <- data.frame(fit_type=type,
                       subid=subid,
                       train_horizon=horizon,
                       day = day,
                       converged=mod_obj$converged,
                       fit_error=mod_obj$error,
                       pred_0=pred_0,
                       act_0=act_0)
  return(output)
}

# Function for forecasting the next step
fr_step <-function(B,drep,x,w){
  xt1 <- B %*% x + drep + w
  return(xt1)
}
# Function for returning the clipped lapse probability
# Future work to add noise here to determine if bias concern is legitimate
lapse_step <- function(A10,c10rep,x){
  #lapse <- A10 %*% x + c10rep + v
  lapse <- A10 %*% x + c10rep
  return(lapse)
}

# Function for clipping predictions for [epsilon, 1-epsilon]
clip <- function(lapse){
  eps <- 1*10^(-12)
  if(lapse < eps){
    return(eps)
  }
  if(lapse > 1-eps){
    return(1-eps)
  }
  return(lapse)
}

# Function for returning a vector of lapse probabilities for the given window
calc_lapse_window <- function(mat){
  prod_val <- apply(mat,2,prod)
  p_window_lapse <- 1-prod_val
  return(p_window_lapse)
}

# Function for generating window predictions using simulated trajectories
get_window_preds<-function(model_pars,input_data,meanTT,varTT,type,TT,Tfinal,window_length){
  # Number of samples for MC averaging
  n<-20000
  
  # Create lists for storing forecast results
  fr_x <- list()
  fr_lapse <- list()
  
  # Get the model coefficients
  if(type=='MAP_mle'||type=='MAP_marss' || type=='MLE'){
    B <- model_pars$B
    d <- model_pars$d
    A <- model_pars$A
    c <- model_pars$c
  } else {
    B <- model_pars$B
    d <- model_pars$U
    A <- model_pars$Z
    c <- model_pars$A
  }
  c10 <- c[10]
  A10 <- A[10,]
  
  # Tile vectors for batch processing of samples
  drep <- t(matrix(d,nrow=n,ncol=length(d), byrow=TRUE))
  c10rep <- t(matrix(c10, nrow=n, ncol=length(c10),byrow=TRUE))
  
  # Generate a set of initiate state estimates
  xt0 <- t(rmvnorm(n,mean=meanTT,sigma=varTT))
  
  # Store the initial states as the first entry of the fr_x list
  fr_x[[paste(0)]] <- xt0
  # Store same day prediction values (complement of probability of no lapse)
  fr_lapse[[paste(0)]] <- 1-vapply(lapse_step(A10,c10rep,xt0),clip,numeric(1))
  
  # Loop over desired prediction horizon
  for(i in 1:window_length){
    if((i+TT)<=Tfinal){
      # Use the previous day's x value (each trajectory is independent)
      x_iter <- fr_x[[paste(i-1)]]
      # Generate identity noise for every trajectory
      w_iter <- t(rmvnorm(n,rep(0,2),sigma=diag(2)))
      # Create current step's hidden state using the fr_step function (just a step forward in the transition equation)
      xt <- fr_step(B,drep,x_iter,w_iter)
      # Store the current step's hidden state
      fr_x[[paste(i)]] <- xt
      # Noiseless evaluation of the unclipped lapse probability
      lapset <- lapse_step(A10,c10rep,xt)
      # Clip the lapse probability
      clipped_lapset <- vapply(lapset,clip,numeric(1))
      # Store the probability of no lapse
      comp_lapset <- 1-clipped_lapset
      fr_lapse[[paste(i)]] <- comp_lapset
    } else {
      fr_x[[paste(i)]] <- rep(NA,n)
      fr_lapse[[paste(i)]] <- rep(NA,n)
    }
  }
  final <- list()
  lapse_mat <- do.call("rbind",fr_lapse)
  for(i in 1:length(fr_lapse)){
    mat <- matrix(lapse_mat[1:i,],nrow=i)
    final[[i]] <- calc_lapse_window(mat)
  }
  return(final)
}

# Function for running the kalman filter, built directly from MARSS code
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
  
  # Lag-1 version (block identity matrix propagates the augmented hidden state correctly)
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
  x00 <- matrix(0,m,1) # dummy?
  V10 <- par.1$V0
  V00 <- diag(1,m) # dummy? 
  
  # a1 vector is mean of the initial state
  a1 <- rbind(x10,1)
  stack.a1 <- rbind(x10,1,x00,1)
  
  # Set up the initial state covariance matrices (P1)
  # Rest of the variance covariance is determined by P1 and system equations
  # KFAS accepts diffuse and non-diffuse inputs, inf corresponds to part of diffuse handling
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
  #ks.out$alphahat <- t(ks.out$alphahat)
  ks.out$att <- t(ks.out$att)
  model[["ks_raw"]]<-ks.out
  
  # Processed KF quantities directly available from KFS
  proc <-list()
  
  # VtT is the state variance covariance matrix based on the complete dataset y1:yT
  VtT <-ks.out$V[1:m, 1:m, , drop=FALSE]
  
  Vtt1T <- ks.out$V[1:m, (m + 2):(2 * m + 1), , drop = FALSE]
  if (model$settings$tinit == 1) Vtt1T[, , 1] <- matrix(NA, m, m)
  
  xtT<-t(ks.out$alphahat)[1:m, , drop = FALSE]
  
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
      # Recall YM is boolean matrix where !na = TRUE, na = FALSE
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
  
  return(model)
} 

# Run the expectation-maximization loop until model fit reaches improvement tolerance
run_em <- function(model,maxit=10000,tol=0.0001,zero_priors_bool=FALSE,bound=FALSE){
  # Based on the tinit chosen, choose the relevant string
  kf.x0 <- ifelse(model$settings$tinit == 1, "x10", "x00")
  
  # Pull in the relevant model/data as local variables
  y <- model$data
  
  # Relevant initializations
  kf <- NULL
  loglik_old <- -Inf
  
  # Run the EM loop
  # EM parameters, can later be shifted as model settings
  converged <- 0
  conv_check_min_iter <- 15
  error_status<-0
  try_safe_update<-FALSE
  out_of_bound_b <- FALSE
  degen_r <- FALSE
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
    kf_loglik <- kf$loglik
    
    # Some numerical instability encountered for short/heavily missing trajectories and small convergence tolerances
    # Attempt to address by updating one parameter at a time (safe updates)
    # This is the 'correct' way to do so in general, but is slower and typically unnecessary
    # References to 'try_safe_update' refer to updates/reporting when in this mode
    if(try_safe_update){
      if(zero_priors_bool){
        prior_loglik <- 0
      } else {
        prior_loglik <- calc_prior_loglik(model)
      }
      update_loglik<-kf_loglik + prior_loglik
      print(paste(update_loglik,"A update"))
    }
    
    # Log likelihood reporting (accommodates both MLE and MAP based on whether priors are zeroed)
    if(zero_priors_bool){
      prior_loglik <- 0
    } else {
      prior_loglik <- calc_prior_loglik(model)
    }
    loglik<-kf_loglik + prior_loglik
    if((iter %% 500) == 0){
      cat(paste("Iter ", iter," loglik: ",loglik,"\n",sep=''))
    }
    # Convergence reporting
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
    if(iter > 2){
      if(cvg < -sqrt(.Machine$double.eps)){
        if(!try_safe_update){
          print(paste("Warning: log likelihood dropped by ",cvg," on iter: ",iter,". Trying safe updates.",sep=''))
          try_safe_update<-TRUE
        }else{
          print(paste("Fatal error: Log likelihood continues to drop. Stopping on iter ",iter,sep=''))
          error_status<-1
        }
      }
    }
    
    if(error_status==1){
      cat("Error. Log-likelihood dropped on iter:",iter)
      break
    }
    
    if(iter>conv_check_min_iter){
      if(cvg>=0 && cvg<tol){
        cat(paste("EM satisfied convergence criterion on iteration: ",iter,"\n",sep=''))
        converged<-TRUE
      }
    }
    if(iter>maxit){
      iter <- maxit
      break
    }
    loglik_old <- loglik
    
    # Maximization/M Step
    
    # Update order R,c,d,B,A
    # R update (assumes diagonal structure)
    alpha_tilde <- model$priors$r_alpha_tilde
    beta_tilde <- model$priors$r_beta_tilde
    sum_term <- 0
    for(i in 1:TT){
      sum_term_it <- (yytT[,,i]
                      - 2*(yxtT[,,i] %*% t(model$par$A))
                      - 2*(ytT[,i,drop=FALSE] %*% t(model$par$c))
                      + (model$par$A %*% (VtT[,,i]+tcrossprod(xtT[,i,drop=FALSE])) %*% t(model$par$A))
                      + 2*(model$par$A %*% xtT[,i,drop=FALSE] %*% t(model$par$c))
                      + (model$par$c %*% t(model$par$c))
      )
      sum_term <- sum_term + sum_term_it
    }
    sum_term <- diag(diag(sum_term),n,n)
    
    R_raw <- solve(2*alpha_tilde + TT*diag(1,n)) %*% (sum_term + 2*beta_tilde)
    R_new <- diag(diag(R_raw),n,n)
    # Variance assignment (fixed lapse variance case)
    R_new[n,n]<-0.25
    negativity_check <- any(R_new<=-1*10^(-10))
    if(negativity_check){
      print("Negative variance")
      print(R_new)
      break
    }
    # Some respondents had degenerate variances for specific questions (i.e., always answering with the same value)
    # Avoid numerical instability by enforcing a small tolerance that R cannot drop below
    R_tol <- 1*10^(-4)
    if(any(R_new <= R_tol)){
      for(i in 1:n){
        param <- R_new[i,i]
        if(param<R_tol){
          R_new[i,i] <- R_tol
        }
      }
      degen_r <- TRUE
    }
    model$par$R <- R_new
    
    if(try_safe_update){
      model <- run_kf(model)
      kf <- model$kf_proc
      # Get specific expectations and likelihood
      xtT <- kf$xtT
      VtT <- kf$VtT
      Vtt1T <- kf$Vtt1T
      ytT <- kf$ytT
      yytT <- kf$yytT
      yxtT <- kf$yxtT
      kf_loglik <- kf$loglik
      if(zero_priors_bool){
        prior_loglik <- 0
      } else {
        prior_loglik <- calc_prior_loglik(model)
      }
      print(paste(kf_loglik,prior_loglik))
      update_loglik<-kf_loglik + prior_loglik
      print(paste(update_loglik,"R update"))
    }
    
    # c update
    mu_c <- model$priors$mu_c
    sig_c <- model$priors$sig_c
    sum_term <- 0
    sum_term_it <- 0
    for(i in 1:TT){
      sum_term_it <- ytT[,i,drop=FALSE]-(model$par$A %*% xtT[,i,drop=FALSE])
      sum_term<-sum_term +sum_term_it
    }
    c_new <- ginv(TT*diag(1,n) + model$par$R %*% sig_c) %*% (sum_term + model$par$R %*% sig_c %*% mu_c)
    
    c_tol <- 1*10^(-12)
    for(i in 1:n){
      param <- c_new[i]
      if(abs(param) < c_tol){
        if(param>0){ c_new[i]<-c_tol} else {c_new[i]<- -c_tol}
        # if(try_safe_update){
        #   print("c update blocked")
        # }
      }
    }
    model$par$c <- c_new
    
    if(try_safe_update){
      model <- run_kf(model)
      kf <- model$kf_proc
      # Get specific expectations and likelihood
      xtT <- kf$xtT
      VtT <- kf$VtT
      Vtt1T <- kf$Vtt1T
      ytT <- kf$ytT
      yytT <- kf$yytT
      yxtT <- kf$yxtT
      kf_loglik <- kf$loglik
      if(zero_priors_bool){
        prior_loglik <- 0
      } else {
        prior_loglik <- calc_prior_loglik(model)
      }
      print(paste(kf_loglik,prior_loglik))
      update_loglik <- kf_loglik + prior_loglik
      print(paste(update_loglik,"c update"))
    }
    
    # d update
    mu_d <- model$priors$mu_d
    sig_d <- model$priors$sig_d
    
    sum_term<-0
    for(i in 2:TT){
      sum_term_it <- xtT[,i,drop=FALSE]-model$par$B %*% xtT[,i-1,drop=FALSE]
      sum_term<-sum_term +sum_term_it
    }
    d_new <- solve((TT-1)*diag(1,m) + sig_d) %*% (sum_term + sig_d %*% mu_d)
    model$par$d <- d_new
    
    if(try_safe_update){
      model <- run_kf(model)
      kf <- model$kf_proc
      # Get specific expectations and likelihood
      xtT <- kf$xtT
      VtT <- kf$VtT
      Vtt1T <- kf$Vtt1T
      ytT <- kf$ytT
      yytT <- kf$yytT
      yxtT <- kf$yxtT
      kf_loglik <- kf$loglik
      if(zero_priors_bool){
        prior_loglik <- 0
      } else {
        prior_loglik <- calc_prior_loglik(model)
      }
      print(paste(kf_loglik,prior_loglik))
      update_loglik<-kf_loglik + prior_loglik
      print(paste(update_loglik,"d update"))
    }
    
    # Check for need to flip columns to maintain convention (first column of A has largest norm)
    need_flip<-flip_check(model)
    
    # B update (relies on diagonal structure)
    sig_B <- model$priors$sig_B
    mu_B <- model$priors$mu_B
    lhs_sum_term <- 0
    rhs_sum_term <- 0
    
    for(i in 2:TT){
      lhs_sum_term <- lhs_sum_term + (VtT[,,i-1]+tcrossprod(xtT[,i-1,drop=FALSE]))
      rhs_sum_term <- rhs_sum_term + (-model$par$d %*% t(xtT[,i-1,drop=FALSE])) + Vtt1T[, , i] + tcrossprod(xtT[,i,drop=FALSE], xtT[,i-1,drop=FALSE])
    }
    lhs_sum_term <- diag(diag(lhs_sum_term),m,m)
    rhs_sum_term <- diag(diag(rhs_sum_term),m,m)
    B_raw <- (rhs_sum_term + sig_B %*% mu_B) %*% solve(sig_B + lhs_sum_term)
    B_new <- diag(diag(B_raw),m,m)
    
    # Enforce stability (maintain eigenvectors in unit circle)
    for(k in 1:m){
      if(B_new[k,k] > 1){
        B_new[k,k]<-1
        out_of_bound_b<-TRUE
      }
      if(B_new[k,k] < -1){
        B_new[k,k] <- -1
        out_of_bound_b<-TRUE
      }
    }

    model$par$B <- B_new
    if(try_safe_update){
      model <- run_kf(model)
      kf <- model$kf_proc
      # Get specific expectations and likelihood
      xtT <- kf$xtT
      VtT <- kf$VtT
      Vtt1T <- kf$Vtt1T
      ytT <- kf$ytT
      yytT <- kf$yytT
      yxtT <- kf$yxtT
      kf_loglik <- kf$loglik
      if(zero_priors_bool){
        prior_loglik <- 0
      } else {
        prior_loglik <- calc_prior_loglik(model)
      }
      print(paste(kf_loglik,prior_loglik))
      update_loglik<-kf_loglik + prior_loglik
      print(paste(update_loglik,"B update"))
    }
    
    # A update (fully estimated), note assumption that R is diagonal
    sig_A <- model$priors$sig_A
    mu_A <- model$priors$mu_A
    comb_A <- model$priors$comb_A
    
    Cmat <- 0
    Dmat <- model$par$R %*% sig_A
    Emat <- 0
    for(i in 1:TT){
      Cmat <- Cmat + VtT[,,i] + tcrossprod(xtT[,i,drop=FALSE])
      Emat <- Emat + yxtT[,,i]-(model$par$c %*% t(xtT[,i,drop=FALSE]))
    }

    Emat <- Emat + model$par$R %*% comb_A
    A_raw <- matrix(0,n,m)
    for(i in 1:n){
      d_i <- Dmat[i,]
      e_i <- Emat[i,]
      diag_d_i <- diag(d_i,m)
      a_raw_i <- e_i %*% ginv(Cmat + diag_d_i)
      A_raw[i,]<- a_raw_i
    }
    A_new <- A_raw
    a_tol<- 1*10^(-10)
    for(i in 1:n){
      for(j in 1:m){
        param <- A_new[i,j]
        if(abs(param)<=a_tol){
          if(param>0){
            A_new[i,j] <- a_tol
          } else {
            A_new[i,j]<- -a_tol
          }
        }
      }
    }
    if(try_safe_update){
      test <- Emat %*% ginv(Cmat)
    }
    old_A <- model$par$A
    model$par$A <- A_new
    
  }

  returns <- list()
  returns[['model']]<-model
  returns[['converged']]<-converged
  returns[['error']] <- error_status
  returns[['tol']] <- tol
  return(returns)
}

# Export coefficients of model as a frame
make_coef_frame <- function(model,need_flip,fit_lapse=FALSE){
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
  out<-cbind(data.frame(A),data.frame(c),data.frame(R),data.frame(B),data.frame(d))
  return(out)
}

# Assemble data frame of actual outcomes for comparison to predicted values
get_act_pred_frame <- function(window_length,TT, Tfinal,full_raw_datamat,mod_pars,raw_datamat,masked_datamat,meanTT,varTT,type){
  # Get the actual lapse outcomes for the different windows
  actuals_mat <- matrix(NA, nrow = 1, ncol = (window_length+1))
  actuals_df <- data.frame(actuals_mat)
  for(j in 0:(window_length)){
    if((TT+j)<=Tfinal){
      relevant_lapse_data <- full_raw_datamat[10,TT:(TT+j)]
      lapse_bool <- as.numeric(1 %in% relevant_lapse_data)
      val <- lapse_bool
    } else {
      val <- NA
    }
    actuals_df[(j+1)]<- val
  }
  # Make window predictions
  raw_window_preds <- get_window_preds(mod_pars,masked_datamat,meanTT,varTT,type,TT,Tfinal,window_length)
  proc_window_preds <- vapply(raw_window_preds,mean,numeric(1))
  probs_df <- data.frame(t(proc_window_preds))
  # Adjust column names
  colnames(probs_df)<-prob_col_names
  colnames(actuals_df)<-act_col_names
  return(cbind(probs_df,actuals_df))
}

# Process output frame specifically for the rolling 'data availability' processing scheme
get_rolling_out_frame <- function(type,subid,train_width,day,mod_obj,pred_0,act_0){
  output <- data.frame(fit_type=type,
                       subid=subid,
                       train_width= train_width,
                       day = day,
                       converged=mod_obj$converged,
                       fit_error=mod_obj$error,
                       pred_0=pred_0,
                       act_0=act_0)
  return(output)
}

