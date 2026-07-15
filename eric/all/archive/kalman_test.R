# Adapted from MARSS
library("KFAS")
# Set key constants

# Observed variable dimension
n <- 10
# Total time steps
TT <- 90
# Dimension of hidden states
m <- 2
# UNCLEAR - Dimension of hidden state variance covariance?
g1 <- 2
# UNCLEAR - Dimension of observed variance covariance?
h1 <- 10
# UNCLEAR - Dimension of initial state variance covariance?
l1 <- 2
# # # # # # # # # # # # # # # # # # #
# Parameter matrix                  #
# # # # # # # # # # # # # # # # # # #
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
par.1[["H"]] <- diag(1,m,m)

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

# INITIAL STATES

# B transpose (not really relevant for the diagonal case)
t.B <- matrix(par.1$B,m,m,byrow=TRUE)

# Data import
# Import the data
path <- '/Users/eric/repos/aud/data/subject_marss_60_90_covar.csv'
data <- read.csv(path)
subid_list <- unique(data$subid)

debug_id <- 268
current_data <- data[data$subid==debug_id,]
datamat <- t(as.matrix(current_data))
datamat <- datamat[-c(1,2,3,14,15,16,17,18,19,20),1:TT]

# Check for NA values and replace in the yt matrix accordingly
YM <- matrix(as.numeric(!is.na(datamat)),n,TT)
yt <- datamat
yt[!YM] <- as.numeric(NA)

# Convert orientation to rows
yt <- t(yt)

# Construction of At matrix - assumes hidden state has an additional dimension with value 1 for offset
# This is Z in the KFAS notation
At <- cbind(par.1$A,par.1$c)
# UNCLEAR - At is put into the first 'space' for it, not sure what the zeros are for yet
stack.At <- matrix(0,n,2*(m+1))
stack.At[1:n, 1:(m+1)] <- At

# Construction of Tt matrix - similar to above, incorporation of offset term into matrix 
# Equivalent of B in AUD notation
Tt <- cbind(rbind(par.1$B,matrix(0,1,m)),matrix(c(par.1$d,1),m+1,1))
stack.Tt <- matrix(0,2*(m+1),2*(m+1))
# UNCLEAR - not sure what the purpose of putting an identity matrix in the lower left block section of this is
stack.Tt[(m+2):(2*m+2),1:(m+1)] <- diag(1,m+1)

# Build the Ht matrix
# Equivalent of (G^T R G) in AUD notation -- observation process noise with noise coefficient matrix
Ht <- tcrossprod(par.1$G %*% par.1$R, par.1$G)

# Build the Qt matrix
# Same in both notations, the transition eqn noise variance-covariance
Qt <- par.1$Q
stack.Qt <- matrix(0,2*g1,2*g1)
# UNCLEAR - not sure what the purpose of putting an identity matrix in the upper left block section of this is
stack.Qt[1:g1,1:g1]<-Qt

# Build the Rt matrix
# Equivalent of H matrix in AUD notation (coefficient matrix on state noise)
Rt <- diag(1,m+1,g1)
Rt[1:m, 1:g1] <- par.1$H
stack.Rt <- matrix(0, 2 * (m + 1), 2 * g1)
stack.Rt[1:(m + 1), 1:g1] <- Rt

# Start with tinit=1
x10 <- par.1$x0
x00 <- matrix(10,m,1) # dummy?
V10 <- par.1$V0
V00 <- diag(0,m) # dummy? 
# UNCLEAR - need to read more about init=0 vs. 1 philosophy

# a1 vector is mean of the initial state
a1 <- rbind(x10,1)
stack.a1 <- rbind(x10,1,x00,1)

# Non-diffuse
P1inf <- matrix(0, m + 1, m + 1)
stack.P1inf <- matrix(0, 2 * (m + 1), 2 * (m + 1))
P1 <- matrix(0, m + 1, m + 1)
P1[1:m, 1:m] <- V10
# P1 is the var-cov matrix for the stacked x1,x0
# it is matrix(c(V10,B*V00,V00*t(B),V00),2,2)
# x1=B*x0+U+w; in the var-cov mat looks like
# E[x1*t(x1)] E[(Bx0+U+w1)*t(x0)] -   E[x1]*E[x1]     E[(Bx0+U+w1)]E[t(x0)]
# E[x0*t(Bx0+U+w1)] E[x0*t(x0)]       E[(Bx0+U+w1)]E[t(x0)] - E[x0]E[t(x0)]
stack.P1 <- matrix(0, 2 * (m + 1), 2 * (m + 1))
stack.P1[1:m, 1:m] <- V10
stack.P1[(m + 2):(2 * m + 1), (m + 2):(2 * m + 1)] <- V00
stack.P1[1:m, (m + 2):(2 * m + 1)] <- par.1$B %*% V00
stack.P1[(m + 2):(2 * m + 1), 1:m] <- tcrossprod(V00, par.1$B)

# Case for no lag 1 (sanity checking)
kfas.model <- SSModel(yt ~ -1 + SSMcustom(Z = At, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), H = Ht)
kfas.model$tol <- 0

# Case for lag 1 (needed for EM)
#kfas.model <- SSModel(yt ~ -1 + SSMcustom(Z = stack.At, T = stack.Tt, R = stack.Rt, Q = stack.Qt, a1 = stack.a1, P1 = stack.P1, P1inf = stack.P1inf), H = Ht)

# For n>1
diag.R <- unname(par.1$R)[1 + 0:(n - 1) * (n + 1)]

# Actual call to the filter/smoother?
if (any(diag.R == 0)) { # because KFAS 1.0.4 added a warning message
  ks.out <- suppressWarnings(KFS(kfas.model, simplify = FALSE))
} else {
  ks.out <- KFS(kfas.model, simplify = FALSE)
}

VtT <- ks.out$V[1:m, 1:m, , drop = FALSE]
Vtt1 <- ks.out$P[1:m, 1:m, 1:TT, drop = FALSE]
Vtt <- ks.out$Ptt[1:m, 1:m, 1:TT, drop = FALSE]






































































