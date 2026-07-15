# Fitting priors
#library('MASS',exclude=c('select'))
library('invgamma')
library('truncnorm')
#library('univariateML')
library('fitdistrplus')
#library('survival')
prior_path = 'data/mle_coef_fits_lapse_not_fitted.csv'
prior_dist = read.csv(prior_path)
#prior_dist <- prior_dist[prior_dist$id != 1,]

# Normal fits for A (Z in MARSS notation), c (a in MARSS notation), and d (u in MARSS notation)
m<-2
n<-10

# Normal fits are just the closed form MLE expressions, using package for simplicity
mu_A <- matrix(NA,n,m)
sig_A <- matrix(NA,n,m)
mu_B <- matrix(0,m,m)
sig_B <- matrix(0,m,m)
mu_c <- matrix(NA,n,1)
sig_c <- matrix(NA,n,1)
mu_d <- matrix(NA,m,1)
sig_d <- matrix(NA,m,1)
al_r <- matrix(NA,n,1)
bet_r <- matrix(NA,n,1)
r_threshold <- 0.0001

# for(i in 1:n){
#   for(j in 1:m){
#     # Normal fitting for A values (Z in MARSS notation)
#     iter_label <- paste('z',i,j,sep='')
#     iter_dist <- prior_dist[[iter_label]]
#     est <- fitdist(iter_dist,'norm')
#     mu_A[i,j] <- est$estimate[1]
#     sig_A[i,j] <- 1/est$estimate[2]^2
#     if(i==1){
#       # Normal fitting for d values (u in MARSS notation)
#       iter_label <- paste('u',j,sep='')
#       iter_dist <- prior_dist[[iter_label]]
#       est <- fitdist(iter_dist,'norm')
#       mu_d[j] <- est$estimate[1]
#       sig_d[j] <- 1/est$estimate[2]^2
# 
#       # Truncated normal fitting for B (same in MARSS notation)
#       iter_label <- paste('b',j,j,sep='')
#       iter_dist <- prior_dist[[iter_label]]
#       lb <- -.5
#       ub <- 1
#       iter_dist <- iter_dist[iter_dist<ub & iter_dist>lb]
#       est <- fitdist(iter_dist,dtruncnorm,method='mle',start=list(mean=max(iter_dist),sd=0.5),fix.arg=list(a=lb,b=ub))
#       mu_B[j,j]=est$estimate[1]
#       sig_B[j,j] <- 1/est$estimate[2]^2
#     }
#   }
#   # Normal fitting for c values (a in MARSS notation)
#   iter_label <- paste('a',i,sep='')
#   iter_dist <- prior_dist[[iter_label]]
#   est <- fitdist(iter_dist,'norm')
#   mu_c[i] <- est$estimate[1]
#   sig_c[i] <- 1/est$estimate[2]^2
# 
#   # # Inverse gamma fitting for R (same in MARSS notation)
#   iter_label <- paste('r',i,sep='')
#   iter_dist <- prior_dist[[iter_label]]
#   lb <- 0.001
#   iter_dist <- iter_dist[iter_dist > lb]
#   est <- fitdist(r1,dinvgamma,method='mle',start=list(shape=2,rate=.1),lower=0)
#   al_r[i] <- est$estimate[1]
#   bet_r[i] <- 1/est$estimate[2]
# }

# print(mu_A)
# print(sig_A)
# print(mu_c)
# print(sig_c)
# print(mu_d)
# print(sig_d)
# if(any(is.na(mu_A))|any(is.na(sig_A))){
#   print("Fitting error, NAs remain.")
# }

# A test
var <- 'z51'
iter_dist <- prior_dist[[var]]
print(paste(min(iter_dist),max(iter_dist)))
est <- fitdist(iter_dist,'norm')
x<-seq(from=-3,to=3,by=0.01)
plot(x,dtruncnorm(x,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(iter_dist,breaks=100,freq=FALSE,add=TRUE)

# d test
var <- 'u2'
iter_dist <- prior_dist[[var]]
ub <- 10
lb <- -10
#lb <- distmean - 3*distsd
print(paste(min(iter_dist),max(iter_dist)))
print(iter_dist)
iter_dist <- iter_dist[iter_dist <ub]
iter_dist <- iter_dist[iter_dist >lb]
est <- fitdist(iter_dist,'norm')
x<-seq(from=-20,to=20,by=0.01)
plot(x,dtruncnorm(x,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,2))
iter_dist <- prior_dist[[var]]
hist(iter_dist,breaks=100,freq=FALSE,add=TRUE)
hist(iter_dist,breaks=100,freq=FALSE)

# c test
var <- 'a10'
iter_dist <- prior_dist[[var]]
distmean <- mean(iter_dist)
distsd <- sqrt(var(iter_dist))
#ub <- distmean + 3*distsd
ub <- 1
lb <- -1
#lb <- distmean - 3*distsd
print(paste(min(iter_dist),max(iter_dist)))
print(paste(lb,distmean,ub))
#iter_dist <- iter_dist[iter_dist <ub]
#iter_dist <- iter_dist[iter_dist >lb]
est <- fitdist(iter_dist,'norm')
x<-seq(from=-5,to=5,by=0.01)
plot(x,dtruncnorm(x,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(iter_dist,breaks=100,freq=FALSE,add=TRUE)

# Inverse gamma for variances
# Example case
r1 <- prior_dist$r9
x<-seq(from=0,to=0.08,by=0.0001)
lb <-r_threshold
#r1 <- r1[r1>lb]
r1[r1<lb]<-r_threshold
est <- fitdist(r1,dinvgamma,method='mle',start=list(shape=2,rate=.1),lower=0)
print(est)
plot(x,dinvgamma(x,shape=est$estimate[1],rate=est$estimate[2]),col='blue')
hist(r1,breaks=500,freq=FALSE,add=TRUE)

shape<-est$estimate[1]
rate<- est$estimate[2]
y <- 0.01
print(scale ** shape / gamma(shape) * y ** (-shape - 1) * exp(-rate / y))
print(dinvgamma(y,shape=shape,rate=rate))



# Truncated normal
# Mimic CHTC implementation
eps <- 1*10^(-10)
iter_dist <- prior_dist[["b22"]]
lb <--1+eps
ub <- 1-eps
iter_dist[iter_dist<lb]<-lb
iter_dist[iter_dist>ub]<-ub
#iter_dist <- iter_dist[iter_dist<ub & iter_dist>lb]
# hist(iter_dist,breaks=100,freq=FALSE)
est <- fitdist(iter_dist,dtruncnorm,method='mle',start=list(mean=max(iter_dist),sd=0.5),fix.arg=list(a=lb-eps,b=ub+eps))
print(est)
x<-seq(from=-1,to=1.25,by=0.01)
plot(x,dtruncnorm(x,a=-1,b=1,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(iter_dist,breaks=100,freq=FALSE,add=TRUE)

# fitdistrplus testing
dist <- prior_dist[['b22']]
print(dist)
b_ub<- 1
b_lb<- -1
alim <- b_lb
blim <- b_ub
dist[dist<b_lb]<-b_lb
dist[dist>b_ub]<-b_ub
#dist<-dist[dist<=b_ub & dist>=b_lb]
est <- fitdist(dist,dtruncnorm,method='mle',start=list(mean=0.75,sd=0.5),fix.arg=list(a=alim,b=blim))
print(est)
x<-seq(from=-1,to=1.25,by=0.01)
plot(x,dtruncnorm(x,a=alim,b=blim,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(dist,breaks=100,freq=FALSE,add=TRUE)

dist <- prior_dist[['b22']]
print(dist)
b_ub<- 1
b_lb<- -1
alim <- b_lb
blim <- b_ub
dist<-dist[dist<=b_ub & dist>=b_lb]
est <- fitdist(dist,dtruncnorm,method='mle',start=list(mean=1,sd=0.5),fix.arg=list(a=alim,b=blim))
print(est)
x<-seq(from=-1,to=1.25,by=0.01)
plot(x,dtruncnorm(x,a=alim,b=blim,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(dist,breaks=100,freq=FALSE,add=TRUE)

dist <- prior_dist[['b22']]
print(dist)
b_ub<- 1
b_lb<- -1
alim <- b_lb
blim <- b_ub
dist<-dist[dist<=b_ub & dist>=b_lb]
est <- fitdist(dist,dtruncnorm,method='mle',start=list(mean=1,sd=0.5),fix.arg=list(a=alim,b=blim))
print(est)
x<-seq(from=-1,to=1.25,by=0.01)
plot(x,dtruncnorm(x,a=alim,b=blim,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(dist,breaks=100,freq=FALSE,add=TRUE)

dist <- prior_dist[['b22']]
print(dist)
b_ub<- 1
b_lb<- 0
alim <- b_lb
blim <- b_ub
dist<-dist[dist<=b_ub & dist>=b_lb]
est <- fitdist(dist,'norm')
print(est)
x<-seq(from=-1,to=1.25,by=0.01)
plot(x,dnorm(x,est$estimate[1],est$estimate[2]),col='blue',ylim=c(0,10))
hist(dist,breaks=100,freq=FALSE,add=TRUE)

dist <- prior_dist[['b11']]
print(dist)
b_ub<- 1
b_lb<- -0.5
alim <- b_lb
blim <- b_ub
dist<-dist[dist<b_ub & dist>b_lb]
est <- fitdist(dist,dtruncnorm,method='mle',start=list(mean=max(dist),sd=0.5),fix.arg=list(a=alim,b=blim))
print(est)
x<-seq(from=-1,to=1.25,by=0.01)
plot(x,dtruncnorm(x,a=alim,b=blim,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(dist,breaks=100,freq=FALSE,add=TRUE)

dist <- prior_dist[['b22']]
print(dist)
b_ub<- 1
b_lb<- -0.5
alim <- b_lb
blim <- b_ub
dist<-dist[dist<b_ub & dist>b_lb]
est <- fitdist(dist,dtruncnorm,method='mle',start=list(mean=max(dist),sd=0.5),fix.arg=list(a=alim,b=blim))
print(est)
x<-seq(from=-1,to=1.25,by=0.01)
plot(x,dtruncnorm(x,a=alim,b=blim,mean=est$estimate[1],sd=est$estimate[2]),col='blue',ylim=c(0,10))
hist(dist,breaks=100,freq=FALSE,add=TRUE)
