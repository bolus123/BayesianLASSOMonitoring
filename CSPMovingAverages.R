library(GIGrvg)
library(breakfast)
library(pscl)


###############################

seed <- 1234#

nsim <- 1000#
burnin <- 200#

lambda2 <- 10#

p <- 10#
K <- 10#

dat <- read.csv(file = '/Users/yuihuiyao/Library/CloudStorage/Box-Box/Yuhui R21/Daily opioid-overdose-related ER visits in Walker.csv')#

dat <- dat[1:1461, ]#dat[which(dat[, 1] == "5/1/16"):which(dat[, 1] == "1/1/17"), ]

n <- dim(dat)[1]

#Y <- dat[, 2]

W <- 14

dat1 <- rep(NA, n - W - 1)

k <- 0
for (i in W:n) {
  k <- k + 1
  dat1[k] <- sum(dat[i - (0:(W - 1)), 2])
}

Yorig <- dat1

#Y <- Yorig + 0.5
Y <- Yorig / W

n <- length(Y)

###############################

betadeltamat0 <- matrix(NA, nrow = nsim, ncol = 1 + p)
taumat0 <- matrix(NA, nrow = nsim, ncol = 1 + p)
sigmamat0 <- matrix(NA, nrow = nsim, ncol = 1)
fit0mat0 <- matrix(NA, nrow = nsim, ncol = n)
fit1mat0 <- matrix(NA, nrow = nsim, ncol = n)


betadeltamat1 <- matrix(NA, nrow = nsim, ncol = 1 + p + K)
taumat1 <- matrix(NA, nrow = nsim, ncol = 1 + p + K)
sigmamat1 <- matrix(NA, nrow = nsim, ncol = 1)
fit0mat1 <- matrix(NA, nrow = nsim, ncol = n)
fit1mat1 <- matrix(NA, nrow = nsim, ncol = n)
Vmat1 <- matrix(NA, nrow = nsim, ncol = n)
Uarray1 <- array(NA, c(n, K, nsim))
probarray1 <- array(NA, c(n, K, nsim))

omegamat <- matrix(NA, nrow = nsim, ncol = 1)

lambda2mat0 <- matrix(NA, nrow = nsim, ncol = 1)
lambda2mat1 <- matrix(NA, nrow = nsim, ncol = 1)

fit0mat <- matrix(NA, nrow = nsim, ncol = n)
fit1mat <- matrix(NA, nrow = nsim, ncol = n)

sigma2mat <- matrix(NA, nrow = nsim, ncol = 1)

###############################

set.seed(seed)

###############################

V1 <- Y;
T1 <- getT(V1, p)

X1 <- cbind(1, T1)
m0 <- optim_rcpp(c(mean(V1), acf(V1, plot = FALSE)$acf[2], rep(0, p - 1 )), V1, X1)

beta00 <- m1[1]
beta02 <- m1[2:(p + 1)]

fit1 <- X1 %*% m0

sigma02 <- var(V1 - fit1)

tau0 <- getTau2(m0, sigma02, lambda2)

taubeta00 <- tau0[1]
taubeta02 <- tau0[2:(p + 1)]

###############################

V1 <- Y;
T1 <- getT(V1, p)

U1 <- idetect_rcpp(V1, K)
U1 <- getU(U1, n, K)

gamma1 <- colMeans(U1)

X1 <- cbind(1, T1, U1)
m1 <- optim_rcpp(c(mean(V1), acf(V1, plot = FALSE)$acf[2], rep(0, p - 1 + K)), V1, X1)

beta10 <- m1[1]
beta12 <- m1[2:(p + 1)]
delta10 <- m1[(p + 2):(1 + p + K)]

fit1 <- X1 %*% m1

sigma12 <- var(V1 - fit1)

tau1 <- getTau2(m1, sigma12, lambda2)

taubeta10 <- tau1[1]
taubeta12 <- tau1[2:(p + 1)]
taudelta10 <- tau1[(p + 2):(1 + p + K)]

###############################
kk <- 0
for (i in 1:(nsim + burnin)) {
  
  m0 <- getGaussianPosterior(
    V1, beta00, taubeta00, sigma02, lambda20, matrix(0, 1, 1), T1, matrix(0, 1, 1), matrix(0, 1, 1),
    matrix(0, 1, 1), beta02, matrix(0, 1, 1), matrix(0, 1, 1), matrix(0, 1, 1), taubeta02, matrix(0, 1, 1),
    matrix(0, 1, 1), 0, p, 0, 0)
  
  beta00 <- m0$betadelta[1]
  beta02 <- m0$betadelta[2:(p + 1)]
  
  sigma02 <- var(V1 - m0$fit1)
  
  taubeta00 <- m0$tau2all[1]
  taubeta02 <- m0$tau2all[2:(p + 1)]
  
  lambda20 <- getLambda2EM(c(m0$expectedtau2all))
  
  ###############################
  if (i > burnin) {
    kk <- kk + 1
    
    betadeltamat0[kk, ] <- m0$betadelta
    taumat0[kk, ] <- m0$tau2all
    sigmamat0[kk, ] <- sigma12
    fit0mat0[kk, ] <- m0$fit0
    fit1mat0[kk, ] <- m0$fit1
    
    lambda2mat0[kk, ] <- lambda20
    
    #fit0mat[kk, ] <- (1 - omega) * exp(m1$fit0)
    #fit1mat[kk, ] <- (1 - omega) * exp(m1$fit1)
    
    #sigma2mat[kk, ] <- var(Y - fit1mat[kk, ])
  }
}

###############################
kk <- 0
for (i in 1:(nsim + burnin)) {
 
  ###############################
  #for (j in 1:(burnin + 1)) {
  m1 <- getGaussianPosterior(
    V1, beta10, taubeta10, sigma12, lambda21, matrix(0, 1, 1), T1, U1, matrix(0, 1, 1),
    matrix(0, 1, 1), beta12, delta10, matrix(0, 1, 1), matrix(0, 1, 1), taubeta12, taudelta10,
    matrix(0, 1, 1), 0, p, K, 0)
  
  beta10 <- m1$betadelta[1]
  beta12 <- m1$betadelta[2:(p + 1)]
  delta10 <- m1$betadelta[(p + 2):(1 + p + K)]
  
  sigma12 <- var(V1 - m1$fit1)
  
  taubeta10 <- m1$tau2all[1]
  taubeta12 <- m1$tau2all[2:(p + 1)]
  taudelta10 <- m1$tau2all[(p + 2):(1 + p + K)]
  
  fit1 <- m1$fit1
  #}
  zetadelta1 <- V1 - m1$fit0
  tmpU1 <- getUWithoutW(zetadelta1, K, delta10, sigma12, gamma1)
  U1 <- tmpU1$U
  prob1 <- tmpU1$prob
  ###############################
  ###############################
  lambda21 <- getLambda2EM(c(m1$expectedtau2all))
  ###############################
  if (i > burnin) {
    kk <- kk + 1
    
    
    betadeltamat1[kk, ] <- m1$betadelta
    taumat1[kk, ] <- m1$tau2all
    sigmamat1[kk, ] <- sigma12
    fit0mat1[kk, ] <- m1$fit0
    fit1mat1[kk, ] <- m1$fit1
    Uarray1[,,kk] <- U1
    probarray1[,,kk] <- prob1
    
    lambda2mat1[kk, ] <- lambda21
    
    #fit0mat[kk, ] <- (1 - omega) * exp(m1$fit0)
    #fit1mat[kk, ] <- (1 - omega) * exp(m1$fit1)
    
    #sigma2mat[kk, ] <- var(Y - fit1mat[kk, ])
  }
}

###############################

tol <- 1e-6
perc <- 0.8

newdelta <- t(betadeltamat1[, (1 + p + 1):(1 + p + K)])

perccl <- rep(NA, K - 2)

for (i in 2:(K - 1)) {
  cl <- kmeans(newdelta, centers = i)
  perccl[i - 1] <- cl$betweenss / cl$totss
}

k <- which(perccl >= perc)[1]

cl <- kmeans(newdelta, centers = k + 1)

newU1 <- matrix(0, nrow = n, ncol = K)

for (i in 1:(k + 1)) {
  indx <- which(cl$cluster == i)
  nindx <- length(indx)
  newU1[which(U1[, indx[1]] == 1), indx[1]] <- 1
  if (nindx > 1) {
    for (j in 2:nindx){
      newU1[which(U1[, indx[j]] == 1), indx[1]] <- 1
      delta10[indx[j]] <- 0
    }
  }
}

U1 <- newU1
gamma1 <- colMeans(U1)
gamma1[which(gamma1 == 0)] <- 1e-6
gamma1 <- gamma1 / (sum(gamma1))

###############################

kk <- 0
for (i in 1:(nsim + burnin)) {
  ###############################
  #for (j in 1:(burnin + 1)) {
  
  
  m1 <- getGaussianPosterior(
    V1, beta10, taubeta10, sigma12, lambda2, matrix(0, 1, 1), T1, U1, matrix(0, 1, 1),
    matrix(0, 1, 1), beta12, delta10, matrix(0, 1, 1), matrix(0, 1, 1), taubeta12, taudelta10,
    matrix(0, 1, 1), 0, p, K, 0)
  
  beta10 <- m1$betadelta[1]
  beta12 <- m1$betadelta[2:(p + 1)]
  delta10 <- m1$betadelta[(p + 2):(1 + p + K)]
  
  sigma12 <- var(V1 - m1$fit1)
  
  taubeta10 <- m1$tau2all[1]
  taubeta12 <- m1$tau2all[2:(p + 1)]
  taudelta10 <- m1$tau2all[(p + 2):(1 + p + K)]
  
  fit1 <- m1$fit1
  #}
  zetadelta1 <- V1 - m1$fit0
  tmpU1 <- getUWithoutW(zetadelta1, K, delta10, sigma12, gamma1)
  U1 <- tmpU1$U
  prob1 <- tmpU1$prob
  ###############################
  ###############################
  lambda2 <- getLambda2EM(c(m1$expectedtau2all))
  ###############################
  if (i > burnin) {
    kk <- kk + 1
    
    betadeltamat1[kk, ] <- m1$betadelta
    taumat1[kk, ] <- m1$tau2all
    sigmamat1[kk, ] <- sigma12
    fit0mat1[kk, ] <- m1$fit0
    fit1mat1[kk, ] <- m1$fit1
    Uarray1[,,kk] <- U1
    probarray1[,,kk] <- prob1
    
    lambda2mat1[kk, ] <- lambda2
    
    #fit0mat[kk, ] <- (1 - omega) * exp(m1$fit0)
    #fit1mat[kk, ] <- (1 - omega) * exp(m1$fit1)
    
    #sigma2mat[kk, ] <- var(Y - fit1mat[kk, ])
  }
}

###############################
#
#nsim1 <- 10000
#
#aa <- rep(NA, nsim1) 
#bb1 <- matrix(NA, nsim1, n)
#bb0 <- matrix(NA, nsim1, n)
#
#for (i in 1:nsim1) {
#  samp <- sample(1:nsim, 1)
#  beta10 <- betadeltamat1[samp, 1]
#  taubeta10 <- taumat1[samp, 1]
#  sigma12 <- sigmamat1[samp, ]
#  lambda2 <- lambda2mat[samp, ]
#  U1 <- Uarray1[, , samp]
#  beta12 <- betadeltamat1[samp, 2:(p + 1)]
#  delta10 <- betadeltamat1[samp, (p + 2):(p + 1 + K)]
#  taubeta12 <- taumat1[samp, 2:(p + 1)]
#  taudelta10 <- taumat1[samp, (p + 2):(p + 1 + K)]
#  bb1[i, ] <- loglikelihoodLayer3(Y, beta10, taubeta10,
#                                 sigma12, lambda2, 
#                                 gamma1,
#                                 X=NULL,
#                                 T=T1,
#                                 U=U1,
#                                 W=NULL,
#                                 beta1=NULL, 
#                                 beta2=beta12, 
#                                 delta0=delta10, 
#                                 delta1=NULL,
#                                 tau2beta1=NULL,
#                                 tau2beta2=taubeta12,
#                                 tau2delta0=taudelta10,
#                                 tau2delta1=NULL)
#  
#  bb0[i, ] <- loglikelihoodLayer3(Y, beta10, taubeta10,
#                                  sigma12, lambda2, 
#                                  gamma1,
#                                  X=NULL,
#                                  T=T1,
#                                  U=NULL,
#                                  W=NULL,
#                                  beta1=NULL, 
#                                  beta2=beta12, 
#                                  delta0=NULL, 
#                                  delta1=NULL,
#                                  tau2beta1=NULL,
#                                  tau2beta2=taubeta12,
#                                  tau2delta0=NULL,
#                                  tau2delta1=NULL)
#  aa[i] <- -2 * (sum(bb0[i, ]) - sum(bb1[i, ]))
#}
#
#1 - mean(aa > 0)

###############################

TT <- matrix(NA, nrow = nsim, ncol = 1)

for (i in 1:nsim) {
  TT[i, ] <- sum((Y - fit0mat0[i, ]) ^ 2 / sigmamat0[i, ]) - 
    sum((Y - fit1mat1[i, ]) ^ 2 / sigmamat1[i, ])
    
}

(mean(TT > 0))

###############################


###############################

2 * min(mean(betadeltamat1[, 12] > 0), mean(betadeltamat1[, 12]  < 0))
2 * min(mean(betadeltamat1[, 13] > 0), mean(betadeltamat1[, 13]  < 0))
2 * min(mean(betadeltamat1[, 14] > 0), mean(betadeltamat1[, 14]  < 0))
2 * min(mean(betadeltamat1[, 15] > 0), mean(betadeltamat1[, 15]  < 0))
2 * min(mean(betadeltamat1[, 16] > 0), mean(betadeltamat1[, 16]  < 0))
2 * min(mean(betadeltamat1[, 17] > 0), mean(betadeltamat1[, 17]  < 0))
2 * min(mean(betadeltamat1[, 18] > 0), mean(betadeltamat1[, 18]  < 0))
2 * min(mean(betadeltamat1[, 19] > 0), mean(betadeltamat1[, 19]  < 0))
2 * min(mean(betadeltamat1[, 20] > 0), mean(betadeltamat1[, 20]  < 0))
2 * min(mean(betadeltamat1[, 21] > 0), mean(betadeltamat1[, 21]  < 0))

###############################
probll <- matrix(NA, nrow = n, ncol = K)
U1vec <- rep(NA, n)
for (i in 1:n) {
  for (j in 1:K) {
    probll[i, j] <- median(probarray1[i, j, ], na.rm = T) 
  }
  U1vec[i] <- which.max(probll[i, ])
}

U1vec[!(U1vec %in% c(1, 2, 3, 4, 5, 9))] <- 0

###############################
#plot(Y[500:1000], type = 'l')
#
#col2 <- rep('black', n)
#col2[which(U1vec == 2)] <- NA
#col2[which(U1vec == 4)] <- NA
#points(Y[500:1000], col = col2[500:1000], cex = 0.5)
#  
#col1 <- rep(NA, n)
#col1[which(U1vec == 2)] <- "red"
#col1[which(U1vec == 4)] <- "blue"
#points(Y[500:1000], col = col1[500:1000], cex = 1.5)
#      
###############################
#YY <- dat[1:1461, 2]
#YY <- YY[514:1014]
#
#plot(YY, type = 'l')
#
#col2 <- rep('black', n)
#col2[which(U1vec == 2)] <- NA
#col2[which(U1vec == 4)] <- NA
#points(YY, col = col2[500:1000], cex = 0.5)
#
#col1 <- rep(NA, n)
#col1[which(U1vec == 2)] <- "red"
#  col1[which(U1vec == 4)] <- "blue"
#    points(YY, col = col1[500:1000], cex = 1.5)
#    
###############################
YY <- dat[732:1096, 2]
UU <- c(rep(0, W - 1), U1vec)
UU <- U1vec[732:1096]

UU[UU == 1] <- 0

plot(YY, type = 'l', main = "Walker County, AL 2018", 
     ylab = 'Opioid-overdose-related ER Visits', xaxt = "n", xlab = "")

axis(1, at = c(
  1, 
  32, 
  60,
  91,
  121,
  152,
  182,
  213,
  244,
  274,
  305,
  335
), labels = c(
  "Jan 1",
  "Feb 1",
  "Mar 1",
  "Apr 1",
  "May 1",
  "Jun 1",
  "Jul 1",
  "Aug 1",
  "Sep 1",
  "Oct 1",
  "Nov 1",
  "Dec 1"
))


col2 <- rep('black', n)
col2[which(UU == 2)] <- NA
col2[which(UU == 3)] <- NA
col2[which(UU == 5)] <- NA
col2[which(UU == 9)] <- NA
points(YY, col = col2, cex = 0.7, pch = 16)

col1 <- rep(NA, n)
#col1[which(UU == 1)] <- "green"

col1[which(UU == 3)] <- "red"
  col1[which(UU == 4)] <- "red"
    col1[which(UU == 5)] <- "red"
      
    #col1[which(UU == 2)] <- "green"
    #      col1[which(UU == 9)] <- "green"
    
    points(YY, col = col1, cex = 1.2, pch = 16)
    
    legend("topright", 
           legend = c("High risk", "Mitigated or low risk"), 
           col = c('red', "black"), pch = c(16, 16, 16))
###############################
YY <- c(rep(0, W - 1), dat1)
YY <- YY[732:1096]
UU <- c(rep(0, W - 1), U1vec)
UU <- U1vec[732:1096]

plot(YY, type = 'l', main = "Walker County, AL 2018", 
     ylab = 'Two-week Moving Average', xaxt = "n", xlab = "")

axis(1, at = c(
  1, 
  32, 
  60,
  91,
  121,
  152,
  182,
  213,
  244,
  274,
  305,
  335
), labels = c(
  "Jan 1",
  "Feb 1",
  "Mar 1",
  "Apr 1",
  "May 1",
  "Jun 1",
  "Jul 1",
  "Aug 1",
  "Sep 1",
  "Oct 1",
  "Nov 1",
  "Dec 1"
))

col2 <- rep('black', n)
col2[which(UU == 2)] <- NA
col2[which(UU == 3)] <- NA
col2[which(UU == 5)] <- NA
col2[which(UU == 9)] <- NA
points(YY, col = col2, cex = 0.7, pch = 16)

col1 <- rep(NA, n)
#col1[which(UU == 1)] <- "green"

col1[which(UU == 3)] <- "red"
  col1[which(UU == 4)] <- "red"
  col1[which(UU == 5)] <- "red"
    
#col1[which(UU == 2)] <- "green"
#      col1[which(UU == 9)] <- "green"
        
      points(YY, col = col1, cex = 1.2, pch = 16)
    
    legend("topright", 
           legend = c("High risk", "Mitigated or low risk"), 
           col = c('red', "black"), pch = c(16, 16, 16))
    
###############################
