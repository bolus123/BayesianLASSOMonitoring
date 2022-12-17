library(GIGrvg)
library(breakfast)
library(pscl)
#
tol <- 1e-6
#
nsim <- 100
burnin <- 100
#
lambda2 <- 100
K <- 20
p <- 5

spl <- 0.9

#
#n <- 365
#
a1 <- 1
a2 <- a1
b1 <- a2
b2 <- b1
#
d1 <- 1
d2 <- d1
#
#V1 <- arima.sim(list(ar = 0.5), n = n)
#V1[201:230] <-  V1[201:230] + 2
##Y[1:50] <- Y[1:50] + 1
##Y[51:100] <- Y[51:100] + 5
##Y[101:150] <- Y[101:150] + 10
##Y[151:200] <- Y[151:200] + 20
##Y[201:250] <- Y[201:250] + 30
#
#V1orig <- V1
#
#psi <- 1
#
#mu <- exp(V1)
#
#Y1 <- rnbinom(n, psi, mu = mu)
#
#V2 <- arima.sim(list(ar = 0.5), n = n)
#V2[101:130] <-  V2[101:130] + 20
#
#V2orig <- V2
#
#Y2 <- as.numeric(1 / (1 + exp(- V2)) < 0.5)
#
#Y <- rep(NA, n)
#
#for (i in 1:n) {
#  if (Y2[i] == 1) {
#    Y[i] <- 0
#  } else {
#    Y[i] = Y1[i]
#  }
#}


##########################

dat <- read.csv(file = '/Users/yuihuiyao/Library/CloudStorage/Box-Box/Yuhui R21/Daily opioid-overdose-related ER visits in Walker.csv')

Y <- dat[, 2]

n <- length(Y)

undebug(posteriorNB)
bb <- posteriorNB(Y, X = NULL, W = NULL, p = 5, K = 5, lambda = 10,
                  c1 = 1, c2 = 1,
                  updatelambda = TRUE, 
                  updatepsi = TRUE, 
                  clustering = "kmeans", spl = 0.9, 
                  changepoint = "idetect", 
                  burnin1 = 50, burnin2 = 50, nsim = 1000, tol = 1e-6)


undebug(posteriorZinfNB)
aa <- posteriorZinfNB(Y, X = NULL, W = NULL, p = 5, K = 5, lambda = 10,
                            c1 = 1, c2 = 1,
                            updatelambda = TRUE, 
                            updatepsi = TRUE, 
                            clustering = "kmeans", spl = 0.9, 
                            changepoint = "idetect", 
                            burnin1 = 50, burnin2 = 50, nsim = 1000, tol = 1e-6) 


##########################

V1 <- log(Y + 0.5)
TT1 <- getT(V1, p)

Uvec <- idetect_rcpp(V1, K)
U1 <- getU(Uvec, n, K)
gamma1 <- colMeans(U1)

init1 <- initializeGaussianPosterior(V1, lambda2, T = TT1, U = U1);
psi <- exp(init1$sigma2)

betadelta1 <- init1$betadelta
tau2all1 <- init1$tau2all
sigma21 <- init1$sigma2

fit1 <- init1$fit1

##########################

V2 <- as.numeric(Y == 0)
TT2 <- getT(V2, p)

Uvec <- idetect_rcpp(V2, K)
U2 <- getU(Uvec, n, K)
gamma2 <- colMeans(U2)

dat <- as.data.frame(cbind(V2, TT2, U2))
names(dat)[1] <- 'VV'

init2 <- glm(VV~., data = dat, family = binomial())

betadelta2 <- init2$coefficients
betadelta2[is.na(betadelta2)] <- tol
sigma22 <- sd(init2$residuals)
tau2all2 <- getTau2(betadelta2, sigma22, lambda2)

V2 <- rnorm(n, init2$fitted.values, sqrt(sigma22))

fit2 <- init2$fitted.values
##########################


#TT2 <- getT(V2, p)
#
#init2 <- initializeGaussianPosterior(V2, lambda2, T = TT2, U = U2);
#betadelta2 <- init2$betadelta
#tau2all2 <- init2$tau2all
#sigma22 <- init2$sigma2
#

betadelta1out <- matrix(NA, nrow = nsim, ncol = 1 + p + K)
betadelta2out <- matrix(NA, nrow = nsim, ncol = 1 + p + K)
tau2all1out <- matrix(NA, nrow = nsim, ncol = 1 + p + K)
tau2all2out <- matrix(NA, nrow = nsim, ncol = 1 + p + K)
sigma21out <- rep(NA, nsim)
sigma22out <- rep(NA, nsim)
lambda2out <- rep(NA, nsim)

U1out <- array(NA, dim = c(n, K, nsim))
prob1out <- array(NA, dim = c(n, K, nsim))

U2out <- array(NA, dim = c(n, K, nsim))
prob2out <- array(NA, dim = c(n, K, nsim))

V1out <- matrix(NA, nrow = nsim, ncol = n)
V2out <- matrix(NA, nrow = nsim, ncol = n)

fit10out <- matrix(NA, nrow = nsim, ncol = n)
fit11out <- matrix(NA, nrow = nsim, ncol = n)

fit20out <- matrix(NA, nrow = nsim, ncol = n)
fit21out <- matrix(NA, nrow = nsim, ncol = n)

fit0out <- matrix(NA, nrow = nsim, ncol = n)
fit1out <- matrix(NA, nrow = nsim, ncol = n)

V1out <- matrix(NA, nrow = nsim, ncol = n)
V2out <- matrix(NA, nrow = nsim, ncol = n)

psiout <- rep(NA, nsim)

df0out <- rep(NA, nsim)
df1out <- rep(NA, nsim)

loglik0 <- rep(NA, nsim)
loglik1 <- rep(NA, nsim)

cnt <- 0

test <- rep(NA, nsim)

for (i in 1:(nsim + burnin)) {
  
  V1 <- getV1(Y, V1, V2, psi, fit1, sigma21, burnin)
  
  TT1 <- getT(V1, p)
  
  post1 <- getGaussianPosterior(V1, betadelta1[1], tau2all1[1],
                                sigma21, lambda2, T = TT1, U = U1,
                                beta2=betadelta1[2:6], delta0 = betadelta1[7:(1 + p + K)],
                                tau2beta2=tau2all1[2:6], tau2delta0 = tau2all1[7:(1 + p + K)])
  
  
  betadelta1 <- post1$betadelta
  tau2all1 <- post1$tau2all
  sigma21 <- mean((V1 - post1$fit1) ^ 2)
  fit1 <- post1$fit1
  
  tmpUList1 <- getUWithoutW(V1 - post1$fit0, K, post1$delta0, sigma21, gamma1)
  U1 <- tmpUList1$U
  prob1 <- tmpUList1$prob
  
  
  
  V2 <- getV2(Y, V1, V2, psi, fit2, sigma22, burnin)
  
  
  TT2 <- getT(V2, p)
  
  post2 <- getGaussianPosterior(V2, betadelta2[1], tau2all2[1],
                                sigma22, lambda2, T = TT2, U = U2,
                                beta2=betadelta2[2:6], delta0 = betadelta2[7:(1 + p + K)],
                                tau2beta2=tau2all2[2:6], tau2delta0 = tau2all2[7:(1 + p + K)])
  
  betadelta2 <- post2$betadelta
  tau2all2 <- post2$tau2all
  sigma22 <- mean((V2 - post2$fit1) ^ 2)
  fit2 <- post2$fit1
  
  tmpUList2 <- getUWithoutW(V2 - post2$fit0, K, post2$delta0, sigma22, gamma2)
  U2 <- tmpUList2$U
  prob2 <- tmpUList2$prob
  
  
  psi <- getPsi(Y, V1, V2, psi, burnin, d1, d2)
  
  lambda2 <- getLambda2EM(c(post1$expectedtau2all, post2$expectedtau2all))
  
  
  if (i > burnin) {
    cnt = cnt + 1
    
    lambda2out[cnt] <- lambda2
    psiout[cnt] <- psi
    
    betadelta1out[cnt, ] <- betadelta1
    tau2all1out[cnt, ] <- tau2all1
    sigma21out[cnt] <- sigma21
    
    U1out[, , cnt] <- U1
    prob1out[, , cnt] <- prob1
    
    betadelta2out[cnt, ] <- betadelta2
    tau2all2out[cnt, ] <- tau2all2
    sigma22out[cnt] <- sigma22
    
    U2out[, , cnt] <- U2
    prob2out[, , cnt] <- prob2
    
    fit10out[cnt, ] <- post1$fit0
    fit11out[cnt, ] <- post1$fit1
    
    fit20out[cnt, ] <- post2$fit0
    fit21out[cnt, ] <- post2$fit1
    
    fit0out[cnt, ] <- (1 - 1 / (1 + exp(-post2$fit0))) * exp(post1$fit0)
    fit1out[cnt, ] <- (1 - 1 / (1 + exp(-post2$fit1))) * exp(post1$fit1)
    
    V1out[cnt, ] <- V1
    V2out[cnt, ] <- V2
    
    df0out[cnt] <- sum(abs(betadelta1[1:6]) > tol) + sum(abs(betadelta2[1:6]) > tol)
    df1out[cnt] <- sum(abs(betadelta1) > tol) + sum(abs(betadelta2) > tol)
  }
  
  
}

##########################

perc1 <- rep(NA, length(K - 2))
perc2 <- rep(NA, length(K - 2))

kk <- 0
for (i in 2:(K - 1)) {
  kk <- kk + 1
  cluster1 <- kmeans(t(betadelta1out[, 7:(1 + p + K)]), centers = i)
  cluster2 <- kmeans(t(betadelta2out[, 7:(1 + p + K)]), centers = i)
  perc1[kk] <- cluster1$betweenss / cluster1$totss
  perc2[kk] <- cluster2$betweenss / cluster2$totss
}

target1 <- which(perc1 > spl)[1]
target2 <- which(perc2 > spl)[1]

cluster1 <- kmeans(t(betadelta1out[, 7:(1 + p + K)]), centers = target1 + 1)
cluster2 <- kmeans(t(betadelta2out[, 7:(1 + p + K)]), centers = target2 + 1)

maxgroup1 <- max(cluster1$cluster)

U1new <- matrix(0, nrow = dim(U1)[1], ncol = dim(U1)[2])

for (i in 1:maxgroup1) {
  tmptarget <- which(cluster1$cluster == i)
  for (j in tmptarget) {
    U1new[which(U1[, j] == 1), tmptarget[1]] <- 1
  }
}

maxgroup2 <- max(cluster2$cluster)

U2new <- matrix(0, nrow = dim(U2)[1], ncol = dim(U2)[2])

for (i in 1:maxgroup2) {
  tmptarget <- which(cluster2$cluster == i)
  for (j in tmptarget) {
    U2new[which(U2[, j] == 1), tmptarget[1]] <- 1
  }
}

##########################

U1 <- U1new
U2 <- U2new

gamma1 <- colMeans(U1)
gamma2 <- colMeans(U2)

##########################

cnt <- 0
for (i in 1:(nsim + burnin)) {
  
  V1 <- getV1(Y, V1, V2, psi, fit1, sigma21, burnin)
  
  TT1 <- getT(V1, p)
  
  post1 <- getGaussianPosterior(V1, betadelta1[1], tau2all1[1],
                                sigma21, lambda2, T = TT1, U = U1,
                                beta2=betadelta1[2:6], delta0 = betadelta1[7:(1 + p + K)],
                                tau2beta2=tau2all1[2:6], tau2delta0 = tau2all1[7:(1 + p + K)])
  
  
  betadelta1 <- post1$betadelta
  tau2all1 <- post1$tau2all
  sigma21 <- mean((V1 - post1$fit1) ^ 2)
  fit1 <- post1$fit1
  
  tmpUList1 <- getUWithoutW(V1 - post1$fit0, K, post1$delta0, sigma21, gamma1)
  U1 <- tmpUList1$U
  prob1 <- tmpUList1$prob
  
  
  
  V2 <- getV2(Y, V1, V2, psi, fit2, sigma22, burnin)
  
  
  TT2 <- getT(V2, p)
  
  post2 <- getGaussianPosterior(V2, betadelta2[1], tau2all2[1],
                                sigma22, lambda2, T = TT2, U = U2,
                                beta2=betadelta2[2:6], delta0 = betadelta2[7:(1 + p + K)],
                                tau2beta2=tau2all2[2:6], tau2delta0 = tau2all2[7:(1 + p + K)])
  
  betadelta2 <- post2$betadelta
  tau2all2 <- post2$tau2all
  sigma22 <- mean((V2 - post2$fit1) ^ 2)
  fit2 <- post2$fit1
  
  tmpUList2 <- getUWithoutW(V2 - post2$fit0, K, post2$delta0, sigma22, gamma2)
  U2 <- tmpUList2$U
  prob2 <- tmpUList2$prob
  
  
  psi <- getPsi(Y, V1, V2, psi, burnin, d1, d2)
  
  lambda2 <- getLambda2EM(c(post1$expectedtau2all, post2$expectedtau2all))
  
  
  if (i > burnin) {
    cnt = cnt + 1
    
    lambda2out[cnt] <- lambda2
    psiout[cnt] <- psi
    
    betadelta1out[cnt, ] <- betadelta1
    tau2all1out[cnt, ] <- tau2all1
    sigma21out[cnt] <- sigma21
    
    U1out[, , cnt] <- U1
    prob1out[, , cnt] <- prob1
    
    betadelta2out[cnt, ] <- betadelta2
    tau2all2out[cnt, ] <- tau2all2
    sigma22out[cnt] <- sigma22
    
    U2out[, , cnt] <- U2
    prob2out[, , cnt] <- prob2
    
    fit10out[cnt, ] <- post1$fit0
    fit11out[cnt, ] <- post1$fit1
    
    fit20out[cnt, ] <- post2$fit0
    fit21out[cnt, ] <- post2$fit1
    
    fit0out[cnt, ] <- (1 - 1 / (1 + exp(-post2$fit0))) * exp(post1$fit0)
    fit1out[cnt, ] <- (1 - 1 / (1 + exp(-post2$fit1))) * exp(post1$fit1)
    
    V1out[cnt, ] <- V1
    V2out[cnt, ] <- V2
   
    df0out[cnt] <- sum(abs(betadelta1[1:6]) > tol) + sum(abs(betadelta2[1:6]) > tol)
    df1out[cnt] <- sum(abs(betadelta1) > tol) + sum(abs(betadelta2) > tol)
    
    V10 <- rnorm(n, mean = fit1, sd = sigma21)
    V20 <- rnorm(n, mean = fit2, sd = sigma22)
    
    loglik0[cnt] <- loglikelihoodLayer1(Y, V10, V20, psi)
    loglik1[cnt] <- loglikelihoodLayer1(Y, V1, V2, psi)
  }
  
  
}

plot(Y, type = 'l')

for (i in 1:nsim){
  points(fit1out[i, ], col = 'blue')
}

for (i in 1:nsim){
  points(fit0out[i, ], col = 'red')
}
