library(GIGrvg)
library(breakfast)

tol <- 1e-6

nsim <- 1000
burnin <- 200

lambda2 <- 100
K <- 5
p <- 5

n <- 365

a1 <- 1
a2 <- a1
b1 <- a2
b2 <- b1

Y <- arima.sim(list(ar = 0.5), n = n)
Y[21:30] <- Y[21:30] + 3
#Y[1:50] <- Y[1:50] + 1
#Y[51:100] <- Y[51:100] + 5
#Y[101:150] <- Y[101:150] + 10
#Y[151:200] <- Y[151:200] + 20
#Y[201:250] <- Y[201:250] + 30


TT <- getT(Y, p)

Uvec <- idetect_rcpp(Y, K - 1)
U = getU(Uvec, n, K)
gamma <- colMeans(U)

init0 <- initializeGaussianPosterior(Y, lambda2, T = TT);
betaldelta0 <- init0$betadelta
tau2all0 <- init0$tau2all
sigma20 <- init0$sigma2
sigma20resi <- sigma20
lambda20 <- lambda2


init1 <- initializeGaussianPosterior(Y, lambda2, T = TT, U = U);
betadelta1 <- init1$betadelta
tau2all1 <- init1$tau2all
sigma21 <- init1$sigma2
sigma21resi <- sigma21
lambda21 <- lambda2



betadelta0out <- matrix(NA, nrow = nsim, ncol = 6)
betadelta1out <- matrix(NA, nrow = nsim, ncol = 11)
tau2all0out <- matrix(NA, nrow = nsim, ncol = 6)
tau2all1out <- matrix(NA, nrow = nsim, ncol = 11)
sigma20out <- rep(NA, nsim)
sigma21out <- rep(NA, nsim)
sigma20resiout <- rep(NA, nsim)
sigma21resiout <- rep(NA, nsim)
lambda20out <- rep(NA, nsim)
lambda21out <- rep(NA, nsim)
lambda20simout <- rep(NA, nsim)
lambda21simout <- rep(NA, nsim)

Uout <- array(NA, dim = c(n, K, nsim))
probout <- array(NA, dim = c(n, K, nsim))

rss0out <- matrix(NA, nrow = nsim, ncol = n)
rss1out <- matrix(NA, nrow = nsim, ncol = n)

rssdif <- rep(NA, nsim)

loglik0out <- matrix(NA, nrow = nsim, ncol = n)
loglik1out <- matrix(NA, nrow = nsim, ncol = n)

loglikout <- matrix(NA, nrow = nsim, ncol = n)

loglikratio <- rep(NA, nsim)

fit0out <- matrix(NA, nrow = nsim, ncol = n)
fit1out <- matrix(NA, nrow = nsim, ncol = n)


cnt <- 0
for (i in 1:(nsim + burnin)) {
  
  #cat(i, "\n")
  
  post0 <- getGaussianPosterior(Y, betaldelta0[1], tau2all0[1],
                                sigma20resi, lambda20, T = TT, beta2=betaldelta0[2:6],
                                tau2beta2 = tau2all0[2:6])
  betadelta0 <- post0$betadelta
  tau2all0 <- post0$tau2all
  sigma20 <- getSigma2(Y - post0$fit1, post0$betadelta, post0$tau2all, a1, a2)
  sigma20resi <- mean((Y - post0$fit1) ^ 2)
  lambda20 <- getLambda2EM(post0$expectedtau2all)
  #lambda20 <- getLambda2(post0$tau2all, b1, b2)
  lambda20out[cnt] <- lambda20
  lambda20simout[cnt] <- getLambda2(post0$tau2all, b1, b2)
    
  post1 <- getGaussianPosterior(Y, betadelta1[1], tau2all1[1],
                                sigma21resi, lambda21, T = TT, U = U,
                                beta2=betadelta1[2:6], delta0 = betadelta1[7:11],
                                tau2beta2=tau2all1[2:6], tau2delta0 = tau2all1[7:11])
  betadelta1 <- post1$betadelta
  tau2all1 <- post1$tau2all
  sigma21 <- getSigma2(Y - post1$fit1, post1$betadelta, post1$tau2all, a1, a2)
  sigma21resi <- mean((Y - post1$fit1) ^ 2)
  lambda21 <- getLambda2EM(post1$expectedtau2all)
  #lambda21 <- getLambda2(post1$tau2all, b1, b2)
  
  
  
  tmpUList <- getUWithoutW(Y - post1$fit0, K, post1$delta0, sigma21, gamma)
  U <- tmpUList$U
  
  
  
  if (i > burnin) {
    cnt = cnt + 1
    
    lambda21out[cnt] <- lambda21
    lambda21simout[cnt] <- getLambda2(post1$tau2all, b1, b2)
    
    rss0out[cnt, ] <- (Y - post1$fit0) ^ 2
    rss1out[cnt, ] <- (Y - post1$fit1) ^ 2
    
    rssdif[cnt] <- (sum(rss0out[cnt, ]) - sum(rss1out[cnt, ])) / sigma21
    
    loglik0out[cnt, ] <- loglikelihood(Y - post1$fit0, sigma21)
    loglik1out[cnt, ] <- loglikelihood(Y - post1$fit1, sigma21)
    
    loglikout[cnt, ] <- loglik1out[cnt, ] - loglik0out[cnt, ]
    
    loglikratio[cnt] <- sum(loglikout[cnt, ])
    
    betadelta0out[cnt, ] <- betadelta0
    tau2all0out[cnt, ] <- tau2all0
    sigma20out[cnt] <- sigma20
    sigma20resiout[cnt] <- sigma20resi
    lambda20out[cnt] <- lambda20
    
    betadelta1out[cnt, ] <- betadelta1
    tau2all1out[cnt, ] <- tau2all1
    sigma21out[cnt] <- sigma21
    sigma21resiout[cnt] <- sigma21resi
    lambda21out[cnt] <- lambda21
    
    Uout[, , cnt] <- U
    probout[, , cnt] <- tmpUList$prob
    
    fit0out[cnt, ] <- post1$fit0
    fit1out[cnt, ] <- post1$fit1
  }
  

}
