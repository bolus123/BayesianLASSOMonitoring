library(GIGrvg)
library(breakfast)

assignBetadelta <- function(list, q, p, K, r) {
  beta0 <- list$beta0
  if (q > 0) {
    beta1 <- list$beta1
  } else {
    beta1 <- NULL
  }
  if (p > 0) {
    beta2 <- list$beta2
  } else {
    beta2 <- NULL
  }
  if (K > 0) {
    delta0 <- list$delta0
    if (r > 0) {
      delta1 <- list$delta1
    } else {
      delta1 <- NULL
    }
  } else {
    delta0 <- NULL
    delta1 <- NULL
  }
  out <- list(
    "beta0" = beta0,
    "beta1" = beta1,
    "beta2" = beta2,
    "delta0" = delta0,
    "delta1" = delta1
  )
}

assignTau2 <- function(list, q, p, K, r) {
  tau2beta0 <- list$tau2beta0
  if (q > 0) {
    tau2beta1 <- list$tau2beta1
  } else {
    tau2beta1 <- NULL
  }
  if (p > 0) {
    tau2beta2 <- list$tau2beta2
  } else {
    tau2beta2 <- NULL
  }
  if (K > 0) {
    tau2delta0 <- list$tau2delta0
    if (r > 0) {
      tau2delta1 <- list$tau2delta1
    } else {
      tau2delta1 <- NULL
    }
  } else {
    tau2delta0 <- NULL
    tau2delta1 <- NULL
  }
  out <- list(
    "tau2beta0" = tau2beta0,
    "tau2beta1" = tau2beta1,
    "tau2beta2" = tau2beta2,
    "tau2delta0" = tau2delta0,
    "tau2delta1" = tau2delta1
  )
}

posteriorLayer2 <- function(Y, X = NULL, W = NULL, p = 5, K = 5, lambda = 10,
                            updatelambda = TRUE, 
                            clustering = "kmeans", spl = 0.9, 
                            changepoint = "idetect", 
                            burnin = 50, nsim = 1000, tol = 1e-6) {
  
  lambda2 = lambda ^ 2
  
  if (p > 0) {
    TT <- getT(Y, p)
  } else {
    TT <- NULL
  }
  
  if (K > 0) {
    if (changepoint == "idetect") {
      Uvec <- idetect_rcpp(Y, K - 1)
      U <- getU(Uvec, n, K)
      gamma <- colMeans(U)
    }
  } else {
    U <- NULL
    gamma <- NULL
  }
  
  if (is.null(X)) {
    q <- 0
  } else {
    q <- dim(X)[2]
  }
  
  if (is.null(W)) {
    r <- 0
  } else {
    r <- dim(W)[2]
  }
  
  m <- 1 + q + p + K + K * r
  
  ## get the model for the no-shift situation
  
  initNoShift <- initializeGaussianPosterior(V = Y, lambda2 = lambda2, X = X, T = TT);
  
  beta0 <- initNoShift$betadelta[1]
  tau2beta0 <- initNoShift$tau2all[1]
  sigma2 <- initNoShift$sigma2
  
  tmpbetadelta <- readbetadelta(initNoShift$betadelta, q, p, 0, 0)
  tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, 0, 0)
  
  tmptau2all <- readtau2all(initNoShift$tau2all, q, p, 0, 0)
  tmptau2all <- assignTau2(tmptau2all, q, p, 0, 0)
  
  layer2NoShift <- getPosteriorLayer2NoShift(
    V = Y, beta0 = tmpbetadelta$beta0, tau2beta0 = tmptau2all$tau2beta0, 
    sigma2 = initNoShift$sigma2, lambda2 = lambda2, updatelambda2 = updatelambda, 
    burnin = burnin, nsim = nsim,
    X = X, T = TT, beta1 = tmpbetadelta$beta1, beta2 = tmpbetadelta$beta2, 
    tau2beta1 = tmptau2all$tau2beta1, tau2beta2 = tmptau2all$tau2beta2
  )
  
  ## get the model for the shifted situation
  
  if (K > 0) {
    initShift <- initializeGaussianPosterior(V = Y, lambda2 = lambda2, X = X, 
                                             T = TT, U = U, W = W);
    
    tmpbetadelta <- readbetadelta(initShift$betadelta, q, p, K, r)
    tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, K, r)
    
    tmptau2all <- readtau2all(initShift$tau2all, q, p, K, r)
    tmptau2all <- assignTau2(tmptau2all, q, p, K, r)
    
    layer2Shift <- getPosteriorLayer2Shift(
      V = Y, beta0 = tmpbetadelta$beta0, tau2beta0 = tmptau2all$tau2beta0, 
      sigma2 = initShift$sigma2, lambda2 = lambda2, updatelambda2 = updatelambda, 
      burnin = burnin, nsim = nsim,
      gamma = gamma,
      X = X, T = TT, U = U, W = W, 
      beta1 = tmpbetadelta$beta1, beta2 = tmpbetadelta$beta2, 
      delta0 = tmpbetadelta$delta0, delta1 = tmpbetadelta$delta1,
      tau2beta1 = tmptau2all$tau2beta1, tau2beta2 = tmptau2all$tau2beta2,
      tau2delta0 = tmptau2all$tau2delta0, tau2delta1 = tmptau2all$tau2delta1
    )
    
    ## rebuild the estimated time
    
    perc <- rep(0, K + K * r - 2)
    
    start <- (1 + q + p + 1)
    delta <- layer2Shift$betadelta[, start:m]
    
    if (clustering == "kmeans") {
      #tmpU <- getUMaxProb(layer2Shift$prob)
      
      kk <- 0
      for (i in 2:(K + K * r - 1)) {
        kk <- kk + 1
        cluster <- kmeans(t(delta), centers = i)
        perc[kk] <- cluster$betweenss / cluster$totss
      }
      
      #cat(perc, '\n')
      
      target <- which(perc > spl)
      
      if (length(target) > 0) {
        cluster <- kmeans(t(delta), centers = target[1] + 1)
        maxgroup <- max(cluster$cluster)
        tmpU <- matrix(0, nrow = n, ncol = K)
        
        #cat(maxgroup, '\n')
        
        for (i in 1:maxgroup) {
          tmptarget <- which(cluster$cluster == i)
          for (j in tmptarget) {
            tmpU[which(U[, j] == 1), tmptarget[1]] <- 1
          }
        }
        
        tmpgamma <- colMeans(tmpU)
        tmpgamma[which(tmpgamma == 0)] <- tol
        tmpgamma <- tmpgamma / sum(tmpgamma)
                
        tmpbetadelta <- readbetadelta(layer2Shift$betadelta[nsim, ], q, p, K, r)
        tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, K, r)
        
        tmptau2all <- readtau2all(layer2Shift$tau2all[nsim, ], q, p, K, r)
        tmptau2all <- assignTau2(tmptau2all, q, p, K, r)
        
        layer2Shift <- getPosteriorLayer2Shift(
          V = Y, beta0 = tmpbetadelta$beta0, tau2beta0 = tmptau2all$tau2beta0, 
          sigma2 = initShift$sigma2, lambda2 = lambda2, updatelambda2 = updatelambda, 
          burnin = burnin, nsim = nsim,
          gamma = tmpgamma,
          X = X, T = TT, U = tmpU, W = W, 
          beta1 = tmpbetadelta$beta1, beta2 = tmpbetadelta$beta2, 
          delta0 = tmpbetadelta$delta0, delta1 = tmpbetadelta$delta1,
          tau2beta1 = tmptau2all$tau2beta1, tau2beta2 = tmptau2all$tau2beta2,
          tau2delta0 = tmptau2all$tau2delta0, tau2delta1 = tmptau2all$tau2delta1
        )
        
      } else {
        break
      }
      
    }
    
    
  }
  
  out <- list(
    "NoShift" = layer2NoShift,
    "Shift" = layer2Shift
  )
  return(out)
  
  
}


posteriorLayer2NB <- function(Y, X = NULL, W = NULL, p = 5, K = 5, lambda = 10,
                            c1 = 1, c2 = 1,
                            updatelambda = TRUE, 
                            updatepsi = TRUE, 
                            clustering = "kmeans", spl = 0.9, 
                            changepoint = "idetect", 
                            burnin1 = 50, burnin2 = 50, nsim = 1000, tol = 1e-6) {
  
  n <- length(Y)
  
  lambda2 = lambda ^ 2
  
  V1 <- log(Y + 0.5)
  
  if (p > 0) {
    TT <- getT(V1, p)
  } else {
    TT <- NULL
  }
  
  if (K > 0) {
    if (changepoint == "idetect") {
      Uvec <- idetect_rcpp(Y, K - 1)
      U <- getU(Uvec, n, K)
      gamma <- colMeans(U)
    }
  } else {
    U <- NULL
    gamma <- NULL
  }
  
  if (is.null(X)) {
    q <- 0
  } else {
    q <- dim(X)[2]
  }
  
  if (is.null(W)) {
    r <- 0
  } else {
    r <- dim(W)[2]
  }
  
  m <- 1 + q + p + K + K * r
  
  #########################
  
  initNoShift <- initializeGaussianPosterior(V = V1, lambda2 = lambda2, X = X, T = TT);
  
  beta0 <- initNoShift$betadelta[1]
  tau2beta0 <- initNoShift$tau2all[1]
  sigma2 <- initNoShift$sigma2
  
  psi <- exp(initNoShift$sigma2)
  
  fit0 <- initNoShift$fit0
  
  tmpbetadelta <- readbetadelta(initNoShift$betadelta, q, p, 0, 0)
  tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, 0, 0)
  
  tmptau2all <- readtau2all(initNoShift$tau2all, q, p, 0, 0)
  tmptau2all <- assignTau2(tmptau2all, q, p, 0, 0)
  
  V2 <- rep(-Inf, length(Y))
  
  #########################
  
  betadelta0out <- matrix(NA, nrow = nsim, ncol = (1 + p + q))
  tau2all0out <- matrix(NA, nrow = nsim, ncol = (1 + p + q))
  sigma2out <- rep(NA, nsim)
  lambda2out <- rep(NA, nsim)
  psiout <- rep(NA, nsim)
  
  fit00out <- matrix(NA, nrow = nsim, ncol = n)
  fit01out <- matrix(NA, nrow = nsim, ncol = n)
  
  #########################
  
  cnt = 1;
  for (sim in 1:(nsim + burnin1)) {
    V1 <- getV1(Y, V1, V2, psi, fit0, sigma2, burnin)
    
    TT <- getT(V1, p)
    
    layer1NBNoShift <- getPosteriorLayer2NoShift(
      V = V1, beta0 = tmpbetadelta$beta0, tau2beta0 = tmptau2all$tau2beta0, 
      sigma2 = sigma2, lambda2 = lambda2, updatelambda2 = updatelambda, 
      burnin = burnin2, nsim = 1,
      X = X, T = TT, beta1 = tmpbetadelta$beta1, beta2 = tmpbetadelta$beta2, 
      tau2beta1 = tmptau2all$tau2beta1, tau2beta2 = tmptau2all$tau2beta2
    )
    
    tmpbetadelta <- readbetadelta(layer1NBNoShift$betadelta, q, p, 0, 0)
    tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, 0, 0)
    
    tmptau2all <- readtau2all(layer1NBNoShift$tau2all, q, p, 0, 0)
    tmptau2all <- assignTau2(tmptau2all, q, p, 0, 0)
    
    sigma2 <- layer1NBNoShift$sigma2
    psi <- getPsi(Y, V1, V2, psi, burnin, c1, c2)
    lambda2 <- layer1NBNoShift$lambda2
    
    fit0 <- layer1NBNoShift$fit0
    
    if (sim > burnin1) {
      betadelta0out[cnt, ] = layer1NBNoShift$betadelta;
      tau2all0out[cnt, ] = layer1NBNoShift$tau2all;
      sigma2out[cnt] <- sigma2
      psiout[cnt] <- psi
      lambda2out[cnt] <- lambda2
      fit00out[cnt, ] <- layer1NBNoShift$fit0
      fit01out[cnt, ] <- layer1NBNoShift$fit1
      cnt = cnt + 1;
    }
  }

  NoShift <- list(
    "betadelta" = betadelta0out,
    "tau2all" = tau2all0out,
    "sigma2" = sigma2out,
    "psi" = psiout,
    "lambda2" = lambda2out,
    "fit0" = fit00out,
    "fit1" = fit01out,
    "m" = m,
    "q" = q,
    "p" = p,
    "K" = K,
    "r" = r
  )
  
  #########################
  
   
}




tol <- 1e-6

nsim <- 1000
burnin <- 100

lambda2 <- 100
K <- 20
p <- 10

n <- 730

spl <- 0.95

Y <- arima.sim(list(ar = 0.5), n = n)
Y[100:190] <- Y[100:190] + 3
#Y[1:50] <- Y[1:50] + 1
#Y[51:100] <- Y[51:100] + 5
#Y[101:150] <- Y[101:150] + 10
#Y[151:200] <- Y[151:200] + 20
#Y[201:250] <- Y[201:250] + 30


aa <- posteriorLayer2(Y, X = NULL, W = NULL, p = p, K = K, lambda = 10,
                            updatelambda = TRUE, 
                            clustering = "kmeans", spl = spl, 
                            changepoint = "idetect", 
                            burnin = burnin, nsim = nsim, tol = 1e-6)


TT <- getT(Y, p)

Uvec <- idetect_rcpp(Y, K - 1)
U = getU(Uvec, n, K)
gamma <- colMeans(U)


init0 <- initializeGaussianPosterior(Y, lambda2, T = TT);
layer20 <- getPosteriorLayer20(Y, init0$betadelta[1], init0$tau2all[1], 
                               init0$sigma2, lambda2, burnin = burnin, nsim = nsim,
                               gamma = gamma,
                               T = TT,
                               beta2 = init0$betadelta[2:6],  
                               tau2beta2 = init0$tau2all[2:6])

init1 <- initializeGaussianPosterior(Y, lambda2, T = TT, U = U);
layer21 <- getPosteriorLayer21(Y, init0$betadelta[1], init0$tau2all[1], 
                             init0$sigma2, lambda2, burnin = burnin, nsim = nsim,
                             gamma = gamma,
                             T = TT, U = U,
                             beta2 = init0$betadelta[2:6], delta0 = init0$betadelta[7:11], 
                             tau2beta2 = init0$tau2all[2:6], tau2delta0 = init0$tau2all[7:11])

 
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




plot(Y, type = 'l')

for (i in 1:nsim) {
  points(layer21$fit1[i, ], col = 'blue')
  points(layer20$fit0[i, ], col = 'red')
}

getPosteriorLayer2(arma::vec V, arma::vec beta0, arma::vec tau2beta0,
                   double sigma2, double lambda2, bool updatelambda2 = true, 
                   int burnin = 50, int nsim = 100,
                   Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> T=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> U=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> W=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> beta1=R_NilValue, 
                   Rcpp::Nullable<Rcpp::NumericVector> beta2=R_NilValue, 
                   Rcpp::Nullable<Rcpp::NumericVector> delta0=R_NilValue, 
                   Rcpp::Nullable<Rcpp::NumericVector> delta1=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> tau2beta1=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> tau2beta2=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> tau2delta0=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> tau2delta1=R_NilValue)