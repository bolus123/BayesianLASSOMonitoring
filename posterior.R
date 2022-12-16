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

posteriorGaussian <- function(Y, X = NULL, W = NULL, p = 5, K = 5, lambda = 10,
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
          sigma2 = layer2Shift$sigma2[nsim], lambda2 = layer2Shift$lambda2[nsim], updatelambda2 = updatelambda, 
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


posteriorNB <- function(Y, X = NULL, W = NULL, p = 5, K = 5, lambda = 10,
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
      Uvec <- idetect_rcpp(V, K - 1)
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
  
  #mu <- exp(initNoShift$fit1)
  #EY <- mean(mu)
  #VarY <- mean((Y - EY) ^ 2)
  #
  #pp <- EY / VarY
  #
  #psi <- EY * (pp / (1 - pp))
  psi <- exp(initNoShift$sigma2)
  
  tmpbetadelta <- readbetadelta(initNoShift$betadelta, q, p, 0, 0)
  tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, 0, 0)
  
  tmptau2all <- readtau2all(initNoShift$tau2all, q, p, 0, 0)
  tmptau2all <- assignTau2(tmptau2all, q, p, 0, 0)
  
  #########################
  
  ## get the model for the no-shift situation
  
  layer1NoShift <- getPosteriorLayer1NBNoShift(Y, V1, initNoShift$fit0, 
                                    tmpbetadelta$beta0, tmptau2all$tau2beta0, 
                              initNoShift$sigma2, psi, lambda2, p,
                              updatelambda2 = updatelambda, updatepsi = updatepsi,
                              c1 = c1, c2 = c2,
                              burnin1 = burnin1, burnin2 = burnin2, nsim = nsim,
                              X = X,
                              beta1=tmpbetadelta$beta1, 
                              beta2=tmpbetadelta$beta2, 
                              tau2beta1=tmptau2all$tau2beta1,
                              tau2beta2=tmptau2all$tau2beta2) 
  
  #########################
  
  if (K > 0) {
    initShift <- initializeGaussianPosterior(V = V1, lambda2 = lambda2, X = X, 
                                             T = TT, U = U, W = W);
    
    psi <- exp(initShift$sigma2)
    
    #EY <- exp(initShift$fit1)
    #VarY <- (Y - EY) ^ 2
    #psi <- mean(EY^2 / (VarY - EY))
    
    tmpbetadelta <- readbetadelta(initShift$betadelta, q, p, K, r)
    tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, K, r)
    
    tmptau2all <- readtau2all(initShift$tau2all, q, p, K, r)
    tmptau2all <- assignTau2(tmptau2all, q, p, K, r)
    
    layer1Shift <- getPosteriorLayer1NBShift(Y, V1, initShift$fit1, 
                                             tmpbetadelta$beta0, tmptau2all$tau2beta0, 
                                             initShift$sigma2, psi, lambda2, p,
                                             updatelambda2 = updatelambda, updatepsi = updatepsi,
                                             c1 = c1, c2 = c2,
                                             burnin1 = burnin1, burnin2 = burnin2, nsim = nsim,
                                             gamma=gamma,
                                             X=X,
                                             U=U,
                                             W=W,
                                             beta1=tmpbetadelta$beta1, 
                                             beta2=tmpbetadelta$beta2, 
                                             delta0=tmpbetadelta$delta0, 
                                             delta1=tmpbetadelta$delta1,
                                             tau2beta1=tmptau2all$tau2beta1,
                                             tau2beta2=tmptau2all$tau2beta2,
                                             tau2delta0=tmptau2all$tau2delta0,
                                             tau2delta1=tmptau2all$tau2delta1)
    
    ## rebuild the estimated time
    
    perc <- rep(0, K + K * r - 2)
    
    start <- (1 + q + p + 1)
    delta <- layer1Shift$betadelta[, start:m]
    
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
        
        tmpbetadelta <- readbetadelta(layer1Shift$betadelta[nsim, ], q, p, K, r)
        tmpbetadelta <- assignBetadelta(tmpbetadelta, q, p, K, r)
        
        tmptau2all <- readtau2all(layer1Shift$tau2all[nsim, ], q, p, K, r)
        tmptau2all <- assignTau2(tmptau2all, q, p, K, r)
        
        
        layer1Shift <- getPosteriorLayer1NBShift(Y, layer1Shift$V[nsim, ], layer1Shift$fit1[nsim, ], 
                                  tmpbetadelta$beta0, tmptau2all$tau2beta0, 
                                  layer1Shift$sigma2[nsim], layer1Shift$psi[nsim], layer1Shift$lambda2[nsim], p,
                                  updatelambda2 = updatelambda, updatepsi = updatepsi,
                                  c1 = c1, c2 = c2,
                                  burnin1 = burnin1, burnin2 = burnin2, nsim = nsim,
                                  gamma=tmpgamma,
                                  X=X,
                                  U=tmpU,
                                  W=W,
                                  beta1=tmpbetadelta$beta1, 
                                  beta2=tmpbetadelta$beta2, 
                                  delta0=tmpbetadelta$delta0, 
                                  delta1=tmpbetadelta$delta1,
                                  tau2beta1=tmptau2all$tau2beta1,
                                  tau2beta2=tmptau2all$tau2beta2,
                                  tau2delta0=tmptau2all$tau2delta0,
                                  tau2delta1=tmptau2all$tau2delta1)
        
      } else {
        break
      }
      
    }
    
    
  }
  
  #########################
  
  out <- list(
    "NoShift" = layer1NoShift,
    "Shift" = layer1Shift
  )
  return(out) 
  
  
}


posteriorZinfNB <- function(Y, X = NULL, W = NULL, p = 5, K = 5, lambda = 10,
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
    TT1 <- getT(V1, p)
  } else {
    TT1 <- NULL
  }
  
  if (K > 0) {
    if (changepoint == "idetect") {
      U1vec <- idetect_rcpp(V1, K - 1)
      U1 <- getU(U1vec, n, K)
      gamma1 <- colMeans(U1)
    }
  } else {
    U1 <- NULL
    gamma1 <- NULL
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
  
  initNoShift1 <- initializeGaussianPosterior(V = V1, lambda2 = lambda2, X = X, T = TT1);
  
  tmpbetadelta1 <- readbetadelta(initNoShift1$betadelta, q, p, 0, 0)
  tmpbetadelta1 <- assignBetadelta(tmpbetadelta1, q, p, 0, 0)
  
  tmptau2all1 <- readtau2all(initNoShift1$tau2all, q, p, 0, 0)
  tmptau2all1 <- assignTau2(tmptau2all1, q, p, 0, 0)
  
  #########################
  
  V2 <- as.numeric(Y == 0)

  if (p > 0) {
    TT2 <- getT(V2, p)
  } else {
    TT2 <- NULL
  }
  
  if (K > 0) {
    if (changepoint == "idetect") {
      U2vec <- idetect_rcpp(V2, K - 1)
      U2 <- getU(U2vec, n, K)
      gamma1 <- colMeans(U2)
    }
  } else {
    U2 <- NULL
    gamma2 <- NULL
  }
  
  U2vec <- idetect_rcpp(V2, K)
  U2 <- getU(U2vec, n, K)
  gamma2 <- colMeans(U2)
  
  dat <- as.data.frame(cbind(V2, TT2, U2))
  names(dat)[1] <- 'VV'
  
  init2 <- glm(VV~., data = dat, family = binomial())
  
  betadelta2 <- init2$coefficients
  betadelta2[is.na(betadelta2)] <- tol
  sigma22 <- sd(init2$residuals)
  tau2all2 <- getTau2(betadelta2, sigma22, lambda2)
  
  tmpbetadelta2 <- readbetadelta(betadelta2, q, p, 0, 0)
  tmpbetadelta2 <- assignBetadelta(tmpbetadelta2, q, p, 0, 0)
  
  tmptau2all2 <- readtau2all(tau2all2, q, p, 0, 0)
  tmptau2all2 <- assignTau2(tmptau2all2, q, p, 0, 0)
  
  V2 <- rnorm(n, init2$fitted.values, sqrt(sigma22))
  
  fit2 <- init2$fitted.values
  
  #########################
  
  #mu <- mean(exp(initNoShift1$fit1))
  #pi <- mean(1 / (1 + exp(-fit2)))
  
  #EY <- (1 - pi) * mu
  #VarY <- mean((Y - EY) ^ 2)
  
  #psi <- 1 / ((VarY / EY - 1) / mu - pi)
  
  psi <- (1 - 1 / (1 + exp(-sigma22))) * exp(initNoShift1$sigma2)
  
  #########################
  
  ## get the model for the no-shift situation
  
  layer1NoShift <- getPosteriorLayer1ZinfNBNoShift(Y = Y, 
                                               V1 = V1, V1hat = initNoShift1$fit1, 
                                               V2 = V2, V2hat = fit2, 
                                               beta01 = tmpbetadelta1$beta0, tau2beta01 = tmptau2all1$tau2beta0,
                                               beta02 = tmpbetadelta2$beta0, tau2beta02 = tmptau2all2$tau2beta0, 
                                               sigma21 = initNoShift1$sigma2, sigma22 = sigma22,
                                               psi = psi, lambda2 = lambda2, p = p, 
                                               updatelambda2 = updatelambda, updatepsi = updatepsi,
                                               c1 = c1, c2 = c2,
                                               burnin1 = burnin1, burnin2 = burnin2, nsim = nsim, 
                                               X = X, 
                                               beta11 = tmpbetadelta1$beta1, 
                                               beta21 = tmpbetadelta1$beta2, 
                                               tau2beta11 = tmptau2all1$tau2beta1,
                                               tau2beta21 = tmptau2all2$tau2beta2,
                                               beta12 = tmpbetadelta2$beta1,
                                               beta22 = tmpbetadelta2$beta2,
                                               tau2beta12 = tmptau2all2$tau2beta1,
                                               tau2beta22 = tmptau2all2$tau2beta2) 
  
  
  #########################
  
  if (K > 0) {
    
    initShift1 <- initializeGaussianPosterior(V = V1, lambda2 = lambda2, X = X, 
                                             T = TT1, U = U1, W = W);
    
    initShift2 <- initializeGaussianPosterior(V = V2, lambda2 = lambda2, X = X, 
                                              T = TT2, U = U2, W = W);
    
    psi <- (1 - 1 / (1 + exp(-initShift2$sigma2))) * exp(initShift1$sigma2)
    
    tmpbetadelta1 <- readbetadelta(initShift1$betadelta, q, p, K, r)
    tmpbetadelta1 <- assignBetadelta(tmpbetadelta1, q, p, K, r)
    
    tmptau2all1 <- readtau2all(initShift1$tau2all, q, p, K, r)
    tmptau2all1 <- assignTau2(tmptau2all1, q, p, K, r)
    
    tmpbetadelta2 <- readbetadelta(initShift2$betadelta, q, p, K, r)
    tmpbetadelta2 <- assignBetadelta(tmpbetadelta2, q, p, K, r)
    
    tmptau2all2 <- readtau2all(initShift2$tau2all, q, p, K, r)
    tmptau2all2 <- assignTau2(tmptau2all2, q, p, K, r)
    
    layer1Shift <- getPosteriorLayer1ZinfNBShift(Y, 
                                             V1, initShift1$fit1, 
                                             V2, initShift2$fit1,
                                             tmpbetadelta1$beta0, tmptau2all1$tau2beta0, 
                                             tmpbetadelta2$beta0, tmptau2all2$tau2beta0, 
                                             initShift1$sigma2, initShift2$sigma2, 
                                             psi, lambda2, p, 
                                             updatelambda, updatepsi, 
                                             c1, c2, 
                                             burnin1, burnin2, nsim, 
                                             gamma1, gamma2, 
                                             X, U1, U2, 
                                             W, 
                                             tmpbetadelta1$beta1, tmpbetadelta1$beta2, 
                                             tmptau2all1$tau2beta1, tmptau2all1$tau2beta2, 
                                             tmpbetadelta2$beta1, tmpbetadelta2$beta2, 
                                             tmptau2all2$tau2beta1, tmptau2all2$tau2beta2, 
                                             tmpbetadelta1$delta0, tmpbetadelta1$delta1, 
                                             tmpbetadelta2$delta0, tmpbetadelta2$delta1, 
                                             tmptau2all1$tau2delta0, tmptau2all1$tau2delta1, 
                                             tmptau2all2$tau2delta0, tmptau2all2$tau2delta1)
    
  
    
    ## rebuild the estimated time
    
    perc1 <- rep(0, K + K * r - 2)
    perc2 <- rep(0, K + K * r - 2)
    
    start <- (1 + q + p + 1)
    delta1 <- layer1Shift$betadelta1[, start:m]
    delta2 <- layer1Shift$betadelta2[, start:m]
    
    if (clustering == "kmeans") {
      #tmpU <- getUMaxProb(layer2Shift$prob)
      
      kk <- 0
      for (i in 2:(K + K * r - 1)) {
        kk <- kk + 1
        cluster1 <- kmeans(t(delta1), centers = i)
        perc1[kk] <- cluster1$betweenss / cluster1$totss
      }
      
      kk <- 0
      for (i in 2:(K + K * r - 1)) {
        kk <- kk + 1
        cluster2 <- kmeans(t(delta2), centers = i)
        perc2[kk] <- cluster2$betweenss / cluster2$totss
      }
      
      #cat(perc, '\n')
      
      target1 <- which(perc1 > spl)
      target2 <- which(perc2 > spl)
      
      if ((length(target1) > 0) | (length(target2) > 0)) {
        
        if (length(target1) > 0) {
          cluster1 <- kmeans(t(delta1), centers = target1[1] + 1)
          maxgroup1 <- max(cluster1$cluster)
          tmpU1 <- matrix(0, nrow = n, ncol = K)
          
          #cat(maxgroup, '\n')
          
          for (i in 1:maxgroup1) {
            tmptarget1 <- which(cluster1$cluster == i)
            for (j in tmptarget1) {
              tmpU1[which(U1[, j] == 1), tmptarget1[1]] <- 1
            }
          }
          
          tmpgamma1 <- colMeans(tmpU1)
          tmpgamma1[which(tmpgamma1 == 0)] <- tol
          tmpgamma1 <- tmpgamma1 / sum(tmpgamma1)
          
        } else {
          tmpU1 <- U1;
        }
        
        if (length(target2) > 0) {
          cluster2 <- kmeans(t(delta2), centers = target2[1] + 1)
          maxgroup2 <- max(cluster2$cluster)
          tmpU2 <- matrix(0, nrow = n, ncol = K)
          
          #cat(maxgroup, '\n')
          
          for (i in 1:maxgroup2) {
            tmptarget2 <- which(cluster2$cluster == i)
            for (j in tmptarget2) {
              tmpU2[which(U2[, j] == 1), tmptarget2[1]] <- 1
            }
          }
          
          tmpgamma2 <- colMeans(tmpU2)
          tmpgamma2[which(tmpgamma2 == 0)] <- tol
          tmpgamma2 <- tmpgamma2 / sum(tmpgamma2)
          
        } else {
          tmpU2 <- U2;
        }
        
        
        
        tmpbetadelta1 <- readbetadelta(layer1Shift$betadelta1[nsim, ], q, p, K, r)
        tmpbetadelta1 <- assignBetadelta(tmpbetadelta1, q, p, K, r)
        
        tmptau2all1 <- readtau2all(layer1Shift$tau2all1[nsim, ], q, p, K, r)
        tmptau2all1 <- assignTau2(tmptau2all1, q, p, K, r)
        
        tmpbetadelta2 <- readbetadelta(layer1Shift$betadelta2[nsim, ], q, p, K, r)
        tmpbetadelta2 <- assignBetadelta(tmpbetadelta2, q, p, K, r)
        
        tmptau2all2 <- readtau2all(layer1Shift$tau2all2[nsim, ], q, p, K, r)
        tmptau2all2 <- assignTau2(tmptau2all2, q, p, K, r)
        
        
        layer1Shift <- getPosteriorLayer1ZinfNBShift(Y, 
                                                     layer1Shift$V1[nsim, ], layer1Shift$fit11[nsim, ], 
                                                     layer1Shift$V2[nsim, ], layer1Shift$fit12[nsim, ],
                                                     tmpbetadelta1$beta0, tmptau2all1$tau2beta0, 
                                                     tmpbetadelta2$beta0, tmptau2all2$tau2beta0, 
                                                     layer1Shift$sigma21[nsim], layer1Shift$sigma22[nsim], 
                                                     layer1Shift$psi[nsim], layer1Shift$lambda2[nsim], p, 
                                                     updatelambda, updatepsi, 
                                                     c1, c2, 
                                                     burnin1, burnin2, nsim, 
                                                     tmpgamma1, tmpgamma2, 
                                                     X, tmpU1, tmpU2, 
                                                     W, 
                                                     tmpbetadelta1$beta1, tmpbetadelta1$beta2, 
                                                     tmptau2all1$tau2beta1, tmptau2all1$tau2beta2, 
                                                     tmpbetadelta2$beta1, tmpbetadelta2$beta2, 
                                                     tmptau2all2$tau2beta1, tmptau2all2$tau2beta2, 
                                                     tmpbetadelta1$delta0, tmpbetadelta1$delta1, 
                                                     tmpbetadelta2$delta0, tmpbetadelta2$delta1, 
                                                     tmptau2all1$tau2delta0, tmptau2all1$tau2delta1, 
                                                     tmptau2all2$tau2delta0, tmptau2all2$tau2delta1)
        
      } else {
        break
      }
      
    }
    
    
  }
  
  #########################
  
  out <- list(
    "NoShift" = layer1NoShift,
    "Shift" = layer1Shift
  )
  return(out) 
  
  
}




tol <- 1e-6

nsim <- 1000
burnin <- 100

lambda2 <- 100
K <- 20
p <- 10

n <- 730

spl <- 0.7

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



plot(Y)

for (i in 1:nsim) {
  points(aa$Shift$fit1[i, ], col = 'blue')
  points(aa$NoShift$fit0[i, ], col = 'red')
}

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