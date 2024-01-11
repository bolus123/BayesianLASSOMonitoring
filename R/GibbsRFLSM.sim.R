#' Get a simulation transformed data using Random Flexible Level Shift Model in the retrospective phase under H0
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param Y1 is the transformed vector
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param X is the matrix of X
#' @param Beta is the matrix of coefficients
#' @param Kappa is the matrix of triggering the corresponding beta
#' @param H is the matrix of X
#' @param Gamma is the matrix of coefficients
#' @param Tau is the matrix of triggering the corresponding beta
#' 
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
GibbsRFLSM.sim.ph1 <- function(nsim, Y1, Phi, muq, sigma2, 
                                  X = NULL, Beta = NULL, Kappa = NULL, 
                                  H = NULL, Gamma = NULL, Tau = NULL) {
  
  TT <- length(Y1)
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  Y1.sim <- matrix(Y1, nrow = TT, ncol = nsim)
  fit.sim <- matrix(NA, nrow = TT, ncol = nsim)
  muX <- matrix(0, nrow = TT, ncol = m)
  
  if (!is.null(X)) {
    muX <- X %*% (Beta * Kappa)
  }
  
  muH <- matrix(0, nrow = TT, ncol = m)
  if (!is.null(H)) {
    muH <- H %*% (Gamma * Tau)
  }
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpmu0 <- muq[tmpsel] + muX[, tmpsel] + muH[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    for (i in (q + 1):TT) {
      fit.sim[i, j] <- tmpmu0[i] + (Y1[(i - 1):(i - q)] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi
      Y1.sim[i, j] <- fit.sim[i, j] + rnorm(1, mean = 0, sd = sqrt(tmpsigma2))
    }
    
  }
  
  list("fit" = fit.sim[-c(1:q), ], "Y.tr" = Y1.sim[-c(1:q), ])
  
}

ewma.filter <- function (x, ratio) {
  c(filter(x * ratio, 1 - ratio, "recursive", init = x[1]))
}

#' Get a simulation transformed data using Random Flexible Level Shift Model in the retrospective phase under H0
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param Y1 is the transformed vector
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param X is the matrix of X
#' @param Beta is the matrix of coefficients
#' @param Kappa is the matrix of triggering the corresponding beta
#' @param H is the matrix of X
#' @param Gamma is the matrix of coefficients
#' @param Tau is the matrix of triggering the corresponding beta
#' 
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
GibbsRFLSM.sim.ph2 <- function(h, nsim, Phi, muq, sigma2, 
                               X = NULL, Beta = NULL, Kappa = NULL, 
                               H = NULL, Gamma = NULL, Tau = NULL, 
                               Y1 = rep(0, dim(Phi)[1]), X1 = NULL, H1 = NULL) {
  
  TT <- length(Y1)
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  Y1.sim <- matrix(Y1, nrow = TT, ncol = nsim)
  
  muX <- matrix(0, nrow = h, ncol = m)
  if (!is.null(X)) {
    muX <- X %*% (Beta * Kappa)
  }
  
  muX1 <- matrix(0, nrow = q, ncol = m)
  if (!is.null(X1)) {
    muX1 <- X1 %*% (Beta * Kappa)
  }
  muX <- rbind(muX1, muX)
  
  muH <- matrix(0, nrow = h, ncol = m)
  if (!is.null(H)) {
    muH <- H %*% (Gamma * Tau)
  }
  
  muH1 <- matrix(0, nrow = q, ncol = m)
  if (!is.null(H1)) {
    muH1 <- H1 %*% (Gamma * Tau)
  }
  muH <- rbind(muH1, muH)
  
  Y2.sim <- matrix(0, nrow = h, ncol = nsim)
  Y2.sim <- rbind(Y1.sim, Y2.sim)
  fit.sim <- matrix(NA, nrow = h + TT, ncol = nsim)
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpmu0 <- muq[tmpsel] + muX[, tmpsel] + muH[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    tmperr <- rnorm(TT + h, mean = 0, sd = sqrt(tmpsigma2))
   
    for (i in (q + 1):(TT + h)) {
      fit.sim[i, j] <- tmpmu0[i] + (Y2.sim[(i - 1):(i - q), j] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi
      Y2.sim[i, j] <- fit.sim[i, j] + tmperr[i]
    }
    
  }
  
  list("fit" = fit.sim[-c(1:q), ], "Y.tr" = Y2.sim[-c(1:q), ])
  
}

GibbsRFLSM.loglik <- function(Y.tr, Phi, muq, sigma2, 
                              X = NULL, Beta = NULL, Kappa = NULL, 
                              H = NULL, Gamma = NULL, Tau = NULL, 
                              Y1 = rep(0, dim(Phi)[1]), X1 = NULL, H1 = NULL) {
  
  YY <- c(Y1, Y.tr)
  TT <- length(YY)
  
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  
  muX <- matrix(0, nrow = TT - q, ncol = m)
  if (!is.null(X)) {
    muX <- X %*% (Beta * Kappa)
  }
  
  muX1 <- matrix(0, nrow = q, ncol = m)
  if (!is.null(X1)) {
    muX1 <- X1 %*% (Beta * Kappa)
  }
  muX <- rbind(muX1, muX)
  
  muH <- matrix(0, nrow = TT - q, ncol = m)
  if (!is.null(H)) {
    muH <- H %*% (Gamma * Tau)
  }
  
  muH1 <- matrix(0, nrow = q, ncol = m)
  if (!is.null(H1)) {
    muH1 <- H1 %*% (Gamma * Tau)
  }
  muH <- rbind(muH1, muH)
  fit.sim <- matrix(NA, nrow = TT - q, ncol = m)
  loglik.sim <- matrix(NA, nrow = TT - q, ncol = m)
  
  for (j in 1:m) {
    tmpPhi <- Phi[, j]
    tmpmu0 <- muq[j] + muX[, j] + muH[, j]
    tmpsigma2 <- sigma2[j]
    
    for (i in (q + 1):(TT)) {
      fit.sim[i - q, j] <- tmpmu0[i] + (YY[(i - 1):(i - q)] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi
      loglik.sim[i - q, j] <- dnorm(YY[i], fit.sim[i - q, j], sqrt(tmpsigma2), log = TRUE)
    }
    
  }
  
  rowMeans(loglik.sim)
  
}


