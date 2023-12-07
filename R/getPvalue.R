#' obtain the fits under H0
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param muq is the mean
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
#' Fit0(Y, result$Phi, result$muq)
#' 
fit.ph1.H0 <- function(Y, Phi, muq, sigma2, 
                       X = NULL, Beta = NULL, Kappa = NULL,
                       bt = TRUE, log = TRUE, const = 1, sta = TRUE, 
                       meanY = 0, sdY = 1, method = 'median') {
  
  TT <- length(Y)
  q <- dim(Phi)[1]
  nsim <- dim(Phi)[2]
  ff <- matrix(NA, nrow = TT - q, ncol = nsim)
  
  if (method == "median") {
    Phihat <- apply(Phi, 1, median)
    muqhat <- median(muq)
    sigma2hat <- median(sigma2)
  }
  
  muXhat <- rep(0, TT)
  if (!is.null(NULL)) {
    if (method == "median") {
      Betahat <- apply(Beta, 1, median)
      Kappahat <- apply(Kappa, 1, median)
    }
    muXhat <- X %*% (Betahat * Kappahat)
  }
  
  ff <- rep(muqhat, TT)
  
  for (i in (q + 1):TT) {
    ff[i] <- muqhat + muXhat[i] + 
      (Y[(i - 1):(i - q)] - muqhat - muXhat[(i - 1):(i - q)] ) %*% Phihat 
  }
  
  out <- ff
  
  if (bt == TRUE) {
    out <- backtrans(out, log, const, sta, meanY, sdY)
  }
  
  out
  
}

#' simulates the time series using Draws from MCMC in Phase I under H0
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param muq is the mean
#' @param sigma2 is the variance of errors
#' @param logcc is the log
#' @param standardization is the standardization
#' @param nsim is the nsim
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
#'
#' GibbsRFLSM.sim0.ph1(Y, result$Phi, result$muq, result$sigma2) 
#' 
GibbsRFLSM.sim0.ph1 <- function(Y, Phi, muq, sigma2, logcc = FALSE, 
                           standardization = FALSE, nsim = 1000) {
  
  TT <- length(Y)
  q <- dim(Phi)[1]
  nsim.par <- length(muq)
  
  sim <- matrix(NA, nrow = TT, ncol = nsim)
  
  Y.tr <- Y
  
  if (logcc == TRUE) {
    Y.tr <- log(Y.tr + 1) 
  }
  
  if (standardization == TRUE) {
    meanY <- mean(Y.tr)
    sdY <- sd(Y.tr)
    Y.tr <- (Y.tr - meanY) / sdY
  }
  
  for (i in 1:nsim) {
    sel <- sample(1:nsim.par, 1)
    tmpPhi <- Phi[, sel]
    tmpmuq <- muq[sel]
    tmpsigma2 <- sigma2[sel]
    
    fit0 <- Fit0(Y.tr, matrix(tmpPhi, ncol = 1), matrix(tmpmuq, ncol = 1))
    tmp <- fit0 + rnorm(TT, 0, sqrt(tmpsigma2))
  
    sim[, i] <- tmp
  }
  
  sim0 <- sim
  
  if (standardization == TRUE) {
    sim <- sim * sdY + meanY
  }
  
  if (logcc == TRUE) {
    sim <- exp(sim) - 1
  }
  
  list("Y.tr" = sim0, "Y" = sim)
  
}

#' simulates the time series using Draws from MCMC in Phase II under H0
#' 
#' @param n is the number of time series
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param muq is the mean
#' @param sigma2 is the variance of errors
#' @param logcc is the log
#' @param standardization is the standardization
#' @param nsim is the nsim
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
#' GibbsRFLSM.sim0.ph2(Y, result$Phi, result$muq, result$sigma2) 
#' 
GibbsRFLSM.sim0.ph2 <- function(n, Phi, muq, sigma2, Y0, logcc = FALSE, 
                                standardization = FALSE, nsim = 1000) {
  
  TT <- length(Y0)
  q <- dim(Phi)[1]
  nsim.par <- length(muq)
  
  meanY <- mean(Y0)
  sdY <- sd(Y0)
  
  sim.tr.ph1 <- matrix(NA, nrow = TT, ncol = nsim)
  sim.ph1 <- matrix(NA, nrow = TT, ncol = nsim)
  
  sim <- matrix(NA, nrow = n, ncol = nsim)

  tmp <- rep(NA, n + q)
  
  for (i in 1:nsim) {
    
    sel <- sample(1:nsim.par, 1)
    tmpPhi <- Phi[, sel]
    tmpmuq <- muq[sel]
    tmpsigma2 <- sigma2[sel]
    
    Y0.sim <- GibbsRFLSM.sim0.ph1(Y0, 
                                  matrix(tmpPhi, ncol = 1), 
                                  matrix(tmpmuq, ncol = 1), 
                                  matrix(tmpsigma2, ncol = 1), 
                                  logcc, standardization, 1)
    
    sim.tr.ph1[, i] <- Y0.sim$Y.tr
    sim.ph1[, i] <- Y0.sim$Y
    
    tmp[1:q] <- Y0.sim$Y.tr[(TT - q + 1):TT, 1]
    
    for (j in (q + 1):(n + q)) {
      tmp[j] <- tmpmuq + (tmp[(j - 1):(j - q)] - tmpmuq) %*% 
        tmpPhi + rnorm(1, 0, sqrt(tmpsigma2))
    }
    sim[, i] <- tmp[(q + 1):(n + q)]
  }
  
  
  sim0 <- sim
  
  if (standardization == TRUE) {
    sim <- sim * sdY + meanY
  }
  
  if (logcc == TRUE) {
    sim <- exp(sim) - 1
  }
  
  list("Y1.tr" = sim.tr.ph1, "Y1" = sim.ph1, "Y2.tr" = sim0, "Y2" = sim)
  
}

#' simulates the maximums of residuals
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Muhat is the mean
#' @param sigma2hat is the variance of errors
#' @param nsim is the number of simulations
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
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' Muhat <- rep(NA, T)
#' for (i in 1:T) {
#'   Muhat[i] <- median(result$Mu[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#' 
#' GibbsRFLSM.simmax.residual(Y, result$Phi, result$Mu, result$sigma2, 
#' Phihat, Muhat, sigma2hat)
#' 
GibbsRFLSM.max.residual <- function(Y.tr, Phihat, muqhat, sigma2hat) {
  
  q <- length(Phihat)
  
  fit0hat <- Fit0(Y.tr, matrix(Phihat, ncol = 1), matrix(muqhat, ncol = 1))
    
  tmp <- (Y.tr[-c(1:q)] - fit0hat[-c(1:q), ]) ^ 2 / sigma2hat
    
  maxZ2 <- max(tmp)
  
  list("maxZ2" = maxZ2, "Z2" = tmp)
  
}

#' simulates the maximums of residuals
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Muhat is the mean
#' @param sigma2hat is the variance of errors
#' @param nsim is the number of simulations
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
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' Muhat <- rep(NA, T)
#' for (i in 1:T) {
#'   Muhat[i] <- median(result$Mu[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#' 
#' GibbsRFLSM.simmax.residual(Y, result$Phi, result$Mu, result$sigma2, 
#' Phihat, Muhat, sigma2hat)
#' 
GibbsRFLSM.max.residual.sim <- function(Y, Phi, muq, sigma2, 
                                       Phihat, muqhat, sigma2hat, 
                                       logcc = FALSE,
                                       standardization = FALSE,
                                       nsim = 1000) {
  
  TT <- length(Y)
  q <- dim(Phi)[1]
  
  Y1 <- GibbsRFLSM.sim0.ph1(Y, Phi, muq, sigma2, logcc, 
                                  standardization, nsim)
  
  maxZ2 <- rep(NA, nsim)
  
  for (i in seq(nsim)) {
    
    tmp <- GibbsRFLSM.max.residual(Y1$Y.tr[, i], Phihat, muqhat, sigma2hat)
  
    maxZ2[i] <- tmp$maxZ2
  }
  
  maxZ2
  
}

#' simulates the maximums of retrospective scan statistics
#' 
#' @param Y is a vector
#' @param N is the sum of all counts
#' @param w is the moving window
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Muhat is the mean
#' @param nsim is the number of simulations
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
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' Muhat <- rep(NA, T)
#' for (i in 1:T) {
#'   Muhat[i] <- median(result$Mu[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#' 
#' GibbsRFLSM.simmax.residual(Y, result$Phi, result$Mu, result$sigma2, 
#' Phihat, Muhat)
#' 
GibbsRFLSM.max.scan <- function(Y, w,
                                Phihat, muqhat,
                                logcc = FALSE, standardization = FALSE, tol = 1e-10) {
  
  
  
  Y.tr <- Y
  
  if (logcc == TRUE) {
    Y.tr <- log(Y + 1)
  }
  
  if (standardization == TRUE) {
    meanY <- mean(Y.tr)
    sdY <- sd(Y.tr)
    Y.tr <- (Y.tr - meanY) / sdY
  }
  
  fit0hat <- Fit0(Y.tr, matrix(Phihat, ncol = 1), matrix(muqhat, ncol = 1))
  
  if (standardization == TRUE) {
    fit0hat <- fit0hat * sdY + meanY
  }
  
  if (logcc == TRUE) {
    fit0hat <- exp(fit0hat) - 1
  }
  
  fit0hat[fit0hat < 0] <- tol
  
  
  U <- w * Y * log(Y / fit0hat) - w * (Y - fit0hat)
  
  for (j in seq(n)) {
    U[j] <- ifelse(Y[j] > fit0hat[j], U[j], 0)
  }
  
  maxU <- max(U)
  
  list("maxU" = maxU, "U" = U)
  
}


#' simulates the maximums of retrospective scan statistics
#' 
#' @param Y is a vector
#' @param N is the sum of all counts
#' @param w is the moving window
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Muhat is the mean
#' @param nsim is the number of simulations
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
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' Muhat <- rep(NA, T)
#' for (i in 1:T) {
#'   Muhat[i] <- median(result$Mu[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#' 
#' GibbsRFLSM.simmax.residual(Y, result$Phi, result$Mu, result$sigma2, 
#' Phihat, Muhat)
#' 
GibbsRFLSM.max.scan.sim <- function(Y, Phi, muq, sigma2, w,
                                  Phihat, muqhat, 
                                  logcc = FALSE, standardization = FALSE, nsim = 1000, tol = 1e-10) {
  
  
  TT <- length(Y)
  q <- dim(Phi)[1]
  
  Y1 <- GibbsRFLSM.sim0.ph1(Y, Phi, muq, sigma2, logcc, 
                            standardization, nsim)
  
  maxU <- rep(NA, nsim)
  for (i in seq(nsim)) {
    
    tmp <- GibbsRFLSM.max.scan(Y1$Y[, i], w,
                               Phihat, muqhat,
                               logcc, standardization, tol)
    
    maxU[i] <- tmp$maxU
  }
  
  maxU
  
}

#' simulates posterior predictive p values using for the Phase I chart
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu0 is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Mu0hat is the mean
#' @param sigma2hat is the variance of errors
#' @param nsim is the number of simulations
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
#' Mu0 <- matrix(NA, nrow = T, ncol = nsim)
#' for (j in 1:nsim) {
#'   Mu0[, j] <- result$muq[j]
#' }
#' 
#' Mu0hat <- rep(NA, T)
#' for (i in 1:T) {
#'   Mu0hat[i] <- median(Mu0[i, ])
#' }
#' 
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#'
#' GibbsRFLSM.PPP.residual(Y, result$Phi, Mu0, result$sigma2, 
#' Phihat, Mu0hat, sigma2hat)
#' 
GibbsRFLSM.PPP.residual <- function(Y, Phi, Mu0, sigma2, 
                                    Phihat, Mu0hat, sigma2hat, 
                                    nsim = 1000) {
  
  q <- dim(Phi)[1]
  n <- length(Y)
  m <- dim(Phi)[2]
  
  ccrep <- GibbsRFLSM.simmax.residual(Y, Phi, Mu0, sigma2, 
                                      Phihat, Mu0hat, sigma2hat, 
                                      nsim)
  
  tmp <- Y - Mu0hat
  tmpV <- getV(tmp, q)
  tmp <- tmp[-c(1:q)]
  tmpV <- tmpV[-c(1:q), ]
  
  tmpresi <- (tmp - tmpV %*% Phihat) ^ 2 / sigma2hat
  
  tmpOmni <- mean(ccrep > max(tmpresi))
  tmpInd <- rep(NA, n - q)
  for (i in 1:(n - q)) {
    tmpInd[i] <- mean(ccrep > tmpresi[i])
  }
  list("Omni" = tmpOmni, "Ind" = tmpInd, "cs" = max(tmpresi), 'ref' = ccrep) 
}


#' simulates posterior predictive p values using for the Phase I chart
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu0 is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Mu0hat is the mean
#' @param sigma2hat is the variance of errors
#' @param nsim is the number of simulations
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
#' Mu0 <- matrix(NA, nrow = T, ncol = nsim)
#' for (j in 1:nsim) {
#'   Mu0[, j] <- result$muq[j]
#' }
#' 
#' Mu0hat <- rep(NA, T)
#' for (i in 1:T) {
#'   Mu0hat[i] <- median(Mu0[i, ])
#' }
#' 
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#'
#' GibbsRFLSM.PPP.scan(Y, result$Phi, Mu0, result$sigma2, 
#' Phihat, Mu0hat, sigma2hat)
#' 
GibbsRFLSM.PPP.scan <- function(Y, N, w, Phi, Mu0, sigma2, 
                              Phihat, Mu0hat, nsim = 1000,
                              logcc = FALSE, standardization = FALSE, 
                              meanY = NULL, sdY = NULL) {
  
  q <- dim(Phi)[1]
  n <- length(Y)
  m <- dim(Phi)[2]
  
  ccrep <- GibbsRFLSM.simmax.scan(Y, N, w, Phi, Mu0, sigma2, 
                                Phihat, Mu0hat, 
                                nsim, logcc, standardization, meanY, sdY)
  
  tmp <- Y - Mu0hat
  tmpV <- getV(tmp, q)
  tmp <- tmp[-c(1:q)]
  tmpV <- tmpV[-c(1:q), ]
  tmpY <- Mu0hat[-c(1:q)] + tmpV %*% Phihat
  
  if (standardization == TRUE) {
    Yr <- Y * sdY + meanY
    tmpY <- tmpY * sdY + meanY
  }
  
  if (logcc == TRUE) {
    Yr <- exp(Yr) - 1
    tmpY <- exp(tmpY) - 1
  }
  
  U <- w * Yr * log(Yr / tmpY) + 
    (N - w * Yr) * log((N - w * Yr) / (N - w * tmpY))
  
  tmpOmni <- mean(ccrep > max(U))
  tmpInd <- rep(NA, n - q)
  for (i in 1:(n - q)) {
    tmpInd[i] <- mean(ccrep > U[i])
  }
  list("Omni" = tmpOmni, "Ind" = tmpInd, "cs" = max(U), 'ref' = ccrep) 
}


#' obtains the Phase I charting constant using posterior predictive p value
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu0 is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Mu0hat is the mean
#' @param sigma2hat is the variance of errors
#' @param FAP0 is the false alarm probability
#' @param nsim is the number of simulations
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
#' Mu0 <- matrix(NA, nrow = T, ncol = nsim)
#' for (j in 1:nsim) {
#'   Mu0[, j] <- result$muq[j]
#' }
#' 
#' Mu0hat <- rep(NA, T)
#' for (i in 1:T) {
#'   Mu0hat[i] <- median(Mu0[i, ])
#' }
#' 
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#'
#' GibbsRFLSM.CC.PPP.residual(Y, result$Phi, Mu0, result$sigma2, 
#' Phihat, Mu0hat, sigma2hat)
#' 
#' 
GibbsRFLSM.CC.PPP.residual <- function(Y, Phi, Mu0, sigma2, 
  Phihat, Mu0hat, sigma2hat, FAP0 = 0.2, nsim = 1000) {
  
  res <- GibbsRFLSM.PPP.residual(Y, Phi, Mu0, sigma2, 
                                      Phihat, Mu0hat, sigma2hat, nsim) 
  
  cc2 <- quantile(res$ref, 1 - FAP0)
  cc <- sqrt(cc2)
  out <- list("cc" = cc,  "Omni" = res$Omni, "Ind" = res$Ind, 
              "cs" = res$cs, 'ref' = res$ref)
  return(out)
  
}


#' obtains the Phase I charting constant using posterior predictive p value
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu0 is the mean
#' @param sigma2 is the variance of errors
#' @param Phihat is the coefficient
#' @param Mu0hat is the mean
#' @param sigma2hat is the variance of errors
#' @param FAP0 is the false alarm probability
#' @param nsim is the number of simulations
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
#' Mu0 <- matrix(NA, nrow = T, ncol = nsim)
#' for (j in 1:nsim) {
#'   Mu0[, j] <- result$muq[j]
#' }
#' 
#' Mu0hat <- rep(NA, T)
#' for (i in 1:T) {
#'   Mu0hat[i] <- median(Mu0[i, ])
#' }
#' 
#' Phihat <- rep(NA, q)
#' for (i in 1:q) {
#'   Phihat[i] <- median(result$Phi[i, ])
#' }
#' 
#' sigma2hat <- median(result$sigma2)
#'
#' GibbsRFLSM.CC.PPP.residual(Y, result$Phi, Mu0, result$sigma2, 
#' Phihat, Mu0hat, sigma2hat)
#' 
#' 
GibbsRFLSM.CC.PPP.scan <- function(Y, Phi, Mu0, sigma2, 
                                  Phihat, Mu0hat, sigma2hat, FAP0 = 0.2, nsim = 1000,
                                  logcc = FALSE, standardization = FALSE, 
                                  meanY = NULL, sdY = NULL) {
  
  res <- GibbsRFLSM.PPP.scan(Y, Phi, Mu0, sigma2, 
                             Phihat, Mu0hat, sigma2hat, nsim, 
                             logcc, standardization, 
                             meanY, sdY)
  
  
  cc <- quantile(res$ref, 1 - FAP0)
  out <- list("cc" = cc,  "Omni" = res$Omni, "Ind" = res$Ind, 
              "cs" = res$cs, 'ref' = res$ref)
  return(out)
  
}