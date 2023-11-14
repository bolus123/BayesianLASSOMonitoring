#' simulates the time series using Draws from MCMC
#' 
#' @param Y is a vector
#' @param Phihat is the coefficient
#' @param Muhat is the mean
#' @param sigma2hat is the variance of errors
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
#' Phihat <- result$Phi[, 1]
#' Muhat <- result$Mu[, 1]
#' sigma2hat <- result$sigma2[1]
#'
#' GibbsRFLSM.sim(Y, Phihat, Muhat, sigma2hat) 
#' 
GibbsRFLSM.sim <- function(Y, Phihat, Muhat, sigma2hat) {
  
  TT <- length(Y)
  q <- length(Phihat)
  
  sim <- rep(NA, TT)
  
  sim[1:q] <- Y[1:q]
  
  for (ii in (q + 1):TT) {
    sim[ii] <- Muhat[ii] + (sim[(ii - 1):(ii - q)] - Muhat[(ii - 1):(ii - q)]) %*% 
      Phihat + rnorm(1, 0, sqrt(sigma2hat))
  }
  
  sim 
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
GibbsRFLSM.simmax.residual <- function(Y, Phi, Mu, sigma2, 
                                       Phihat, Muhat, sigma2hat, 
                                       nsim = 1000) {
  
  q <- dim(Phi)[1]
  n <- length(Y) - q
  out <- rep(NA, nsim)
  m <- dim(Phi)[2]
  
  for (i in seq(nsim)) {
    
    k <- sample(1:m, 1)
    tmpPhi <- Phi[, k]
    tmpMu <- Mu[, k]
    tmpsigma2 <- sigma2[k]
    
    tmp <- GibbsRFLSM.sim(Y, tmpPhi, tmpMu, tmpsigma2)
    tmpV <- tmp - Muhat
    tmpVas <- getV(tmpV, q)
    tmpV <- tmpV[-c(1:q)]
    tmpVas <- tmpVas[-c(1:q), ]
    tmp <- (tmpV - tmpVas %*% Phihat) ^ 2 / sigma2hat
    out[i] <- max(tmp)
  }
  
  out
  
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
  list("Omni" = tmpOmni, "Ind" = tmpInd) 
}