
#' obtain the root squared error
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @export
#' @examples
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, q, diag(nrow = q), 0.1, 0.1, 0.1, 0.1, 
#' 1, 1, 0.1, "MonoALASSO", Inf, 0, 1000, 1, 100, 1e-10, H)
#'
#' Fit(Y, result$Phi, result$Mu)
GibbsRFLSM.sim <- function(Y, Phi, Mu, sigma2) {
  
  TT <- length(Y)
  q <- length(Phi)
  
  sim <- rep(NA, TT)
  
  sim[1:q] <- Y[1:q]
  
  for (ii in (q + 1):TT) {
    sim[ii] <- Mu[ii] + (sim[(ii - 1):(ii - q)] - Mu[(ii - 1):(ii - q)]) %*% 
      Phi + rnorm(1, 0, sqrt(sigma2))
  }
  
  sim 
}

#' obtain the root squared error
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @export
#' @examples
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, q, diag(nrow = q), 0.1, 0.1, 0.1, 0.1, 
#' 1, 1, 0.1, "MonoALASSO", Inf, 0, 1000, 1, 100, 1e-10, H)
#'
#' Fit(Y, result$Phi, result$Mu)
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

#' obtain the root squared error
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @export
#' @examples
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, q, diag(nrow = q), 0.1, 0.1, 0.1, 0.1, 
#' 1, 1, 0.1, "MonoALASSO", Inf, 0, 1000, 1, 100, 1e-10, H)
#'
#' Fit(Y, result$Phi, result$Mu)
GibbsRFLSM.PPP.residual <- function(Y, Phi, muq, sigma2, 
                                    Phihat, muqhat, sigma2hat, 
                                    nsim = 1000) {
  
  q <- dim(Phi)[1]
  n <- length(Y)
  m <- dim(Phi)[2]
  
  Muq <- matrix(muq, nrow = n, ncol = m, byrow = T)
  Muqhat <- rep(muqhat, n)
  
  ccrep <- GibbsRFLSM.simmax.residual(Y, Phi, Muq, sigma2, 
                                      Phihat, Muqhat, sigma2hat, 
                                      nsim)
  
  tmp <- Y - muqhat
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