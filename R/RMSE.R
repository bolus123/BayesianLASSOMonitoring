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
#' RMSE(Y, result$Phi, result$Mu)
RMSE <- function(Y, Phi, Mu) {
  
  T <- length(Y)
  q <- dim(Phi)[1]
  nsim <- dim(Phi)[2]
  ee <- matrix(NA, nrow = T - q, ncol = nsim)
  rmse <- rep(0, nsim)
  
  for (ii in seq(nsim)) {
    V <- matrix(Y, ncol = 1) - Mu[, ii]
    Vas <- getV(V, q)
    V <- V[-c(1:q)]
    Vas <- Vas[-c(1:q), ]
    ee[, ii] <- V - Vas %*% Phi[, ii]
    rmse[ii] <- sqrt(sum(ee[, ii] ^ 2) / (T - q))
  }
  rmse
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
#' Residual(Y, result$Phi, result$Mu)
Residual <- function(Y, Phi, Mu) {
  
  T <- length(Y)
  q <- dim(Phi)[1]
  nsim <- dim(Phi)[2]
  ee <- matrix(NA, nrow = T - q, ncol = nsim)
  
  for (ii in seq(nsim)) {
    V <- matrix(Y, ncol = 1) - Mu[, ii]
    Vas <- getV(V, q)
    V <- V[-c(1:q)]
    Vas <- Vas[-c(1:q), ]
    ee[, ii] <- V - Vas %*% Phi[, ii]
  }
  ee
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
Fit <- function(Y, Phi, Mu) {
  
  T <- length(Y)
  q <- dim(Phi)[1]
  nsim <- dim(Phi)[2]
  ff <- matrix(NA, nrow = T - q, ncol = nsim)
  
  for (ii in seq(nsim)) {
    V <- matrix(Y, ncol = 1) - Mu[, ii]
    Vas <- getV(V, q)
    V <- V[-c(1:q)]
    Vas <- Vas[-c(1:q), ]
    ff[, ii] <- Vas %*% Phi[, ii] + Mu[-c(1:q), ii]
  }
  rbind(Mu[1:q, ], ff)
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
GibbsRFLSM.simmax.Yao <- function(Y, Phi, Mu, sigma2, 
                              nsim = 1000) {
  
  q <- dim(Phi)[1]
  out <- rep(NA, nsim)
  xbar <- mean(Y[-c(1:q)])
  std <- sd(Y[-c(1:q)])
  m <- dim(Phi)[2]
  
  for (i in seq(nsim)) {
    
    k <- sample(1:m, 1)
    tmpPhi <- Phi[, k]
    tmpMu <- Mu[, k]
    tmpsigma2 <- sigma2[k]
    
    tmp <- c(Y, tmpPhi, tmpMu, tmpsigma2)
    out[i] <- max(((tmp - xbar) / std) ^ 2)
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
GibbsRFLSM.simmax.residual <- function(Y, Phi, Mu, sigma2, 
                                  Phihat, Muhat, sigma2hat, 
                                  nsim = 1000) {
  
  q <- dim(Phi)[1]
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
    out[i] <- max((tmpV - tmpVas %*% Phihat) ^ 2 / sigma2hat)
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
#' RSS(Y, result$Phi, result$Mu)
RSS <- function(Y, Phi, Mu) {
  
  T <- length(Y)
  q <- dim(Phi)[1]
  nsim <- dim(Phi)[2]
  ee <- matrix(NA, nrow = T - q, ncol = nsim)
  rss <- rep(0, nsim)
  
  for (ii in seq(nsim)) {
    V <- matrix(Y, ncol = 1) - Mu[, ii]
    Vas <- getV(V, q)
    V <- V[-c(1:q)]
    Vas <- Vas[-c(1:q), ]
    ee[, ii] <- V - Vas %*% Phi[, ii]
    rss[ii] <- sum(ee[, ii] ^ 2)
  }
  rss
}