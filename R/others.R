#' obtain the root mean squared error
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
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
#' RMSE(Y, result$Phi, result$Mu)
#' 
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

#' obtain the residuals
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
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
#' Residual(Y, result$Phi, result$Mu)
#' 
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

#' obtain the fits
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
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
#' Fit(Y, result$Phi, result$Mu)
#' 
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




#' obtain resiudal sum squares
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
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
#' RSS(Y, result$Phi, result$Mu)
#' 
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

#' obtain the log likelihood
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @param sigma2 is the variance
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
#' ll(Y, result$Phi, result$Mu, result$sigma2)
#' 
ll <- function(Y, Phi, Mu, sigma2) {
  
  T <- length(Y)
  q <- dim(Phi)[1]
  nsim <- length(sigma2)
  ll <- matrix(NA, nrow = T - q, ncol = nsim)
  
  for (ii in seq(nsim)) {
    V <- matrix(Y, ncol = 1) - Mu[, ii]
    Vas <- getV(V, q)
    V <- V[-c(1:q)]
    Vas <- Vas[-c(1:q), ]
    resi <- V - Vas %*% Phi[, ii]
    ll[, ii] <- dnorm(resi, mean = 0, sd = sqrt(sigma2[ii]), log = TRUE)
  }
  
  ll
  
}


