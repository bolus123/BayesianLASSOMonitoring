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


fdelta <- function(alpha, lambda0, pi0, lambda1, pi1) {
  
  a <- sqrt(1 - sum(alpha ^ 2)) / (1 - sum(alpha))
  d1 <- (1 - pi1) * lambda1 / 
    sqrt(lambda1 * (1 - pi1) * (1 + pi1 * lambda1))
  d0 <- (1 - pi0) * lambda0 / 
    sqrt(lambda0 * (1 - pi0) * (1 + pi0 * lambda0))
  
  a * (d1 - d0)
  
}

flambda1 <- function(delta, alpha, lambda0, pi0, pi1 = pi0, interval = c(1e-10, 1000)) {
  
  rootfinding <- function(lambda1, delta, alpha, lambda0, pi0, pi1) {
    tmp <- fdelta(alpha, lambda0, pi0, lambda1, pi1) 
    #cat("tmp:", tmp, "lambda1:", lambda1, "\n")
    delta - tmp
  }
  
  uniroot(rootfinding, interval = interval, delta = delta, 
          alpha = alpha, 
          lambda0 = lambda0, pi0 = pi0, pi1 = pi1)$root
  
}



fpi1 <- function(delta, alpha, lambda0, pi0, lambda1 = lambda0, 
                 interval = c(1e-10, 1 - 1e-10)) {
  
  rootfinding <- function(pi1, delta, alpha, lambda0, pi0, lambda1) {
    tmp <- fdelta(alpha, lambda0, pi0, lambda1, pi1) 
    #cat("tmp:", tmp, "pi1:", pi1, "\n")
    delta - tmp
  }
  
  uniroot(rootfinding, interval = interval, delta = delta, 
          alpha = alpha, 
          lambda0 = lambda0, lambda1 = lambda1, pi0 = pi0)$root
  
}


#' simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
#' 
#' @param n is the length
#' @param alpha is the alpha
#' @param lambda is the mean of poisson mixture
#' @param pi is the proportion of zeros
#' @param h is the start point of shift
#' @param delta is the value of the standardized shift
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rzinpoisinar3 <- function(n, alpha, lambda, pi, h, delta, burnin = 100) {
  
  q <- length(alpha)
  out <- rep(NA, n + burnin + q)
  out[1:q] <- VGAM::rzipois(q, lambda, pi)
  
  
  k <- 0
  
  lambda1 <- flambda1(delta, alpha, lambda, pi, pi1 = pi, interval = c(1e-10, 1000))
  
  for (i in (q + 1):(n + burnin + q)) {
    for (j in 1:q) {
      out[i] <- rbinom(1, out[i - j], alpha[j])
    }
    
    if (i >= (q + 1 + burnin)) {
      k <- k + 1
    }
    
    if (k >= h) {
      out[i] <- out[i] + VGAM::rzipois(1, lambda1, pi)
    } else {
      out[i] <- out[i] + VGAM::rzipois(1, lambda, pi)
    }
    
    
  }
  
  out[(burnin + q + 1):(n + burnin + q)]
  
}


