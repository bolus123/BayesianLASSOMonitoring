#' Get a simulation transformed data using Random Flexible Level Shift Model in the retrospective phase under H0 without censoring
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param n is the length of simulated process
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param X is the matrix of X
#' @param Beta is the matrix of coefficients
#' @param Kappa is the matrix of triggering the corresponding beta
#' @param Y1 is the transformed vector
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
GibbsRFLSM.sim.ph1.H0 <- function(nsim, n, Phi, muq, sigma2, 
                                  X = NULL, Beta = NULL, Kappa = NULL,
                                  Y1 = rep(median(muq), n)) {
  
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  Y.sim <- matrix(Y1, nrow = n, ncol = nsim)
  muX <- matrix(0, nrow = n, ncol = m)
  
  if (!is.null(X)) {
    muX <- X %*% (Beta * Kappa)
  }
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpmu0 <- muq[tmpsel] + muX[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    for (i in (q + 1):n) {
      Y.sim[i, j] <- tmpmu0[i] + rnorm(1, 
                 mean = (Y.sim[(i - 1):(i - q), j] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi, 
                 sd = sqrt(tmpsigma2))
    }
    
  }
  
  Y.sim[-c(1:q), ]
  
}

lower.bound.ma <- function(lower.ma = 0, log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1) {
  
  out <- lower.ma
  
  if (log == TRUE) {
    out <- log(out + const)
  }
  
  if (sta == TRUE) {
    out <- (out - meanY) / sdY
  }
  
  out
  
}

#' Get a simulated moving averages using Random Flexible Level Shift Model in the retrospective phase under H0 
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param n is the length of simulated process
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param X is the matrix of X
#' @param Beta is the matrix of coefficients
#' @param Kappa is the matrix of triggering the corresponding beta
#' @param Y1 is the transformed vector
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization
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
GibbsRFLSM.sim.ma.ph1.H0 <- function(nsim, n, Phi, muq, sigma2, 
                                  X = NULL, Beta = NULL, Kappa = NULL, 
                                  Y1 = rep(median(muq), n),
                                  log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1) {
  
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  Y.sim <- matrix(Y1, nrow = n, ncol = nsim)
  muX <- matrix(0, nrow = n, ncol = m)
  
  lower.bound <- lower.bound.ma(lower.ma = 0, log = log, const = const, sta = sta, 
                                meanY = meanY, sdY = sdY)
  
  if (!is.null(X)) {
    muX <- X %*% (Beta * Kappa)
  }
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpmu0 <- muq[tmpsel] + muX[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    for (i in (q + 1):n) {
      Y.sim[i, j] <- tmpmu0[i] + rnorm(1, 
                                       mean = (Y.sim[(i - 1):(i - q), j] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi, 
                                       sd = sqrt(tmpsigma2))
      Y.sim[i, j] <- ifelse(Y.sim[i, j] < lower.bound, lower.bound, Y.sim[i, j]) 
    }
    
  }
  
  out <- backtrans(Y.sim, log = log, const = const, sta = sta, 
            meanY = meanY, sdY = sdY)
    
  out[-c(1:q), ]
}

lower.bound.count <- function(ma, Y0, w, log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1) {
  
  out <- ma - Y0 / w
  
  if (log == TRUE) {
    out <- log(out + const)
  }
  
  if (sta == TRUE) {
    out <- (out - meanY) / sdY
  }
  
  out
  
}

#' Get a simulated count using Random Flexible Level Shift Model in the retrospective phase under H0 
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param n is the length of simulated process
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param X is the matrix of X
#' @param Beta is the matrix of coefficients
#' @param Kappa is the matrix of triggering the corresponding beta
#' @param Y1 is the transformed vector
#' @param Y0 is the raw count
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization
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
GibbsRFLSM.sim.count.ph1.H0 <- function(nsim, n, w, Phi, muq, sigma2,
                                     X = NULL, Beta = NULL, Kappa = NULL, 
                                     Y1 = rep(median(muq), n),
                                     Y0 = rep(0, n + w - 1),
                                     log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1) {
  
  n.Y0 <- length(Y0)
  Y0 <- Y0[(n.Y0 - (n + w - 1) + 1):n.Y0]
  ma <- movaver(Y0)[w:(n + w - 1)]
  
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  Y.sim <- matrix(Y1, nrow = n, ncol = nsim)
  muX <- matrix(0, nrow = n, ncol = m)
  
  if (!is.null(X)) {
    muX <- X %*% (Beta * Kappa)
  }
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpmu0 <- muq[tmpsel] + muX[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    for (i in (q + 1):n) {
      Y.sim[i, j] <- tmpmu0[i] + rnorm(1, 
                                       mean = (Y.sim[(i - 1):(i - q), j] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi, 
                                       sd = sqrt(tmpsigma2))
      lower.bound <- lower.bound.count(ma[i - 1], Y0[i - 1], w, log = log, const = const, sta = sta, 
                                       meanY = meanY, sdY = sdY)
      Y.sim[i, j] <- ifelse(Y.sim[i, j] < lower.bound, lower.bound, Y.sim[i, j]) 
    }
    
  }
  
  out <- backtrans(Y.sim, log = log, const = const, sta = sta, 
            meanY = meanY, sdY = sdY)
  
  for (i in (q + 1):n) {
    out[i, ] <- round(w * out[i, ] - (w * ma[i - 1] - Y0[i - 1]))
  }
  
  out[-c(1:q), ]
  
}