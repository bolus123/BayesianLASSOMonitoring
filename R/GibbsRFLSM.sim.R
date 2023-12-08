#' Get a simulation transformed data using Random Flexible Level Shift Model in the retrospective phase under H0
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
#' @param bt is the flag of triggering the back transformation
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization after the log transformation
#' @param meanY is the mean of the log Y
#' @param sdY is the standard deviation of the log Y
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
                                  Y1 = rep(median(muq), n), 
                                  bt = TRUE, log = TRUE, const = 1, sta = TRUE, 
                                  meanY = 0, sdY = 1, 
                                  lower.bound = -Inf, rounding = FALSE) {
  
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
      Y.sim[i, j] <- tmpmu0[i] + (Y1[(i - 1):(i - q)] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi + 
                  rnorm(1, mean = 0, sd = sqrt(tmpsigma2))
    }
    
  }
  
  if (bt == TRUE) {
    Y.sim <- backtrans(Y.sim, log = log, const = const, sta = sta, 
                       meanY = meanY, sdY = sdY)
    
    
  }
  
  Y.sim[Y.sim < lower.bound] <- lower.bound
  
  if (rounding == TRUE) {
    Y.sim <- round(Y.sim)
  }
  
  Y.sim[-c(1:q), ]
  
}

#' Get a simulation transformed data using Random Flexible Level Shift Model in the prospective phase under H0
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param h is the steps of forward predictions
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param X is the matrix of X
#' @param Beta is the matrix of coefficients
#' @param Kappa is the matrix of triggering the corresponding beta
#' @param Y1 is the transformed vector
#' @param bt is the flag of triggering the back transformation
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization after the log transformation
#' @param meanY is the mean of the log Y
#' @param sdY is the standard deviation of the log Y
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
GibbsRFLSM.sim.ph2.H0 <- function(nsim, h, Phi, muq, sigma2, 
                                  X = NULL, Beta = NULL, Kappa = NULL,
                                  Y1 = rep(median(muq), dim(Phi)[1]), 
                                  bt = TRUE, log = TRUE, const = 1, sta = TRUE, 
                                  meanY = 0, sdY = 1) {
  
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  n.Y1 <- length(Y1)
  Y1 <- Y1[(n.Y1 - q + 1):n.Y1]
  
  Y1.sim <- matrix(Y1, nrow = q, ncol = nsim)
  Y2.sim <- matrix(NA, nrow = h, ncol = nsim)
  Y.sim <- rbind(Y1.sim, Y2.sim)
  
  muX <- matrix(0, nrow = q + h, ncol = m)
  
  if (!is.null(X)) {
    n.X <- dim(X)[1]
    X <- X[(n.X - q - h + 1):n.X, ]
    muX <- X %*% (Beta * Kappa)
  }
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpmu0 <- muq[tmpsel] + muX[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    for (i in (q + 1):(h + q)) {
      Y.sim[i, j] <- tmpmu0[i] + (Y.sim[(i - 1):(i - q)] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi + 
        rnorm(1, mean = 0, sd = sqrt(tmpsigma2))
    }
    
  }
  
  if (bt == TRUE) {
    Y.sim <- backtrans(Y.sim, log = log, const = const, sta = sta, 
                       meanY = meanY, sdY = sdY)
    
    
  }
  
  
  Y.sim[-c(1:q), ]
  
}

#' Get a simulation transformed data using Random Flexible Level Shift Model in the retrospective phase under H1
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param n is the length of simulated process
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param Y1 is the transformed vector
#' @param bt is the flag of triggering the back transformation
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization after the log transformation
#' @param meanY is the mean of the log Y
#' @param sdY is the standard deviation of the log Y
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
GibbsRFLSM.sim.ph1.H1 <- function(nsim, n, Phi, Mu, sigma2, 
                                  Y1 = apply(model$Mu, 1, median), 
                                  bt = TRUE, log = TRUE, const = 1, sta = TRUE, 
                                  meanY = 0, sdY = 1) {
  
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  Y.sim <- matrix(Y1, nrow = n, ncol = nsim)
 
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpMu <- Mu[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    for (i in (q + 1):n) {
      Y.sim[i, j] <- tmpMu[i] + (Y1[(i - 1):(i - q)] - tmpMu[(i - 1):(i - q)]) %*% tmpPhi + 
        rnorm(1, mean = 0, sd = sqrt(tmpsigma2))
    }
    
  }
  
  if (bt == TRUE) {
    Y.sim <- backtrans(Y.sim, log = log, const = const, sta = sta, 
                       meanY = meanY, sdY = sdY)
    
    
  }
  
  
  Y.sim[-c(1:q), ]
  
}


#' Get a simulation transformed data using Random Flexible Level Shift Model in the prospective phase under H1
#' 
#' gets the simulated data
#' @param nsim is the number of simulations
#' @param h is the steps of forward predictions
#' @param Phi is the matrix of laggy coefficients
#' @param muq is the vector of grand mean
#' @param sigma2 is the vector of variance
#' @param H is the matrix of H
#' @param Gamma is the matrix of shift coefficients
#' @param Tau is the matrix of triggering the corresponding gamma
#' @param X is the matrix of X
#' @param Beta is the matrix of coefficients
#' @param Kappa is the matrix of triggering the corresponding beta
#' @param Y1 is the transformed vector
#' @param bt is the flag of triggering the back transformation
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization after the log transformation
#' @param meanY is the mean of the log Y
#' @param sdY is the standard deviation of the log Y
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
GibbsRFLSM.sim.ph2.H1 <- function(nsim, h, Phi, muq, sigma2, 
                                  H = NULL, Gamma = NULL, Tau = NULL,
                                  X = NULL, Beta = NULL, Kappa = NULL,
                                  Y1 = rep(median(muq), dim(Phi)[1]), 
                                  bt = TRUE, log = TRUE, const = 1, sta = TRUE, 
                                  meanY = 0, sdY = 1) {
  
  m <- dim(Phi)[2]
  q <- dim(Phi)[1]
  n.Y1 <- length(Y1)
  Y1 <- Y1[(n.Y1 - q + 1):n.Y1]
  
  Y1.sim <- matrix(Y1, nrow = q, ncol = nsim)
  Y2.sim <- matrix(NA, nrow = h, ncol = nsim)
  Y.sim <- rbind(Y1.sim, Y2.sim)
  
  muH <- matrix(0, nrow = q + h, ncol = m)
  
  if (!is.null(H)) {
    n.H <- dim(H)[1]
    H <- H[(n.H - q - h + 1):n.H, ]
    muH <- H %*% (Gamma * Tau)
  }
  
  muX <- matrix(0, nrow = q + h, ncol = m)
  
  if (!is.null(X)) {
    n.X <- dim(X)[1]
    X <- X[(n.X - q - h + 1):n.X, ]
    muX <- X %*% (Beta * Kappa)
  }
  
  for (j in 1:nsim) {
    tmpsel <- sample(1:m, 1)
    tmpPhi <- Phi[, tmpsel]
    tmpmu0 <- muq[tmpsel] + muX[, tmpsel] + muH[, tmpsel]
    tmpsigma2 <- sigma2[tmpsel]
    
    for (i in (q + 1):(h + q)) {
      Y.sim[i, j] <- tmpmu0[i] + (Y.sim[(i - 1):(i - q)] - tmpmu0[(i - 1):(i - q)]) %*% tmpPhi + 
        rnorm(1, mean = 0, sd = sqrt(tmpsigma2))
    }
    
  }
  
  if (bt == TRUE) {
    Y.sim <- backtrans(Y.sim, log = log, const = const, sta = sta, 
                       meanY = meanY, sdY = sdY)
    
    
  }
  
  
  Y.sim[-c(1:q), ]
  
}