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


invert.q <- function(coef) {
  out <- 1
  
  if (all(abs(coef) < 1)) {
    minmod <- min(Mod(polyroot(c(1, coef))))
    
    if (minmod <= 1) {
      out <- 0
    }
  } else {
    out <- 0
  }
  
  return(out)
}

pars.mat <- function(n, parsVec, norder = 1) {
  Check <- invert.q(parsVec)
  if (Check == 0) {
    NULL
  } else {
    Mat <- diag(n)
    for (i in 1:norder) {
      Mat <- Mat + Diag(rep(parsVec[i], n - i), k = -i)
    }
    Mat
  }
}


sigma.mat <- function(n, order = c(1, 0, 0), phi.vec = 0.5, theta.vec = NULL, sigma2 = 1, burn.in = 50) {
  if (order[1] == 0) {
    phiMat <- diag(n + burn.in)
  } else {
    phiMat <- pars.mat(n + burn.in, -phi.vec, norder = order[1])
  }
  
  if (order[3] == 0) {
    thetaMat <- diag(n + burn.in)
  } else {
    thetaMat <- pars.mat(n + burn.in, theta.vec, norder = order[3])
  }
  
  out <- solve(phiMat) %*% thetaMat %*% t(thetaMat) %*% t(solve(phiMat)) * sigma2
  
  gamma0 <- out[dim(out)[1], dim(out)[2]]
  
  if (burn.in > 0) {
    out <- out[-c(1:burn.in), -c(1:burn.in)]
  }
  
  list(sigma.mat = out, sqrtsigma.mat = sqrtm(out)$B, gamma0 = gamma0)
}

#' simulate realizations using ARMA(p, q) and one sustained shift
#' 
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
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
rarma <- function(object, n, h, delta, xreg = NULL, nsim = 100, burnin = 50, lowerbound = 0) {
  
  nphi <- grep("ar", names(object$coef))
  ntheta <- grep("ma", names(object$coef))
  
  if (length(nphi) == 0) {
    phi <- NULL
  } else {
    phi <- object$coef[nphi]
  }
  
  if (length(ntheta) == 0) {
    theta <- NULL
  } else {
    theta <- object$coef[ntheta]
  }
  
  gamma0 <- sigma.mat(n = nsim, order = c(length(nphi), 0, length(ntheta)), 
                      phi.vec = phi, theta.vec = theta, sigma2 = object$sigma2, burn.in = burnin)$gamma0
  
  tmpint <- grep("intercept", names(object$coef))
  
  mu <- rep(ifelse(length(tmpint) == 0, 0, object$coef[tmpint]), n)
  mu[h:n] <- mu[h:n] + sqrt(gamma0) * delta
  
  innov <- rnorm(n, mu, sqrt(object$sigma2))

  ts <- simulate(object, nsim = n, future = FALSE, innov = innov, xreg = xreg)

  ts[which(ts < lowerbound)] <- lowerbound
  ts
}

#' Caculate the moving averages
#' 
#' gets the moving averages
#' @param Y is the input
#' @param w is the length of moving window
#' 
#' @export
#' @examples
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = 100)
#' ma <- movaver(Y, w)
movaver <- function(Y, w = 5){filter(Y, rep(1 / w, w), sides = 1)}

#' Transform the data
#' 
#' gets the transformed input
#' @param Y is the input
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization
#' 
#' @export
#' @examples
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = 100)
#' ma <- movaver(Y, w)
#' ma.tr <- trans(ma, TRUE, 1, TRUE)
trans <- function(Y, log = TRUE, const = 1, sta = TRUE, meanY = NULL, sdY = NULL){
  out <- Y
  if (log == TRUE) {
    out <- log(out + const)
  }
  if (sta == TRUE) {
    meanY <- ifelse(!is.null(meanY), meanY, mean(out))
    sdY <- ifelse(!is.null(sdY), sdY, sd(out))
    out <- (out - meanY) / sdY
  }
  
  list(
    "Y" = out,
    "meanY" = meanY,
    "sdY" = sdY
  )
}

#' Back-transform the data
#' 
#' gets the back-transformed input
#' @param Y is the input
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization
#' @param meanY is the mean of the original Y
#' @param sdY is the standard deviation of the original Y
#' 
#' @export
#' @examples
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = 100)
#' ma <- movaver(Y, w)
#' ma.tr <- trans(ma, TRUE, 1, TRUE)
#' ma.batr <- backtrans(ma.tr$Y, TRUE, 1, TRUE, ma.tr$meanY, ma.tr$sdY)
backtrans <- function(Y, log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1){
  out <- Y
  
  if (sta == TRUE) {
    out <- out * sdY + meanY
  }
  
  if (log == TRUE) {
    out <- exp(out) - const
  }
  
  out
}

#' simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
#' 
#' @param Ph1BayesianLASSO.model is the length
#' @param log is the log
#' @param const is the constant
#' @param sta is the sta
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
RMSE.ph1 <- function(Ph1BayesianLASSO.model, log = TRUE, const = 1, sta = TRUE, lowerbound = 0) {
  nsim <- dim(Ph1BayesianLASSO.model$Phi)[2]
  q <- dim(Ph1BayesianLASSO.model$Phi)[1]
  TT <- length(Ph1BayesianLASSO.model$Y.tr)
  
  tmpfit.tr <- rep(NA, TT)
  tmpfit.ma <- rep(NA, TT)
  tmpresi.tr <- rep(NA, TT)
  tmpresi.ma <- rep(NA, TT)
  RMSE.tr <- rep(NA, nsim)
  RMSE.ma <- rep(NA, nsim)
  
  Y.tr <- Ph1BayesianLASSO.model$Y.tr
  Y.ma <- Ph1BayesianLASSO.model$Y.ma
  meanY <- Ph1BayesianLASSO.model$meanY
  sdY <- Ph1BayesianLASSO.model$sdY
  
  X <- Ph1BayesianLASSO.model$X
  H <- Ph1BayesianLASSO.model$H
  
  for (j in 1:nsim) {
    tmpmuq <- Ph1BayesianLASSO.model$muq[j]
    tmpPhi <- Ph1BayesianLASSO.model$Phi[, j]
    tmpV <- Y.tr - tmpmuq
    if (!is.null(X)) {
      tmpBeta <- Ph1BayesianLASSO.model$Beta[, j]
      tmpKappa <- Ph1BayesianLASSO.model$Kappa[, j]
      tmpV <- tmpV - X %*% (tmpBeta * tmpKappa) 
    }
    if (!is.null(H)) {
      tmpGamma <- Ph1BayesianLASSO.model$Gamma[, j]
      tmpTau <- Ph1BayesianLASSO.model$Tau[, j]
      tmpV <- tmpV - H %*% (tmpGamma * tmpTau)
    }
    for (i in (q + 1):TT) {
      tmpfit.tr[i] <- tmpmuq + tmpV[(i - 1):(i - q)] %*% tmpPhi
      
      if (!is.null(X)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + X[i, ] %*% (tmpBeta * tmpKappa)
      }
      if (!is.null(H)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + H[i, ] %*% (tmpGamma * tmpTau)
      }
      
      tmpresi.tr[i] <- Y.tr[i] - tmpfit.tr[i]
      tmpfit.ma[i] <- backtrans(tmpfit.tr[i], log, const, sta, meanY, sdY)
      tmpfit.ma[i] <- ifelse(tmpfit.ma[i] < lowerbound, lowerbound, tmpfit.ma[i])
      tmpresi.ma[i] <- Y.ma[i] - tmpfit.ma[i]
    }
    tmpresi.tr <- tmpresi.tr[(q + 1):TT]
    RMSE.tr[j] <- sqrt(t(tmpresi.tr) %*% tmpresi.tr / (TT - q))
    
    tmpresi.ma <- tmpresi.ma[(q + 1):TT]
    RMSE.ma[j] <- sqrt(t(tmpresi.ma) %*% tmpresi.ma / (TT - q))
  }
  
  list("RMSE.tr" = RMSE.tr, "RMSE.ma" = RMSE.ma)
  
}

#' simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
#' 
#' @param Ph1BayesianLASSO.model is the length
#' @param log is the log
#' @param const is the constant
#' @param sta is the sta
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
RMSE.ph2 <- function(Y, Ph1BayesianLASSO.model, X = NULL, H = NULL, 
                     log = TRUE, const = 1, sta = TRUE, lowerbound = 0) {
  
  TT <- length(Y)
  YY <- c(Ph1BayesianLASSO.model$Y, Y)
  nn <- length(YY)
  q <- dim(Ph1BayesianLASSO.model$Phi)[1]
  
  Y.ma <- movaver(YY, Ph1BayesianLASSO.model$w)[(nn - TT - q + 1):nn]
  tt <- length(Y.ma)
  
  meanY <- Ph1BayesianLASSO.model$meanY
  sdY <- Ph1BayesianLASSO.model$sdY
  
  Y.tr <- trans(Y.ma, log, const, sta, meanY, sdY)$Y
  
  nsim <- dim(Ph1BayesianLASSO.model$Phi)[2]

  
  tmpfit.tr <- rep(NA, tt)
  tmpfit.ma <- rep(NA, tt)
  tmpresi.tr <- rep(NA, tt)
  tmpresi.ma <- rep(NA, tt)
  RMSE.tr <- rep(NA, nsim)
  RMSE.ma <- rep(NA, nsim)
  
  if (!is.null(X)) {
    X <- as.matrix(rbind(Ph1BayesianLASSO.model$X, X)[(nn - TT - q + 1):nn, ])
  }
  if (!is.null(H)) {
    H <- as.matrix(rbind(Ph1BayesianLASSO.model$H, H)[(nn - TT - q + 1):nn, ])
  }
  
  for (j in 1:nsim) {
    tmpmuq <- Ph1BayesianLASSO.model$muq[j]
    tmpPhi <- Ph1BayesianLASSO.model$Phi[, j]
    tmpV <- Y.tr - tmpmuq
    if (!is.null(X)) {
      tmpBeta <- Ph1BayesianLASSO.model$Beta[, j]
      tmpKappa <- Ph1BayesianLASSO.model$Kappa[, j]
      tmpV <- tmpV - X %*% (tmpBeta * tmpKappa) 
    }
    if (!is.null(H)) {
      tmpGamma <- Ph1BayesianLASSO.model$Gamma[, j]
      tmpTau <- Ph1BayesianLASSO.model$Tau[, j]
      tmpV <- tmpV - H %*% (tmpGamma * tmpTau)
    }
    
    for (i in (q + 1):tt) {
      tmpfit.tr[i] <- tmpmuq + tmpV[(i - 1):(i - q)] %*% tmpPhi
      
      if (!is.null(X)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + X[i, ] %*% (tmpBeta * tmpKappa)
      }
      if (!is.null(H)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + H[i, ] %*% (tmpGamma * tmpTau)
      }
      
      tmpresi.tr[i] <- Y.tr[i] - tmpfit.tr[i]
      tmpfit.ma[i] <- backtrans(tmpfit.tr[i], log, const, sta, meanY, sdY)
      tmpfit.ma[i] <- ifelse(tmpfit.ma[i] < lowerbound, lowerbound, tmpfit.ma[i])
      tmpresi.ma[i] <- Y.ma[i] - tmpfit.ma[i]
    }
    tmpresi.tr <- tmpresi.tr[(q + 1):tt]
    RMSE.tr[j] <- sqrt(t(tmpresi.tr) %*% tmpresi.tr / (tt - q))
    
    tmpresi.ma <- tmpresi.ma[(q + 1):tt]
    RMSE.ma[j] <- sqrt(t(tmpresi.ma) %*% tmpresi.ma / (tt - q))
  }
  
  list("RMSE.tr" = RMSE.tr, "RMSE.ma" = RMSE.ma)
  
}
