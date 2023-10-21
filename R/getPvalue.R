#' get the P value for RFLSM with kernel smoothing
#' 
#' @param TauGamma is the distributions of shifts.
#' @param tail is type of the test.
getPvalueKSRFLSM <- function(TauGamma, tail = "2-sided") {
  tmp <- TauGamma
  
  nn <- dim(tmp)[1]
  pvalue <- rep(0, nn)
  
  for (i in 1:nn) {
    
    tmpkde <- density(tmp[i, ])
    rtmpdens <- spatstat.explore::CDF(tmpkde)
    tmpp <- rtmpdens(0)
    
    if (tail == "2-sided") {
      pvalue[i] <- 2 * min(1 - tmpp, tmpp)
    } else if (tail == "left-sided") {
      pvalue[i] <- tmpp
    } else if (tail == "right-sided") {
      pvalue[i] <- 1 - tmpp
    }
    
  }
  return(pvalue)
}


#' get the P value for RFLSM
#' 
#' @param TauGamma is the distributions of shifts.
#' @param tail is type of the test.
getPvalueRFLSM <- function(TauGamma, tail = "2-sided") {
  tmp <- TauGamma
  
  nn <- dim(tmp)[1]
  pvalue <- rep(0, nn)
  
  for (i in 1:nn) {
    
    if (tail == "2-sided") {
      pvalue[i] <- 2 * min(mean(tmp[i, ] < 0), mean(tmp[i, ] > 0))
    } else if (tail == "left-sided") {
      pvalue[i] <- mean(tmp[i, ] < 0)
    } else if (tail == "right-sided") {
      pvalue[i] <- mean(tmp[i, ] > 0)
    }
    
  }
  return(pvalue)
}


#' get a vector of p values
#' 
#' @param TauGamma is the distributions of shifts.
#' @param tail is type of the test.
#' @param method get p values with or without kernel smoothing.
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
#' getPvalue(result$Tau * result$Gamma)
getPvalue <- function(TauGamma, tail = "2-sided", method = "raw") {
  
  if (method == "raw") {
    pvalue <- getPvalueRFLSM(TauGamma, tail)
  } else if (method == "ks") {
    pvalue <- getPvalueKSRFLSM(TauGamma, tail) 
  }
  
  pvalue
  
}


#' caculate posterior predictve pvalue
#' 
#' @param Y is Y
#' @param Phi is the lagged coefficients.
#' @param Mu is the individual means
#' @param TauGamma is the distributions of shifts.
#' @param sigma2 is the variance of errors.
#' @param muq is the intercept.
#' @param method is the method to estimate the parameters.
#' @param nsim is the number of simulation for the simulation
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
#' PPP(Y0, result$Phi, result$Mu, H, result$Tau * result$Gamma, result$sigma2, result$muq)
#' 
PPP <- function(Y, Phi, Mu, H, TauGamma, sigma2, muq, method = "median", nsim = 3000) {
  
  nnsim <- dim(Phi)[2]
  q <- dim(Phi)[1]
  T <- length(Y)
  m <- T - q
  
  YSim <- matrix(NA, nrow = T, ncol = nsim)
  RSS0 <- rep(NA, nsim)
  RSS1 <- rep(NA, nsim)
  RSS0Sim <- rep(NA, nsim)
  RSS1Sim <- rep(NA, nsim)
  out <- rep(NA, nsim)

  tmpPhi <- matrix(NA, nrow = q)
  tmpGamma <- matrix(NA, nrow = m)
  tmpMu <- matrix(NA, nrow = T)
  
  if (method == "median") {
    tmpmuq <- median(muq)
    tmpsigma2 <- median(sigma2)
    for (i in 1:q) {
      tmpPhi[i, 1] <- median(Phi[i, ]) 
    }
    
    for (i in 1:m) {
      tmpGamma[i, 1] <- median(TauGamma[i, ])
    }
  }
  
  tmpMu <- tmpmuq + H %*% tmpGamma
  
  RSS0 <- RSS(Y, matrix(tmpPhi, ncol = 1), 
              matrix(tmpmuq, ncol = 1))
  RSS1 <- RSS(Y, matrix(tmpPhi, ncol = 1), 
              matrix(tmpMu, ncol = 1))
  
  
  for (i in 1:nsim) {
    
    tmpsel <- sample(1:nnsim, 1)
    YSim[1:q, i] <- Y[1:q]
    
    for (j in (q + 1):T) {
      YSim[j, i] <- Mu[j, tmpsel] + 
        (YSim[(j - 1):(j - q), i] - Mu[(j - 1):(j - q), tmpsel]) %*% Phi[, tmpsel] + 
        rnorm(1, sd = sqrt(sigma2[tmpsel]))
    }
    
    
    RSS0Sim[i] <- RSS(YSim[, i], matrix(tmpPhi, ncol = 1), 
                      matrix(tmpmuq, ncol = 1))
    RSS1Sim[i] <- RSS(YSim[, i], matrix(tmpPhi, ncol = 1), 
                      matrix(tmpMu, ncol = 1))
    
  }
  
  out <- ((RSS0Sim - RSS1Sim) / tmpsigma2) >= 
    ((RSS0 - RSS1) / tmpsigma2)
  mean(out)
  
}