#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param Y is a vector.
#' @param H is the design matrix for shifts.
#' @param X is the input matrix
#' @param Y0 is the initial Y
#' @param q is the number of lags.
#' @param A is a given variance-covariance matrix in MT and regression for the slab-and-spike coefficients.
#' @param a is a given shape of the prior gamma distribution for sigma2.
#' @param b is a given scale of the prior gamma distribution for sigma2.
#' @param alpha is a given shape of the prior gamma distribution for lambda2.
#' @param beta is a given scale of the prior gamma distribution for lambda2.
#' @param theta1 is a given shape1 of the prior beta distribution for the probability of Tau and Kappa.
#' @param theta2 is a given shape2 of the prior beta distribution for the probability of Tau and Kappa.
#' @param xi2 is a given variance of the prior normal distribution for shifts.
#' @param method is a choice of methods including MT(McCulloch-Tsay), regression, LASSO, ALASSO(Adaptive LASSO), MonoLASSO(LASSO with Monotonicity constrains), MonoALASSO(Adaptive LASSO with Monotonicity constrains).
#' @param bound0 is an upper bound of the methods with Monotonicity constrains.
#' @param boundqplus1 is  a lower bound of the methods with Monotonicity constrains.
#' @param nsim is the number of draws from MCMC.
#' @param by is the interval of systematic sampling for the draws from MCMC.
#' @param burnin is the length of burn-in period.
#' @param tol is the tolerance level.
#' @param standardized is the flag triggering the standardization for the time series
#' @param logcc is the log transformation with continuity correction
#' @param FAP0 is the given false alarm probability
#' @param estimation.PPP is the estimation for Mu0, Phi and sigma2
#' @param nsim.PPP is the number of draws for PPP
#' 
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
#' result <- Ph1BayesianLASSO(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
Ph1BayesianLASSO <- function(Y, w = 28, H = NULL, X = NULL, Y0 = rep(mean(Y), w - 1), q = 5, 
                             A = diag(nrow = q + ifelse(is.null(X), 0, dim(X)[2])), 
                             a = 0.1, b = 0.1, alpha = 0.1, beta = 0.1, 
                             theta1 = 1, theta2 = 1, xi2 = 0.1,
                             method = "MonoALASSO", bound0 = Inf, boundqplus1 = 0,
                             nsim = 1000, by = 1, burnin = 1000, tol = 1e-10,
                             log = TRUE, const = 1, sta = TRUE, 
                             Y.hat.method = "median",
                             cc.method = "adjusted alpha",
                             FAP0 = 0.3, side = "two-sided", 
                             nsim.chart = 10000, tol.chart = 1e-6, 
                             plot = TRUE) {
  
  TT <- length(Y)
  
  model <- GibbsRFLSM.count(Y, w, H, X, Y0, q, 
    A, a, b, alpha, beta, 
    theta1, theta2, xi2,
    method, bound0, boundqplus1,
    nsim, by, burnin, tol,
    log, const, sta) 
  
  
  Y.tr.sim <- GibbsRFLSM.sim.ph1(nsim.chart, model$Y.tr, 
                                 model$Phi, model$muq, model$sigma2, 
                                 X, model$Beta, model$Kappa, 
                                 NULL, NULL, NULL)$Y.tr
  
  Y.hat <- rep(NA, TT - q)
  if (Y.hat.method == "median") {
    for (i in 1:(TT - q)) {
      Y.hat[i] <- median(Y.tr.sim[i, ])
    }
    #sigma2hat <- median(model$sigma2)
  } else if (Y.hat.method == "mean") {
    Y.hat <- rowMeans(Y.tr.sim)
    #sigma2hat <- mean(model$sigma2)
  }
  
  sigma2hat <- var(Y[(q + 1):TT] - Y.hat)
  
  lim.tr <- matrix(NA, nrow = TT, ncol = 2)
  
  if (cc.method == "adjusted alpha") {
    adjalpha <- adjalpha.ph1(Y.hat, sigma2hat, Y.tr.sim, FAP0, side, tol.chart)
    
    for (i in (q + 1):TT) {
      if (side == "two-sided") {
        lim.tr[i, ] <- quantile(adjalpha$resi[i - q, ], c(adjalpha$adjalpha / 2, 1 - adjalpha$adjalpha / 2)) * sigma2hat + Y.hat[i - q]
      } else if (side == "right-sided") {
        lim.tr[i, 1] <- -Inf
        lim.tr[i, 2] <- quantile(adjalpha$resi[i - q, ], c(1 - adjalpha$adjalpha)) * sigma2hat + Y.hat[i - q]
      } else if (side == "left-sided") {
        lim.tr[i, 1] <- quantile(adjalpha$resi[i - q, ], c(adjalpha$adjalpha)) * sigma2hat + Y.hat[i - q]
        lim.tr[i, 2] <- Inf
      }
    }
  } else if (cc.method == "classic") {
    cc <- cc.ph1(Y.hat, sigma2hat, Y.tr.sim, FAP0, side, tol.chart)
    
    for (i in (q + 1):TT) {
      if (side == "two-sided") {
        lim.tr[i, 1] <- Y.hat[i - q] - cc$cc * sigma2hat
        lim.tr[i, 2] <- Y.hat[i - q] + cc$cc * sigma2hat
      } else if (side == "right-sided") {
        lim.tr[i, 1] <- -Inf
        lim.tr[i, 2] <- Y.hat[i - q] + cc$cc * sigma2hat
      } else if (side == "left-sided") {
        lim.tr[i, 1] <- Y.hat[i - q] - cc$cc * sigma2hat
        lim.tr[i, 2] <- Inf
      }
    }
  }
  
  
  sig.tr <- (lim.tr[, 1] <= model$Y.tr) & (model$Y.tr <= lim.tr[, 2])
  
  if (plot == TRUE) {
    
    if (side == "two-sided") {
      Ylim <- c(min(lim.tr, model$Y.tr, na.rm = TRUE), max(lim.tr, model$Y.tr, na.rm = TRUE))
    } else if (side == "right-sided") {
      Ylim <- c(min(model$Y.tr, na.rm = TRUE), max(lim.tr, model$Y.tr, na.rm = TRUE))
    } else if (side == "left-sided") {
      Ylim <- c(min(lim.tr, model$Y.tr, na.rm = TRUE), max(model$Y.tr, na.rm = TRUE))
    }
    
    plot(c(1, TT), Ylim, type = 'n',
         main = "Phase I Chart for Transformed Moving Averages", 
         ylab = "Transformed Moving Averages", 
         xlab = "")
    points(model$Y.tr, type = 'o')
    points((1:TT)[which(sig.tr == FALSE)], model$Y.tr[which(sig.tr == FALSE)], col = 'red', pch = 16)
    points(lim.tr[, 1], type = 'l', lty = 2, col = 'red')
    points(lim.tr[, 2], type = 'l', lty = 2, col = 'red')
    
  }
  
  out <- list("model" = model, "cc" = cc$cc, "lim.tr" = lim.tr, 
              "sig.tr" = sig.tr, "sigma2hat" = sigma2hat) 
  out
} 


#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param Y is a vector.
#' @param H is the design matrix for shifts.
#' @param X is the input matrix
#' @param Y0 is the initial Y
#' @param q is the number of lags.
#' @param A is a given variance-covariance matrix in MT and regression for the slab-and-spike coefficients.
#' @param a is a given shape of the prior gamma distribution for sigma2.
#' @param b is a given scale of the prior gamma distribution for sigma2.
#' @param alpha is a given shape of the prior gamma distribution for lambda2.
#' @param beta is a given scale of the prior gamma distribution for lambda2.
#' @param theta1 is a given shape1 of the prior beta distribution for the probability of Tau and Kappa.
#' @param theta2 is a given shape2 of the prior beta distribution for the probability of Tau and Kappa.
#' @param xi2 is a given variance of the prior normal distribution for shifts.
#' @param method is a choice of methods including MT(McCulloch-Tsay), regression, LASSO, ALASSO(Adaptive LASSO), MonoLASSO(LASSO with Monotonicity constrains), MonoALASSO(Adaptive LASSO with Monotonicity constrains).
#' @param bound0 is an upper bound of the methods with Monotonicity constrains.
#' @param boundqplus1 is  a lower bound of the methods with Monotonicity constrains.
#' @param nsim is the number of draws from MCMC.
#' @param by is the interval of systematic sampling for the draws from MCMC.
#' @param burnin is the length of burn-in period.
#' @param tol is the tolerance level.
#' @param standardized is the flag triggering the standardization for the time series
#' @param logcc is the log transformation with continuity correction
#' @param FAP0 is the given false alarm probability
#' @param estimation.PPP is the estimation for Mu0, Phi and sigma2
#' @param nsim.PPP is the number of draws for PPP
#' 
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
#' result <- Ph1BayesianLASSO(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
Ph2BayesianLASSO.EWMA <- function(Y, Ph1BayesianLASSO.model, sigma2hat = 1, lambda = 0.05, w = 28, H = NULL, X = NULL,
                             Y1 = rep(0, dim(Ph1BayesianLASSO.model$Phi)[1] + w), X1 = NULL, H1 = NULL,
                             log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1,
                             Y.hat.method = "median",
                             cc.method = "adjusted alpha",
                             ARL0 = 360, side = "two-sided", max.length = 5000,
                             nsim.chart = 10000, tol.chart = 1e-6, 
                             plot = TRUE) {
  
  TT2 <- length(Y)
  q <- dim(Ph1BayesianLASSO.model$Phi)[1]
  
  Y2 <- c(Y1, Y)
  TT <- length(Y2)
  
  Y2.ma <- movaver(Y2, w)[(TT - TT2 - q + 1):TT]
  
  if (log == TRUE) {
    Y2.tr <- log(Y2.ma + const)
  }
  
  if (sta == TRUE) {
    Y2.tr <- (Y2.tr - meanY) / sdY
  }
  
  Y2.tr0 <- Y2.tr[1:q]
  Y2.tr1 <- Y2.tr[-c(1:q)]
  
  Y2.tr.sim <- GibbsRFLSM.sim.ph2(max.length, nsim.chart, 
                                  Ph1BayesianLASSO.model$Phi, Ph1BayesianLASSO.model$muq, Ph1BayesianLASSO.model$sigma2, 
                                  X, Ph1BayesianLASSO.model$Beta, Ph1BayesianLASSO.model$Kappa, 
                                  H, Ph1BayesianLASSO.model$Gamma, Ph1BayesianLASSO.model$Tau, 
                                  Y2.tr0, X1, H1)$Y.tr

  
  Y.hat <- rep(NA, TT2)
  if (Y.hat.method == "median") {
    for (i in 1:TT2) {
      Y.hat[i] <- median(Y2.tr.sim[i, ])
    }
    #sigma2hat <- median(Ph1BayesianLASSO.model$sigma2)
  } else if (Y.hat.method == "mean") {
    Y.hat <- rowMeans(Y2.tr.sim)
    #sigma2hat <- mean(Ph1BayesianLASSO.model$sigma2)
  }
  
  ewma <- (Y - Y.hat[1:TT2]) / sqrt(sigma2hat)
  for (i in 2:TT2) {
    ewma[i] <- lambda * ewma[i] + (1 - lambda) * ewma[i - 1]
  }
  
  
  lim.tr <- matrix(NA, nrow = TT2, ncol = 2)
  
  if (cc.method == "adjusted alpha") {
    adjalpha <- adjalpha.ph2(Y.hat, sigma2hat, Y2.tr.sim, ARL0, side, tol.chart)
    
    for (i in 1:TT2) {
      if (side == "two-sided") {
        lim.tr[i, ] <- quantile(adjalpha$ewma[i, ], c(adjalpha$adjalpha / 2, 1 - adjalpha$adjalpha / 2))
      } else if (side == "right-sided") {
        lim.tr[i, 1] <- -Inf
        lim.tr[i, 2] <- quantile(adjalpha$ewma[i, ], c(1 - adjalpha$adjalpha))
      } else if (side == "left-sided") {
        lim.tr[i, 1] <- quantile(adjalpha$ewma[i, ], c(adjalpha$adjalpha))
        lim.tr[i, 2] <- Inf
      }
    }
    
  } else if (cc.method == "classic") {
    cc <- cc.ph2(Y.hat, sigma2hat, Y2.tr.sim, ARL0, side, tol.chart)
    
    
    for (i in 1:TT2) {
      if (side == "two-sided") {
        lim.tr[i, 1] <- - cc$cc
        lim.tr[i, 2] <- cc$cc
      } else if (side == "right-sided") {
        lim.tr[i, 1] <- -Inf
        lim.tr[i, 2] <- cc$cc
      } else if (side == "left-sided") {
        lim.tr[i, 1] <- - cc$cc
        lim.tr[i, 2] <- Inf
      }
    }
    
  }
  
  sig.tr <- (lim.tr[, 1] <= ewma) & (ewma <= lim.tr[, 2])
  
  if (plot == TRUE) {
    
    if (side == "two-sided") {
      Ylim <- c(min(lim.tr, ewma, na.rm = TRUE), max(lim.tr, ewma, na.rm = TRUE))
    } else if (side == "right-sided") {
      Ylim <- c(min(ewma, na.rm = TRUE), max(lim.tr, ewma, na.rm = TRUE))
    } else if (side == "left-sided") {
      Ylim <- c(min(lim.tr, ewma, na.rm = TRUE), max(ewma, na.rm = TRUE))
    }
    
    plot(c(1, TT2), Ylim, type = 'n',
         main = "Phase II Chart for Transformed Moving Averages", 
         ylab = "EWMA", 
         xlab = "")
    points(ewma, type = 'o')
    points((1:TT2)[which(sig.tr == FALSE)], ewma[which(sig.tr == FALSE)], col = 'red', pch = 16)
    points(lim.tr[, 1], type = 'l', lty = 2, col = 'red')
    points(lim.tr[, 2], type = 'l', lty = 2, col = 'red')
    
  }
  
  out <- list("EWMA" = ewma, "cc" = cc$cc, "lim.tr" = lim.tr, 
              "sig.tr" = sig.tr) 
  out
} 
