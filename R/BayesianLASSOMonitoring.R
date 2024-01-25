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
Ph1BayesianLASSO <- function(Y, w = 7, H = NULL, X = NULL, Y0 = rep(mean(Y), w - 1), q = 5, 
                             A = diag(nrow = q), 
                             a = 0.1, b = 0.1, alpha = 0.1, beta = 0.1, 
                             theta1 = 1, theta2 = 1, xi2 = 0.1,
                             method = "MonoALASSO", bound0 = Inf, boundqplus1 = 0,
                             nsim = 1000, by = 1, burnin = 1000, tol = 1e-10,
                             log = TRUE, const = 1, sta = TRUE, 
                             sign.method = "DM",
                             adj.method = "holm",
                             FAP0 = 0.3, side = "two-sided",  
                             plot = TRUE) {
  
  TT <- length(Y)
  
  model <- GibbsRFLSM.count(Y, w, H, X, Y0, q, 
    A, a, b, alpha, beta, 
    theta1, theta2, xi2,
    method, bound0, boundqplus1,
    nsim, by, burnin, tol,
    log, const, sta) 
  
  mu <- model$H %*% (model$Gamma * model$Tau)
  
  mu.dif <- matrix(NA, nrow = TT - q, ncol = nsim)
  
  for (i in (q + 1):TT) {
    if (i == q + 1) {
      mu.dif[i - q, ] <- mu[i, ]
    } else {
      mu.dif[i - q, ] <- mu[i, ] - mu[i - 1, ]
    }
  }
  
  sign <- matrix(NA, nrow = TT - q, ncol = 3)
  sign[, 1] <- rowSums(mu.dif > 0)
  sign[, 2] <- rowSums(mu.dif == 0)
  sign[, 3] <- rowSums(mu.dif < 0)
  
  pvalue <- rep(NA, TT - q)
  
  if (sign.method == "DM") {
    for (i in 1:(TT - q)) {
      pvaluetmp <- pbinom(sign[i, 1] + sign[i, 2] / 2, nsim, 0.5)
      
      if (side == "left-sided") {
        pvalue[i] <- pvaluetmp
      } else if (side == "right-sided") {
        pvalue[i] <- 1 - pvaluetmp
      } else if (side == "two-sided") {
        pvalue[i] <- 2 * min(1 - pvaluetmp, pvaluetmp)
      }
      
    }
  } else if (sign.method == "trinomial") {
    
    tmp <- cbind(sign[, 1] - sign[, 3], nsim, sign[, 2] / nsim)
    
    for (i in 1:(TT - q)) {
      if (tmp[i, 3] == 1) {
        pvalue[i] <- 1
      } else {
        pvaluetmp <- ptrinomial(tmp[i, 1], tmp[i, 2], tmp[i, 3])
        if (side == "left-sided") {
          pvalue[i] <- pvaluetmp
        } else if (side == "right-sided") {
          pvalue[i] <- 1 - pvaluetmp
        } else if (side == "two-sided") {
          pvalue[i] <- 2 * min(1 - pvaluetmp, pvaluetmp)
        }
      }
    }
  }
  
  adj.pvalue <- p.adjust(pvalue, method = adj.method)
  
  sig <- adj.pvalue < FAP0
  
  if (plot == TRUE) {   
    
    #if (side == "two-sided") {
    #  Ylim <- c(min(lim.tr, model$Y.tr, na.rm = TRUE), max(lim.tr, model$Y.tr, na.rm = TRUE))
    #} else if (side == "right-sided") {
    #  Ylim <- c(min(model$Y.tr, na.rm = TRUE), max(lim.tr, model$Y.tr, na.rm = TRUE))
    #} else if (side == "left-sided") {
    #  Ylim <- c(min(lim.tr, model$Y.tr, na.rm = TRUE), max(model$Y.tr, na.rm = TRUE))
    #}
    
    Ylim <- c(min(FAP0, adj.pvalue, na.rm = TRUE), 
              max(FAP0, adj.pvalue, na.rm = TRUE))
    
    plot(c(1, TT), Ylim, type = 'n',
         main = "Phase I Chart", 
         ylab = "Adjusted P-Value", 
         xlab = "")
    points(adj.pvalue, type = 'o')
    points((1:(TT - q))[which(sig == TRUE)], adj.pvalue[which(sig == TRUE)], col = 'red', pch = 16)
    abline(h = FAP0, col = 'red')
    
  }
  
  out <- list("model" = model, "adj.pvalue" = adj.pvalue, 
              "adj.pvalue" = adj.pvalue) 
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
Ph2BayesianLASSO.EWMA <- function(Y, Ph1BayesianLASSO.chart, lambda = 0.05, H = NULL, X = NULL,
                             Y1 = rep(0, dim(Ph1BayesianLASSO.model$Phi)[1] + w), X1 = NULL, H1 = NULL,
                             log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1,
                             Y.hat.method = "median",
                             cc.method = "adjusted alpha",
                             ARL0 = 360, side = "two-sided", max.length = 5000,
                             nsim.chart = 10000, tol.chart = 1e-6, 
                             plot = TRUE) {
  
  Ph1BayesianLASSO.model <- Ph1BayesianLASSO.chart$model
  
  w <- Ph1BayesianLASSO.model$w
  
  TT2 <- length(Y)
  q <- dim(Ph1BayesianLASSO.model$Phi)[1]
  
  Y2 <- c(Y1, Y)
  TT <- length(Y2)
  
  Y2.ma <- movaver(Y2, w)[(TT - TT2 - q + 1):TT]
  Y2.tr <- Y2.ma
  if (log == TRUE) {
    Y2.tr <- log(Y2.tr + const)
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
                                  Y2.tr0, X1, H1)

  Y2.hat.sim <- Y2.tr.sim$fit
  Y2.tr.sim <- Y2.tr.sim$Y.tr
  
  Y.hat.overall <- rep(NA, max.length)
  Y.hat <- rep(NA, TT2)
  if (Y.hat.method == "median") {
    for (i in 1:max.length) {
      Y.hat.overall[i] <- median(Y2.hat.sim[i, ])
    }
    sigma2hat <- median(Ph1BayesianLASSO.model$sigma2)
  } else if (Y.hat.method == "mean") {
    Y.hat.overall <- rowMeans(Y2.hat.sim)
    sigma2hat <- mean(Ph1BayesianLASSO.model$sigma2)
  }
  Y.hat <- Y.hat.overall[1:TT2]
  sigmahat <- sqrt(sigma2hat)
  
  lim.tr <- matrix(NA, nrow = TT2, ncol = 2)
  sig.tr <- lim.tr[, 1]
  
  ewma <- (Y - Y.hat) / sigmahat
  ewma <- (ewma - Ph1BayesianLASSO.chart$cs.mean) / Ph1BayesianLASSO.chart$cs.sd
  for (i in 2:TT2) {
    ewma[i] <- lambda * ewma[i] + (1 - lambda) * ewma[i - 1]
  }
  
  if (plot == TRUE) {
  
    
  if (cc.method == "adjusted alpha") {
    adjalpha <- adjalpha.ph2(Y.hat.overall, sigma2hat, Y2.tr.sim, 
                             Ph1BayesianLASSO.chart$cs.mean, Ph1BayesianLASSO.chart$cs.sd, 
                             ARL0, side, tol.chart)
    
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
    cc <- cc.ph2(Y.hat.overall, sigma2hat, Y2.tr.sim, 
                 Ph1BayesianLASSO.chart$cs.mean, Ph1BayesianLASSO.chart$cs.sd, 
                 ARL0, side, tol.chart)
    
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
  
    if (side == "two-sided") {
      Ylim <- c(min(lim.tr, ewma, na.rm = TRUE), max(lim.tr, ewma, na.rm = TRUE))
    } else if (side == "right-sided") {
      Ylim <- c(min(ewma, na.rm = TRUE), max(lim.tr, ewma, na.rm = TRUE))
    } else if (side == "left-sided") {
      Ylim <- c(min(lim.tr, ewma, na.rm = TRUE), max(ewma, na.rm = TRUE))
    }
    
    plot(c(1, TT2), Ylim, type = 'n',
         main = "Phase II Chart", 
         ylab = "EWMA", 
         xlab = "")
    points(ewma, type = 'o')
    points((1:TT2)[which(sig.tr == FALSE)], ewma[which(sig.tr == FALSE)], col = 'red', pch = 16)
    points(lim.tr[, 1], type = 'l', lty = 2, col = 'red')
    points(lim.tr[, 2], type = 'l', lty = 2, col = 'red')
    
  }
  
  out <- list("EWMA" = ewma, "lim.tr" = lim.tr, 
              "sig.tr" = sig.tr, "Y.hat" = Y.hat) 
  out
} 
