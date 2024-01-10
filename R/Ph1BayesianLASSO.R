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
                             FAP0 = 0.3, side = "two-sided", 
                             interval = c(0.00001, 0.4), nsim.chart = 10000, tol.chart = 1e-6, 
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
  
  Y.ma.sim <- backtrans(Y.tr.sim, log, const, sta, model$meanY, model$sdY)
  Y.ma.sim[Y.ma.sim < 0] <- 0
  
  adja.tr <- adjalp(Y.tr.sim, FAP0, side, interval, tol.chart)
  adja.ma <- adjalp(Y.ma.sim, FAP0, side, interval, tol.chart)
  
  lim.tr <- matrix(NA, nrow = TT, ncol = 2)
  lim.ma <- lim.tr
  
  for (i in (q + 1):TT) {
    if (side == "two-sided") {
      lim.tr[i, ] <- quantile(Y.tr.sim[i - q, ], c(adja.tr / 2, 1 - adja.tr / 2))
      lim.ma[i, ] <- quantile(Y.ma.sim[i - q, ], c(adja.ma / 2, 1 - adja.ma / 2))
    } else if (side == "right-sided") {
      lim.tr[i, 1] <- -Inf
      lim.ma[i, 1] <- -Inf
      lim.tr[i, 2] <- quantile(Y.tr.sim[i - q, ], c(1 - adja.tr))
      lim.ma[i, 2] <- quantile(Y.ma.sim[i - q, ], c(1 - adja.ma))
    } else if (side == "left-sided") {
      lim.tr[i, 1] <- quantile(Y.tr.sim[i - q, ], c(adja.tr))
      lim.ma[i, 1] <- quantile(Y.ma.sim[i - q, ], c(adja.ma))
      lim.tr[i, 2] <- Inf
      lim.ma[i, 2] <- Inf
    }
  }
  
  sig.tr <- (lim.tr[, 1] <= model$Y.tr) & (model$Y.tr <= lim.tr[, 2])
  sig.ma <- (lim.ma[, 1] <= model$Y.ma) & (model$Y.ma <= lim.ma[, 2])
  
  if (plot == TRUE) {
    plot(c(1, TT), c(min(lim.tr, model$Y.tr), max(lim.tr, model$Y.tr)), type = 'n',
         main = "Phase I Chart for Transformed Moving Averages", 
         ylab = "Transformed Moving Averages", 
         xlab = "")
    points(model$Y.tr, type = 'o')
    points(lim.tr[, 1], type = 'l', lty = 2, col = 'red')
    points(lim.tr[, 2], type = 'l', lty = 2, col = 'red')
    
    plot(c(1, TT), c(min(lim.ma, model$Y.ma), max(lim.ma, model$Y.ma)), type = 'n',
         main = "Phase I Chart for Moving Averages", 
         ylab = "Moving Averages", 
         xlab = "")
    points(model$Y.ma, type = 'o')
    points(lim.ma[, 1], type = 'l', lty = 2, col = 'red')
    points(lim.ma[, 2], type = 'l', lty = 2, col = 'red')
  }
  
  out <- list("model" = model, "sig.tr" = sig.tr, "sig.ma" = sig.ma) 
  out
} 
