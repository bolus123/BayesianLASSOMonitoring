#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param Y is a vector.
#' @param H is the design matrix for shifts.
#' @param X is the input matrix
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
#' result <- BayesianLASSOPh1(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
BayesianLASSOPh1 <- function(Y, H = NULL, X = NULL, q = 5, 
                             A = diag(nrow = q + ifelse(is.null(X), 0, dim(X)[2])), 
                             a = 0.1, b = 0.1, alpha = 0.1, beta = 0.1, 
                             theta1 = 1, theta2 = 1, xi2 = 0.1,
                             method = "MonoALASSO", bound0 = Inf, boundqplus1 = 0,
                             nsim = 1000, by = 1, burnin = 1000, tol = 1e-10, 
                             standardized = TRUE, logcc = FALSE,
                             FAP0 = 0.2, estimation.PPP = "median", nsim.PPP = 1000) {
  
  TT <- length(Y)
  
  if (logcc == TRUE) {
    Y1 <- log(Y + 0.5)
  } else {
    Y1 <- Y
  }
  
  if (standardized == TRUE) {
    meanY1 <- mean(Y1)
    sdY1 <- sd(Y1)
    Y1 <- (Y1 - meanY1) / sdY1
  } 
  
  model <- GibbsRFLSM(Y1, H, X, q,      
                      A, 
                      a, b, alpha, beta, 
                      theta1, theta2, xi2,
                      method, bound0, boundqplus1,
                      nsim, by, burnin, tol)
  
  Mu0 <- matrix(NA, nrow = TT, ncol = nsim)
  for (j in 1:nsim) {
    Mu0[, j] <- model$muq[j]
    if (!is.null(X)) {
      Mu0[, j] <- X %*% (model$Beta[, j] * model$Kappa[, j])
    }
  }
  
  if (estimation.PPP == "median") {
    Mu0hat <- rep(NA, TT)
    for (i in 1:TT) {
      Mu0hat[i] <- median(Mu0[i, ])
    }
    
    Phihat <- rep(NA, q)
    for (i in 1:q) {
      Phihat[i] <- median(model$Phi[i, ])
    }
    
    sigma2hat <- median(model$sigma2)
  } else if (estimation.PPP == "mean") {
    Mu0hat <- rep(NA, TT)
    for (i in 1:TT) {
      Mu0hat[i] <- mean(Mu0[i, ])
    }
    
    Phihat <- rep(NA, q)
    for (i in 1:q) {
      Phihat[i] <- mean(model$Phi[i, ])
    }
    
    sigma2hat <- mean(model$sigma2)
  }
  
  chart <- GibbsRFLSM.CC.PPP.residual(Y1, model$Phi, Mu0, model$sigma2, 
                                         Phihat, Mu0hat, sigma2hat, FAP0, nsim.PPP)
  
  if (standardized == TRUE) {
    lowerbound <- chart$lowerbound * sdY1 + meanY1
    upperbound <- chart$upperbound * sdY1 + meanY1
  } else {
    lowerbound <- chart$lowerbound
    upperbound <- chart$upperbound
  }
  
  if (logcc == TRUE) {
    lowerbound <- exp(lowerbound) - 0.5
    upperbound <- exp(upperbound) - 0.5
  }
  
  sig <- rep(NA, TT)
  
  for (i in seq(TT)) {
    sig[i] <- ifelse((Y[i] < lowerbound[i]) || (upperbound[i] < Y[i]), 1, 0)
  }
  
  out <- list("lowerbound" = lowerbound, "upperbound" = upperbound, "sig" = sig, 
              "Omni" = chart$Omni, "Ind" = chart$Ind, "model" = model, "Obs" = Y1, 
              "lowerboundtr" = chart$lowerbound, "upperboundtr" = chart$upperbound)
  
} 
