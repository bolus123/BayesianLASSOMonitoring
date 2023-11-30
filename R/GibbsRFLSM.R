#' Random Flexible Level Shift Model
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
#' @references McCulloch, R. E., & Tsay, R. S. (1993). Bayesian inference and prediction for mean and variance shifts in autoregressive time series. Journal of the american Statistical association, 88(423), 968-978.
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
#' result <- GibbsRFLSM(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
GibbsRFLSM <- function(Y, H = NULL, X = NULL, q = 5, 
                       A = diag(nrow = q + ifelse(is.null(X), 0, dim(X)[2])), 
                       a = 0.1, b = 0.1, alpha = 0.1, beta = 0.1, 
                       theta1 = 1, theta2 = 1, xi2 = 0.1,
                       method = "MonoALASSO", bound0 = Inf, boundqplus1 = 0,
                       nsim = 1000, by = 1, burnin = 1000, tol = 1e-10) {
  
  TT <- length(Y)
  
  if (is.null(H) && is.null(X)) {
    model <- GibbsRFLSMcpp(Y, q, 
                       A, a, b, alpha, beta, 
                       theta1, theta2, xi2,
                       method, bound0, boundqplus1,
                       nsim, by, burnin,
                       tol)
  } else {
    H1 <- cbind(H, X)
    model <- GibbsRFLSMcpp(Y, q, 
                        A, a, b, alpha, beta, 
                        theta1, theta2, xi2,
                        method, bound0, boundqplus1,
                        nsim, by, burnin,
                        tol, H1)
    
  }
  
  if (is.null(H)) {
    m <- 0
    Gamma <- NA
    Tau <- NA
    pGamma <- NA
  } else {
    m <- dim(H)[2]
    Gamma <- model$Gamma[1:m, ]
    Tau <- model$Tau[1:m, ]
    pGamma <- model$p[1:m, ]
  }
  
  if (is.null(X)) {
    p <- 0
    Beta <- NA
    Kappa <- NA
    pBeta <- NA
  } else {
    p <- dim(X)[2]
    Beta <- model$Gamma[(m + 1):(m + p), ]
    Kappa <- model$Tau[(m + 1):(m + p), ]
    pBeta <- model$p[(m + 1):(m + p), ]
  }
  
  
  fit <- matrix(NA, nrow = TT, ncol = nsim)
  resi <- matrix(NA, nrow = TT, ncol = nsim)
  stdresi <- matrix(NA, nrow = TT, ncol = nsim)
  
  for (i in 1:nsim) {
    
    fit[, i] <- model$Mu[, i]
    
    tmpresi <- Y - model$Mu[, i]
    tmpV <- getV(tmpresi, q)
    V <- tmpV[-c(1:q), ]
    
    fit[(q + 1):TT, i] <- fit[(q + 1):TT, i] + V %*% model$Phi[, i]
    
    resi[, i] <- Y - fit[, i]

  }
  
  out <- list(
    "Phi" = model$Phi,
    "Beta" = Beta,
    "pBeta" = pBeta,
    "Kappa" = Kappa,
    "Gamma" = Gamma,
    "pGamma" = pGamma,
    "Tau" = Tau,
    "sigma2" = model$sigma2,
    "lambda2" = model$lambda2,
    "muq" = model$muq,
    "Mu" = model$Mu,
    "fit" = fit,
    "resi" = resi
  )
  
  return(out)
  
}

movaver <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 1)}


#' Random Flexible Level Shift Model
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
#' @param logcc is the log transformation
#' @param standardized is the standardization
#' @references McCulloch, R. E., & Tsay, R. S. (1993). Bayesian inference and prediction for mean and variance shifts in autoregressive time series. Journal of the american Statistical association, 88(423), 968-978.
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
#' result <- GibbsRFLSM.count(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
GibbsRFLSM.count <- function(Y, w = 28, H = NULL, X = NULL, Y0 = rep(mean(Y), w - 1), q = 5, 
                       A = diag(nrow = q + ifelse(is.null(X), 0, dim(X)[2])), 
                       a = 0.1, b = 0.1, alpha = 0.1, beta = 0.1, 
                       theta1 = 1, theta2 = 1, xi2 = 0.1,
                       method = "MonoALASSO", bound0 = Inf, boundqplus1 = 0,
                       nsim = 1000, by = 1, burnin = 1000, tol = 1e-10,
                       logcc = TRUE, standardized = TRUE) {
  
  YY <- c(Y0, Y)
  TT <- length(Y)
  nn <- length(YY)
  ####################################
  
    Y1 <- movaver(YY, w)[(nn - TT + 1):nn]
    Y2 <- Y1
    
  ####################################
  
  if (logcc == TRUE) {
    Y1 <- log(Y1 + 0.5)
  }
  
  if (standardized == TRUE) {
    meanY <- mean(Y1)
    sdY <- sd(Y1)
    Y1 <- (Y1 - meanY) / sdY
  }
  
  model <- GibbsRFLSM(Y1, H, X, q, 
                    A, 
                    a, b, alpha, beta, 
                    theta1, theta2, xi2,
                    method, bound0, boundqplus1,
                    nsim, by, burnin, tol)
    
  
  fit0 <- model$fit
  
  if (standardized == TRUE) {
    fit0 <- fit0 * sdY + meanY
  }
  
  if (logcc == TRUE) {
    fit0 <- exp(fit0) - 0.5
  }
  
  fit <- fit0 * w
 
  for (i in 1:TT) {
    fit[i, ] <- fit[i, ] - sum(YY[(w - 1 + i - 1):((w - 1) + i - (w - 1))])
  }
  
  out <- list(
    "Phi" = model$Phi,
    "Beta" = model$Beta,
    "pBeta" = model$pBeta,
    "Kappa" = model$Kappa,
    "Gamma" = model$Gamma,
    "pGamma" = model$pGamma,
    "Tau" = model$Tau,
    "sigma2" = model$sigma2,
    "lambda2" = model$lambda2,
    "muq" = model$muq,
    "Mu" = model$Mu,
    "fit.tr" = model$fit,
    "resi.tr" = model$resi,
    "Y.tr" = Y1,
    "fit.ma" = fit0,
    "resi.ma" = matrix(Y2, ncol = nsim, nrow = length(Y)) - fit0,
    "Y.ma" = Y2,
    "fit" = fit,
    "resi" = matrix(Y, ncol = nsim, nrow = length(Y)) - fit,
    "Y" = Y
  )
  
  return(out)
  
}