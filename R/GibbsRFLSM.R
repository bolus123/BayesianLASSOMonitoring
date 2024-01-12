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
    "Mu" = model$Mu
  )
  
  return(out)
  
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
trans <- function(Y, log = TRUE, const = 1, sta = TRUE){
  out <- Y
  if (log == TRUE) {
    out <- log(out + const)
  }
  meanY <- 0
  sdY <- 1
  if (sta == TRUE) {
    meanY <- mean(out)
    sdY <- sd(out)
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
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization
#' @references McCulloch, R. E., & Tsay, R. S. (1993). Bayesian inference and prediction for mean and variance shifts in autoregressive time series. Journal of the american Statistical association, 88(423), 968-978.
#' 
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
#' 
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM.count(Y, H = H, q = q, nsim = nsim, burnin = burnin)
#' 
GibbsRFLSM.count <- function(Y, w = 7, H = NULL, X = NULL, Y0 = rep(mean(Y), w - 1), q = 5, 
                       A = diag(nrow = q + ifelse(is.null(X), 0, dim(X)[2])), 
                       a = 0.1, b = 0.1, alpha = 0.1, beta = 0.1, 
                       theta1 = 1, theta2 = 1, xi2 = 0.1,
                       method = "MonoALASSO", bound0 = Inf, boundqplus1 = 0,
                       nsim = 1000, by = 1, burnin = 1000, tol = 1e-10,
                       log = TRUE, const = 1, sta = TRUE) {
  
  YY <- c(Y0, Y)
  TT <- length(Y)
  nn <- length(YY)
  ####################################
  
    Y1 <- movaver(YY, w)[(nn - TT + 1):nn]
    Y1.ma <- Y1
    Y2 <- trans(Y1, log = log, const = const, sta = sta)
    Y1 <- Y2$Y
  
  ####################################
  
  model <- GibbsRFLSM(Y1, H, X, q, 
                    A, 
                    a, b, alpha, beta, 
                    theta1, theta2, xi2,
                    method, bound0, boundqplus1,
                    nsim, by, burnin, tol)
    
  
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
    "Y.tr" = Y2$Y,
    "meanY" = Y2$meanY,
    "sdY" = Y2$sdY,
    "Y.ma" = Y1.ma
  )
  
  return(out)
  
}