#' get a posterior sample using gibbs sampling for Random Flexible Level Shift Model
#'
#' @param Y is a vector.
#' @param H is the design matrix for shifts.
#' @param X is the input matrix
#' @param q is the number of lags.
#' @param A is a given variance-covariance matrix in MT and regression for the coefficients for autoregressors.
#' @param a is a given shape of the prior gamma distribution for sigma2.
#' @param b is a given scale of the prior gamma distribution for sigma2.
#' @param alpha is a given shape of the prior gamma distribution for lambda2.
#' @param beta is a given scale of the prior gamma distribution for lambda22.
#' @param theta1 is a given shape1 of the prior beta distribution for the probability of Tau.
#' @param theta2 is a given shape2 of the prior beta distribution for the probability of Tau.
#' @param xi2 is a given variance of the prior normal distribution for shifts.
#' @param method is a choice of methods including MT(McCulloch-Tsay), regression, LASSO, ALASSO(Adaptive LASSO), MonoLASSO(LASSO with Monotonicity constrains), MonoALASSO(Adaptive LASSO with Monotonicity constrains).
#' @param bound0 is a lower bound of the methods with Monotonicity constrains.
#' @param boundqplus1 is  a upper bound of the methods with Monotonicity constrains.
#' @param nsim is the number of samples draw from MCMC.
#' @param by is the interval of systematic sampling for the draws from MCMC.
#' @param burnin is the length of burn-in period.
#' @param tol is the tolerance.
#' @export
#' @examples
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, H = H)
GibbsRFLSM <- function(Y, H = NULL, X = NULL, q = 5, A = diag(nrow = q), 
                       a = 0.1, b = 0.1, alpha = 0.1, beta = 0.1, 
                       theta1 = 1, theta2 = 1, xi2 = 0.1,
                       method = "MonoALASSO", bound0 = Inf, boundqplus1 = 0,
                       nsim = 1000, by = 1, burnin = 1000, tol = 1e-10) {
  
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
    
  }
  
  out <- list(
    "Phi" = model$Phi,
    "Beta" = Beta,
    "pBeta" = pBeta,
    "Kappa" = Kappa,
    "Gamma" = Gamma,
    "pGamma" = pGamma,
    "Tau" = Tau,
    "lambda2" = model$lambda2,
    "muq" = model$muq,
    "Mu" = model$Mu
  )
  
  return(out)
  
}