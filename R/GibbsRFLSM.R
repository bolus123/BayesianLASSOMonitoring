#' Bayesian LASSO Phase I Monitoring
#' 
#' @export
#' 
#' 
bset <- function(X = NULL, method = "ALASSO", phimono = TRUE, phiq = 5, phiA = diag(nrow = phiq), phibound0 = Inf, phiboundqplus1 = 0, 
                 betaA = ifelse(class(X)[1] == "matrix", dim(X)[2], 1), gammaxi2 = 0.1, tautheta1 = 1, tautheta2 = 1, 
                 sigma2a = 1, sigma2b = 1, lambda2 = NULL, updatelambda2 = TRUE, lambda2alpha = 0.1, lambda2beta = 0.1, 
                 theta = NULL, YJ = TRUE, updateYJ = TRUE, leftcensoring = TRUE, lowerbound = 0, rounding = TRUE) {
  bset <- list(
    "method" = method,
    "phimono" = ifelse(phimono, 1, 0),
    "phiq" = phiq,
    "phiA" = phiA,
    "phibound0" = phibound0,
    "phiboundqplus1" = phiboundqplus1,
    "betaA" = betaA,
    "gammaxi2" = gammaxi2,
    "tautheta1" = tautheta1,
    "tautheta2" = tautheta2,
    "sigma2a" = sigma2a,
    "sigma2b" = sigma2b,
    "lambda2" = lambda2,
    "updatelambda2" = ifelse(updatelambda2, 1, 0),
    "lambda2alpha" = lambda2alpha,
    "lambda2beta" = lambda2beta,
    "theta" = theta,
    "YJ" = ifelse(YJ, 1, 0),
    "updateYJ" = ifelse(updateYJ, 1, 0),
    "leftcensoring" = ifelse(leftcensoring, 1, 0),
    "lowerbound" = lowerbound,
    "rounding" = ifelse(rounding, 1, 0)
  )
  return(bset)
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
GibbsRFLSMX <- function(Y, X = NULL, H = NULL, 
                        bset = bset(betap = ifelse(class(X)[1] == "matrix", diag(nrow = dim(X)[2]), diag(nrow = 1))), 
                        tol = 1e-10, nsim = 300, thin = 10, burnin = 1000, verbose = TRUE) {
  
  model <- GibbsRFLSMXcpp(Y, 
                 bset, tol, 
                 nsim, thin, burnin, 
                 verbose = ifelse(verbose, 1, 0),
                 X = X,
                 H = H,
                 lambda2 = bset$lambda2,
                 theta = bset$theta)
  
  return(model)
  
}