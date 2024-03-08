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
GibbsRFLSMX <- function(Y, bset, X = NULL, H = NULL, 
                        tol = 1e-10, nsim = 300, thin = 10, burnin = 1000, verbose = TRUE) {
  
  model <- GibbsRFLSMXcpp(matrix(Y, ncol = 1), bset, 
                        tol = tol, nsim = nsim, thin = thin, burnin = burnin, verbose = verbose,
                        X = X, H = H, lambda2 = bset$lambda2, theta = bset$theta)
  
  out <- list(
    "Phi" = model$Phi,
    "Beta" = model$Beta,
    "Gamma" = model$Gamma,
    "Tau" = model$Tau,
    "mu0" = model$mu0,
    "sigma2" = model$sigma2,
    "lambda2" = model$lambda2,
    "theta" = model$theta,
    "Z" = model$Z,
    "H" = H,
    "X" = X,
    "Y" = Y,
    "nsim" = nsim
  )
  
  return(out)
  
}