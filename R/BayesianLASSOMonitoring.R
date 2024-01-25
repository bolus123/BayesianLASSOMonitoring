#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param model is model.
#' @param sign.method is .
#' @param adj.method is 
#' @param side is side
#' 
#' 
#' @export
Ph1MultipleTesting <- function(model, sign.method = "DM",
                               adj.method = "holm",
                               side = "two-sided") {
  
  TT <- length(model$Y)
  q <- dim(model$Phi)[1]
  nsim <- dim(model$Phi)[2]
  
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
      pvalue[i] <- pbinom(sign[i, 1] + sign[i, 2] / 2, nsim, 0.5)
    }
  } else if (sign.method == "trinomial") {
    
    tmp <- cbind(sign[, 1] - sign[, 3], nsim, sign[, 2] / nsim)
    
    for (i in 1:(TT - q)) {
      if (tmp[i, 3] == 1) {
        pvalue[i] <- 1
      } else {
        pvalue[i] <- ptrinomial(tmp[i, 1], tmp[i, 2], tmp[i, 3])
      }
    }
  }
  
  if (side == "left-sided") {
    pvalue <- pvalue
  } else if (side == "right-sided") {
    pvalue <- 1 - pvalue
  } else if (side == "two-sided") {
    for (i in 1:(TT - q)) {
      pvalue[i] <- 2 * min(1 - pvalue[i], pvalue[i])
    }
  }
  
  adj.pvalue <- p.adjust(pvalue, method = adj.method)
  adj.pvalue
   
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
  
  model <- GibbsRFLSM.ma(Y, w, H, X, Y0, q, 
    A, a, b, alpha, beta, 
    theta1, theta2, xi2,
    method, bound0, boundqplus1,
    nsim, by, burnin, tol,
    log, const, sta) 
  
  adj.pvalue <- Ph1MultipleTesting(model, sign.method, adj.method, side)
  
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
              "sig" = sig) 
  out
} 

