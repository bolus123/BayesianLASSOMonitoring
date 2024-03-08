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
  
  m <- dim(model$Gamma)[1]
  
  GammaTau <- model$Gamma * model$Tau
  
  
  sign <- matrix(NA, nrow = m, ncol = 3)
  sign[, 1] <- rowSums(GammaTau > 0)
  sign[, 2] <- rowSums(GammaTau == 0)
  sign[, 3] <- rowSums(GammaTau < 0)
  
  pvalue <- rep(NA, m)
  
  if (sign.method == "DM") {
    for (i in 1:m) {
      pvalue[i] <- pbinom(sign[i, 1] + sign[i, 2] / 2, nsim, 0.5)
    }
  } else if (sign.method == "trinomial") {
    
    tmp <- cbind(sign[, 1] - sign[, 3], nsim, sign[, 2] / nsim)
    
    for (i in 1:m) {
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
    for (i in 1:m) {
      pvalue[i] <- 2 * min(1 - pvalue[i], pvalue[i])
    }
  }
  
  adj.pvalue <- p.adjust(pvalue, method = adj.method)
  adj.pvalue
   
}

getlogBF <- function(Y, model) {
  m <- dim(model$Phi)[2]
  q <- dim(model$Phi)[1]
  n <- length(Y) - q
  
  tmp0 <- rep(0, n)
  tmp1 <- rep(0, n)
  
  for (i in 1:m) {
    Y0.sim.fit0 <- fit.GibbsRFLSM(Y, 
                                  model$Phi[, i], model$muq[i], 
                                  X = model$X, Beta = model$Beta[, i], 
                                  Kappa = model$Kappa[, i], 
                                  H = NULL, Gamma = NULL, Tau = NULL)
    
    Y0.sim.fit1 <- fit.GibbsRFLSM(Y, 
                                  model$Phi[, i], model$muq[i], 
                                  X = model$X, Beta = model$Beta[, i], 
                                  Kappa = model$Kappa[, i], 
                                  H = model$H, Gamma = model$Gamma[, i], 
                                  Tau = model$Tau[, i])
    
    tmp0 <- tmp0 + dnorm(Y[-c(1:q)], Y0.sim.fit0, sd = sqrt(model$sigma2[i])) / m
    tmp1 <- tmp1 + dnorm(Y[-c(1:q)], Y0.sim.fit1, sd = sqrt(model$sigma2[i])) / m
  }
  
  b10 <- log(tmp1 / tmp0)
  b10
}


#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param model is model.
#' @param nsim is .
#' @param FAP0 is 
#' 
#' 
#' @export
Ph1MultipleTesting.BF <- function(model, nsim = 1000, FAP0 = 0.2) {
  
  q <- dim(model$Phi)[1]
  n <- length(model$Y.tr) - q
  m <- dim(model$Phi)[2]
  
  b10.matrix <- matrix(NA, nrow = n, ncol = nsim)
  b10.max <- rep(NA, nsim)
  
  for (j in 1:nsim) {
    
    tmpsel <- sample(1:m, 1)
    Y0fit <- fit.GibbsRFLSM(model$Y.tr, 
                            model$Phi[, tmpsel], model$muq[tmpsel], 
                            X = model$X, Beta = model$Beta[, tmpsel], 
                            Kappa = model$Kappa[, tmpsel], 
                            H = NULL, Gamma = NULL, Tau = NULL)
    
    Y0.sim <- rnorm(n, Y0fit, sqrt(model$sigma2[tmpsel]))
    
    tmpY.sim <- c(model$Y.tr[1:q], Y0.sim)
    
    b10 <- getlogBF(tmpY.sim, model)
    
    b10.matrix[, j] <- b10
    b10.max[j] <- max(b10, na.rm =TRUE)
    
  }
  
  cs <- getlogBF(model$Y.tr, model)
  
  p.value <- mean(b10.max > max(cs, na.rm = TRUE), na.rm = TRUE)
  
  lim <- quantile(b10.max, 1 - FAP0, na.rm = TRUE)
  
  list(p.value = p.value, lim = lim, cs = cs, sig = cs > lim)
  
}


#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param model is model.
#' @param nsim is .
#' @param FAP0 is 
#' @param log is model.
#' @param const is .
#' @param sta is 
#' 
#' 
#' @export
Ph1MultipleTesting.Y0 <- function(model, nsim = 1000, FAP0 = 0.2, log = FALSE, const = 1, sta = FALSE, lowerbound = 0, side = "right-sided") {
  
  q <- dim(model$Phi)[1]
  n <- length(model$Y.tr) - q
  m <- dim(model$Phi)[2]
  
  tmpY.sim <- matrix(NA, nrow = n, ncol = nsim)
  tmpY.resi <- matrix(NA, nrow = n, ncol = nsim)
  
  Y.resi <- rep(NA, n)
  
  for (j in 1:nsim) {
    
    tmpsel <- sample(1:m, 1)
    Y0fit <- fit.GibbsRFLSM(model$Y.tr, 
                            model$Phi[, tmpsel], model$muq[tmpsel], 
                            X = model$X, Beta = model$Beta[, tmpsel], 
                            Kappa = model$Kappa[, tmpsel], 
                            H = NULL, Gamma = NULL, Tau = NULL)
    
    Y0.sim <- rnorm(n, Y0fit, sqrt(model$sigma2[tmpsel]))
    
    Y0.sim <- backtrans(Y0.sim, log, const, sta, model$meanY, model$sdY)
    
    Y0.sim[Y0.sim < lowerbound] <- lowerbound
    
    tmpY.sim[, j] <- Y0.sim
    
  }
  
  tmpY.sim.median <- rep(NA, n)
  
  for (i in 1:n) {
    tmpY.sim.median[i] <- median(tmpY.sim[i, ], na.rm = TRUE)
    tmpY.resi[i, ] <- tmpY.sim[i, ] - tmpY.sim.median[i]
    Y.resi[i] <- model$Y.ma[i + q] - tmpY.sim.median[i]
  }
  
  for (j in 1:nsim) {
    tmpY.resi[, j] <- tmpY.resi[, j] / sqrt(mean(tmpY.resi[, j] ^ 2))
      #(tmpY.resi[, j] - mean(tmpY.resi[, j])) / sd(tmpY.resi[, j])
  }
  
  Y.resi <- Y.resi / sqrt(mean(Y.resi ^ 2))
    #(Y.resi - mean(Y.resi)) / sd(Y.resi)
  
  
  tmpY.resi.max <- rep(NA, nsim)
  tmpY.resi.min <- rep(NA, nsim)
  
  for (j in 1:nsim) {
    tmpY.resi.max[j] <- max(tmpY.resi[, j], na.rm = TRUE)
    tmpY.resi.min[j] <- min(tmpY.resi[, j], na.rm = TRUE)
  }
  
  p.value.right <- mean(tmpY.resi.max >= max(Y.resi), na.rm = TRUE)
  p.value.left <- mean(tmpY.resi.min <= min(Y.resi), na.rm = TRUE)
  p.value.two <- mean((tmpY.resi.min <= min(Y.resi)) | (max(Y.resi) <= tmpY.resi.max), na.rm = TRUE)
    
  lower <- NA
  upper <- NA
  
  if (side == "right-sided") {
    p.value <- p.value.right
    lower <- NA
    upper <- quantile(tmpY.resi.max, 1 - FAP0)
    sig.ind <- Y.resi >= upper
  } else if (side == "left-sided") {
    p.value <- p.value.left
    lower <- quantile(tmpY.resi.min, FAP0)
    upper <- NA
    sig.ind <- Y.resi <= lower
  } else {
    p.value <- p.value.two
    lower <- quantile(tmpY.resi.min, FAP0 / 2)
    upper <- quantile(tmpY.resi.max, 1 - FAP0 / 2)
    sig.ind <- (Y.resi <= lower) | (Y.resi >= upper)
  }
  
  sig <- p.value <= FAP0
  
  if (side == "right-sided")
  
  list(sig = sig,
       p.value = p.value,
       sig.ind = sig.ind,
       cs = Y.resi,
       lower = lower,
       upper = upper,
       Y0.median = tmpY.sim.median,
       p.value.right = p.value.right, 
       p.value.left = p.value.left,
       p.value.two = p.value.two)
  
}


#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param model is model.
#' @param nsim is .
#' @param FAP0 is 
#' @param log is model.
#' @param const is .
#' @param sta is 
#' 
#' 
#' @export
Ph1MultipleTesting.resi <- function(model, interval = c(1e-8, 0.5 - 1e-8), log = FALSE, const = 1, sta = FALSE, 
                                    meanY = NULL, sdY = NULL, 
                                    nsim = 10000, FAP0 = 0.2, side = "right-sided", tol = 1e-8) {
  
  root.finding <- function(adj.alpha, sim0, FAP0 = 0.2, side = "right-sided") {
    m <- dim(sim0)[1]
    nsim <- dim(sim0)[2]
    
    tmplower <- rep(-Inf, m)
    tmpupper <- rep(Inf, m)
    
    sig <- matrix(NA, nrow = m, ncol = nsim)
    
    for (i in 1:m) {
      if (side == "right-sided") {
        tmpupper[i] <- quantile(sim0[i, ], 1 - adj.alpha, na.rm = TRUE)
      } else if (side == "left-sided") {
        tmplower[i] <- quantile(sim0[i, ], adj.alpha, na.rm = TRUE)
      } else if (side == "two-sided") {
        tmplower[i] <- quantile(sim0[i, ], adj.alpha/2, na.rm = TRUE)
        tmpupper[i] <- quantile(sim0[i, ], 1 - adj.alpha/2, na.rm = TRUE)
      }
    }
    
    for (j in 1:nsim) {
      sig[, j] <- (tmplower <= sim0[, j]) & (sim0[, j] <= tmpupper)
    }
    
    tmp <- 1 - mean(colSums(sig) == m, na.rm = TRUE)
    
    #cat("tmp:", tmp, "and adj.alpha:", adj.alpha, "\n")
    
    tmp - FAP0
  }
  
  nnsim <- dim(model$Phi)[2]
  q <- dim(model$Phi)[1]
  n <- length(model$Y.tr)
  
  resi0 <- matrix(NA, n - q, nsim)
  resi1 <- matrix(NA, n - q, nsim)
  Y0.ma <- matrix(NA, n - q, nsim)
  
  max.resi0 <- rep(NA, nsim)
  max.resi1 <- rep(NA, nsim)
  
  for (i in 1:nsim) {
    tmpsel <- sample(1:nnsim, 1)
    tmpfit <- fit.GibbsRFLSM(model$Y.tr, Phi = model$Phi[, tmpsel], muq = model$muq[tmpsel],
                          X = model$X, Beta = model$Beta[, tmpsel], Kappa = model$Kappa[, tmpsel], 
                          H = NULL, Gamma = NULL, Tau = NULL) 
    
    tmpresi0 <- model$Y.tr[-c(1:q)] - tmpfit
    tmpsigma2 <- mean(tmpresi0 ^ 2, na.rm = TRUE)
    
    tmpY0.ma <- rnorm(n - q, tmpfit, sqrt(tmpsigma2))
    tmpY0.ma <- backtrans(tmpY0.ma, log, const, sta, model$meanY, model$sdY)
    tmpY0.ma[tmpY0.ma < 0] <- 0
    Y0.ma[, i] <- tmpY0.ma

  }
  
  adj.alpha <- uniroot(root.finding, interval = interval, 
                       sim0 = Y0.ma, FAP0 = FAP0, side = side, tol = tol)$root
  
  lower <- rep(-Inf, n - q)
  upper <- rep(Inf, n - q)
  
  for (i in 1:(n - q)) {
    if (side == "right-sided") {
      upper[i] <- quantile(Y0.ma[i, ], 1 - adj.alpha, na.rm = TRUE)
    } else if (side == "left-sided") {
      lower[i] <- quantile(Y0.ma[i, ], adj.alpha, na.rm = TRUE)
    } else if (side == "two-sided") {
      lower[i] <- quantile(Y0.ma[i, ], adj.alpha/2, na.rm = TRUE)
      upper[i] <- quantile(Y0.ma[i, ], 1 - adj.alpha/2, na.rm = TRUE)
    }
  }
  
  cs <- model$Y.ma[-c(1:q)]
  
  sig <- 1 - ((lower <= cs) & (cs <= upper))
  
  list("cs" = cs, "sig" = sig, "lower" = lower, "upper" = upper, 
       "adj.alpha" = adj.alpha, "Y0.ma" = Y0.ma)
  
}


#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param model is model.
#' @param nsim is .
#' @param FAP0 is 
#' @param log is model.
#' @param const is .
#' @param sta is 
#' 
#' 
#' @export
Ph1MultipleTesting.Y0tr <- function(model, FAP0 = 0.2, side = "right-sided", nsim = 10000, interval = c(0.000001, 0.499999)) {
  
  root.finding <- function(adj.alpha, ph1mat, FAP0, n, nsim, side = "right-sided") {
    
    lim <- matrix(NA, nrow = n, ncol = 2)
    sig <- matrix(NA, nrow = n, ncol = nsim)
    
    for (i in 1:n) {
      if (side == "right-sided") {
        lim[i, 1] <- -Inf
        lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha)
      } else if (side == "left-sided") {
        lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha)
        lim[i, 2] <- infert
      } else if (side == "two-sided") {
        lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha / 2)
        lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha / 2)
      }
    }
    
    for (i in 1:nsim) {
      sig[, i] <- (lim[, 1] <= ph1mat[, i]) & (ph1mat[, i] <= lim[, 2])
    }
    
    tmp <- mean(colSums(sig) == n)
    dif <- tmp - (1 - FAP0)
    ##cat("dif:", dif, "\n")
    return(dif)
  }
  
  q <- dim(model$Phi)[1]
  n <- length(model$Y)
  nnsim <- dim(model$Phi)[2]
  
  ph1mat <- matrix(NA, nrow = n - q, ncol = nsim)
  
  for (i in 1:nsim) {
    tmpsel <- sample(1:nnsim, 1)
    Mu0 <- matrix(rep(model$muq[tmpsel], n))
    if (!is.null(model$X)) {
      Mu0 <- Mu0 + model$X %*% (model$Beta[, tmpsel] * model$Kappa[, tmpsel])
    }
    ph1mat[, i] <- simYph1(matrix(model$Yyj[, tmpsel]), matrix(model$Phi[, tmpsel]), matrix(Mu0), 
                            matrix(model$sigma2[tmpsel]), matrix(model$theta[tmpsel]), 1e-32)
  }
  
  adj.alpha <- uniroot(root.finding, interval, ph1mat = ph1mat, FAP0 = FAP0, n = n - q, nsim = nsim, side = side, 
          tol = 1e-6)$root
  
  lim <- matrix(NA, nrow = n - q, ncol = 2)
  sig <- matrix(NA, nrow = n - q, ncol = 1)
  
  for (i in 1:(n - q)) {
    if (side == "right-sided") {
      lim[i, 1] <- -Inf
      lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha, na.rm = TRUE)
    } else if (side == "left-sided") {
      lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha, na.rm = TRUE)
      lim[i, 2] <- infert
    } else if (side == "two-sided") {
      lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha / 2, na.rm = TRUE)
      lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha / 2, na.rm = TRUE)
    }
  }
  
  sig <- 1 - ((lim[, 1] <= model$Y[-c(1:q)]) & (model$Y[-c(1:q)] <= lim[, 2]))
  
  list("grandsig" = sum(sig) > 0, "sig" = sig, "lim" = lim, "adj.alpha" = adj.alpha, 
       "Yph1" = ph1mat)
  
}

#' Bayesian LASSO Phase I Monitoring
#' 
#' gets a posterior sample using Gibbs sampling for Random Flexible Level Shift Model
#' @param model is model.
#' @param nsim is .
#' @param FAP0 is 
#' @param log is model.
#' @param const is .
#' @param sta is 
#' 
#' 
#' @export
Ph1MultipleTesting.Y01 <- function(model, bset, 
                                   FAP0 = 0.2, side = "right-sided", 
                                   nsim = 10000, interval = c(0.000001, 0.499999)) {
  
  root.finding <- function(adj.alpha, ph1mat, FAP0, n, nsim, side = "right-sided") {
    
    lim <- matrix(NA, nrow = n, ncol = 2)
    sig <- matrix(NA, nrow = n, ncol = nsim)
    
    for (i in 1:n) {
      if (side == "right-sided") {
        lim[i, 1] <- -Inf
        lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha)
      } else if (side == "left-sided") {
        lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha)
        lim[i, 2] <- infert
      } else if (side == "two-sided") {
        lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha / 2)
        lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha / 2)
      }
    }
    
    for (i in 1:nsim) {
      sig[, i] <- (lim[, 1] <= ph1mat[, i]) & (ph1mat[, i] <= lim[, 2])
    }
    
    tmp <- mean(colSums(sig) == n)
    dif <- tmp - (1 - FAP0)
    ##cat("dif:", dif, "\n")
    return(dif)
  }
  
  q <- dim(model$Phi)[1]
  n <- length(model$Y)
  nnsim <- dim(model$Phi)[2]
  
  ph1mat <- matrix(NA, nrow = n - q, ncol = nsim)
  
  for (i in 1:nsim) {
    tmpsel <- sample(1:nnsim, 1)
    Mu0 <- matrix(rep(model$mu0[tmpsel], n))
    if (!is.null(model$X)) {
      Mu0 <- Mu0 + model$X %*% (model$Beta[, tmpsel])
    }
    ph1mat[, i] <- simYph1(matrix(model$Yyj[, tmpsel]), matrix(model$Phi[, tmpsel]), matrix(Mu0), 
                            matrix(model$sigma2[tmpsel]), matrix(model$theta[tmpsel]), 1e-32)
    
    tmpYyj <- yeojohnsontr(xx$Z[, tmpsel] + Y1, xx$theta[tmpsel], 1e-32)
    simYXph1(matrix(tmpYyj, ncol = 1), matrix(xx$Phi[, tmpsel], ncol = 1), matrix(tmpMu1[, tmpsel] + xx$mu0[tmpsel], ncol = 1),  
                   xx$sigma2[tmpsel],  xx$theta[tmpsel], 
                   1e-32, 1, 0, 1)
    
  }
  
  adj.alpha <- uniroot(root.finding, interval, ph1mat = ph1mat, FAP0 = FAP0, n = n - q, nsim = nsim, side = side, 
          tol = 1e-6)$root
  
  lim <- matrix(NA, nrow = n - q, ncol = 2)
  sig <- matrix(NA, nrow = n - q, ncol = 1)
  
  for (i in 1:(n - q)) {
    if (side == "right-sided") {
      lim[i, 1] <- -Inf
      lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha, na.rm = TRUE)
    } else if (side == "left-sided") {
      lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha, na.rm = TRUE)
      lim[i, 2] <- infert
    } else if (side == "two-sided") {
      lim[i, 1] <- quantile(ph1mat[i, ], adj.alpha / 2, na.rm = TRUE)
      lim[i, 2] <- quantile(ph1mat[i, ], 1 - adj.alpha / 2, na.rm = TRUE)
    }
  }
  
  sig <- 1 - ((lim[, 1] <= model$Y[-c(1:q)]) & (model$Y[-c(1:q)] <= lim[, 2]))
  
  list("grandsig" = sum(sig) > 0, "sig" = sig, "lim" = lim, "adj.alpha" = adj.alpha, 
       "Yph1" = ph1mat)
  
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
         main = "Gamma Diagnosis", 
         ylab = "Adjusted P-Value", 
         xlab = "")
    points((q + 1):(TT), adj.pvalue, type = 'o')
    points(((q + 1):(TT))[which(sig == TRUE)], adj.pvalue[which(sig == TRUE)], col = 'red', pch = 16)
    abline(h = FAP0, col = 'red')
    
    occpoint <- which(diff(H %*% sig) == 1) + 1 + q

    mulen <- dim(model$Mu)[1]
    mumedian <- rep(NA, mulen)
    
    for (i in 1:mulen) {
      mumedian[i] <- median(model$Mu[i, ], na.rm = TRUE)
    }
    
    plot(c(1, TT), c(min(mumedian), max(mumedian)), type = 'n',
         main = "Median Mu", 
         ylab = "Magnitude", 
         xlab = "")
    
    points(mumedian, type = 'o')
    points(occpoint, mumedian[occpoint], col = 'red', pch = 16)
    
    plot(c(1, TT), c(min(Y), max(Y)), type = 'n',
         main = "Y", 
         ylab = "Magnitude", 
         xlab = "")
    
    points(Y, type = 'o')
    points(occpoint, Y[occpoint], col = 'red', pch = 16)
    
  }
  
  out <- list("model" = model, "adj.pvalue" = adj.pvalue, 
              "sig" = sig) 
  out
} 

