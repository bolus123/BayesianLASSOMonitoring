library(GIGrvg)
library(breakfast)
library(pscl)
library(glmnet)

#################################################

simAR1 <- function(T, q, psi, sigma2, delta, tau) {
  tmp <- arima.sim(list(ar = c(psi)), T + q, sd = sqrt(sigma2))
  tmp[(tau + q):(T + q)] <- tmp[(tau + q):(T + q)] + 
    sqrt(sigma2 / (1 - psi ^ 2)) * delta
  V <- getT(tmp, q)
  V <- V[-c(1:q), ]
  Y <- tmp[-c(1:q)]
  list(Y = Y, V = V)
}

#################################################

getModel <- function(Y, V, shift = c("Isolated", "Sustained"), 
                     lambda2 = 5, burnin = 50, nsim = 1000, ntry = 10) {
  T <- length(Y)
  
  XShift <- vector()
  if (sum(shift %in% "Isolated") > 0) {
    XShift <- cbind(XShift, IsolatedShift(T))
  }
  
  if (sum(shift %in% "Sustained") > 0) {
    XShift <- cbind(XShift, SustainedShift(T))
  }
  
  q <- dim(V)[2]
  p <- dim(XShift)[2]
  
  betamat1 <- matrix(NA, nrow = nsim, ncol = 1 + q + p)
  ###############################################
  X <- cbind(V, XShift)
  m0 <- glmnet::glmnet(X, Y, lambda = lambda2, intercept = TRUE)
  
  beta0 <- m0$a0
  beta <- as.vector(m0$beta)
  beta1 <- beta[1:q]
  beta2 <- beta[(q + 1):(q + p)]
  ###############################################
  i <- 0
  flg <- 0
  while(i < ntry && flg == 0) {
    i <- i + 1
    m1 <- try(getPosterior(Y, V, XShift, lambda2, 
                           beta0, beta1, beta2, 
                           burnin = burnin, nsim = nsim), TRUE)
    if (class(m1) == "try-error") {
      flg <- 0
    } else {
      flg <- 1
    }
  }
  
  ###############################################
  m1
}

#################################################

getChart <- function(FAP0, Y, V, m1, 
                     shift = c("Isolated", "Sustained"), side = "one-sided",
                     DivType = "Residual", 
                     nsim = 1000,
                     interval = c(0, 10), tol = 1e-6) {
  T <- length(Y)
  
  XShift <- vector()
  if (sum(shift %in% "Isolated") > 0) {
    XShift <- cbind(XShift, IsolatedShift(T))
  }
  
  if (sum(shift %in% "Sustained") > 0) {
    XShift <- cbind(XShift, SustainedShift(T))
  }
  
  q <- dim(V)[2]
  p <- dim(XShift)[2]
  
  ##################################################
  
  ref <- getRefDivergence(Y, V, m1$beta[, 1], 
                          m1$beta[, 2:(q + 1)], V[1, ], DivType, nsim)
  
  cc <- bisection(FAP0, ref, side, 
                  min(interval), max(interval), tol)
  
  cs <- getDivergenceCS(Y, V, XShift, 
                        m1$beta[, 1], 
                        m1$beta[, 2:(q + 1)], 
                        m1$beta[, (q + 2):(1 + q + p)],
                        V[1, ], DivType)
  
  ##################################################
  
  if (side == "one-sided") {
    sig <- cs > cc
  } else if (side == "two-sided") {
    sig <- (cs < -cc) || (cc > cs) 
  }
  ##################################################
  list(cs = cs, cc = cc, sig = sig)
}

#################################################



#################################################