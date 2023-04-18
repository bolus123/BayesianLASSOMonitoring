library(GIGrvg)
library(breakfast)
library(pscl)
library(glmnet)
library(parallel)
Rcpp::sourceCpp(file = "C:/Users/yyao17/Documents/GitHub/BayesianLassoMonitoring/src/FinalSource.cpp")
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
    m0 <- try(getPosteriorH0(Y, V, lambda2, 
                             beta0, beta1, 
                             burnin = burnin, nsim = nsim), TRUE)
    if (class(m0) == "try-error") {
      flg <- 0
    } else {
      flg <- 1
    }
  }
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
  list(m0 = m0, m1 = m1)
}

#################################################

getChart <- function(FAP0, Y, V, m0, m1, 
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
  
  ref <- getRefDivergence(Y, V, m0$beta[, 1], 
                          m0$beta[, 2:(q + 1)], V[1, ], DivType, nsim)
  
  cc <- bisection(FAP0, ref, side, 
                  min(interval), max(interval), tol)
  
  cs <- getDivergenceCS(Y, V, XShift, 
                        m0$beta[, 1], 
                        m0$beta[, 2:(q + 1)],
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

wrap <- function(X, pars, tau, 
                 shift = c("Isolated", "Sustained"), 
                 lambda2 = 5, burnin = 50, nbeta = 1000, ntry = 10, 
                 side = "one-sided",
                 DivType = "Residual", 
                 nref = 1000,
                 interval = c(0, 10), tol = 1e-6, seed = 12345) {
  
  #Rcpp::sourceCpp(file = "C:/Users/yyao17/Documents/GitHub/BayesianLassoMonitoring/src/FinalSource.cpp")
  
  set.seed(seed + X)
  
  s <- dim(pars)[1]
  
  out <- rep(NA, s)
  
  for (i in 1:s) {
    FAP0 <- pars[i, 1]
    T <- pars[i, 2]
    q <- pars[i, 3]
    psi <- pars[i, 4]
    sigma2 <- pars[i, 5]
    delta <- pars[i, 6]
    
    realization <- simAR1(T, q, psi, sigma2, delta, tau)
    model <- getModel(realization$Y, realization$V, shift = shift, 
                   lambda2 = lambda2, burnin = burnin, nsim = nbeta, ntry = ntry)
    chart <- getChart(FAP0, realization$Y, realization$V, model$m0, model$m1, 
                      shift = shift, side = side,
                      DivType = DivType, 
                      nsim = nref,
                      interval = interval, tol = tol)
    out[i] <- sum(chart$sig) > 0
  }
  out
  
}



#parFAP <- function(cl, pars, tau, 
#                   shift = c("Isolated", "Sustained"), 
#                   lambda2 = 5, burnin = 50, nbeta = 1000, ntry = 10, 
#                   side = "one-sided",
#                   DivType = "Residual", 
#                   nref = 1000,
#                   interval = c(0, 10), tol = 1e-6,
#                   nsim = 1000, seed = 12345) {
#  
#  
#  outMat <- vector()
#  
#  
#  for (i in 1:nsim) {
#    out <- parLapplyLB(cl = cl, 1:nsim, fun = wrap, 
#                       pars = pars, tau = tau, 
#                       shift = shift, 
#                       lambda2 = lambda2, burnin = burnin, 
#                       nbeta = nbeta, ntry = ntry, 
#                       side = side,
#                       DivType = DivType, 
#                       nref = nref,
#                       interval = interval, tol = tol, seed = seed)
#    outMat <- cbind(outMat, unlist(out))
#  }
#  
#  
#  outMat
#  
#}



#################################################

FAP0vec <- c(0.1)
Tvec <- c(365)
qvec <- c(5)
psivec <- c(0, 0.5)
sigma2vec <- c(1)
deltavec <- c(0, 1, 2, 3)
pars <- expand.grid(FAP0vec, Tvec, qvec, psivec, sigma2vec, deltavec)


no <- 6
addr <- paste("C:/Users/yyao17/Documents/GitHub/BayesianLassoMonitoring/out", no, '.Rdat', sep = "")


out <- vector()

for (i in 1:100) {
  undebug(wrap)
  undebug(getModel)
  tmp <- wrap(i, pars, 183, 
              shift = c("Isolated", "Sustained"), 
              lambda2 = 5, burnin = 100, nbeta = 1000, ntry = 10, 
              side = "one-sided",
              DivType = "Residual", 
              nref = 1000,
              interval = c(0, 10), tol = 1e-6, seed = 12345 + no)
  out <- rbind(out, tmp)
  save(out, file = addr)
}

