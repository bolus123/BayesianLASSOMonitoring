# library(pscl) # rigamma
#library(mnormt) # rmnorm
library(VGAM) # rinv.gaussian
#library(miscTools) # colMeans
library(mvtnorm)
library(truncnorm)



rmvnormMono <- function(mean, varcov) {
  
  m <- length(mean)
  out <- rep(0, m)
  
  for (i in 1:m) {
    if (i > 1) {
      tmpinv <- solve(varcov[1:(i-1), 1:(i-1)])
      tmpvarcov <- varcov[1:(i-1), i]
      mu <- mean[i] + tmpvarcov %*% tmpinv %*% (out[1:(i-1)] - mean[1:(i-1)])
      sigma <- sqrt(varcov[i, i] - tmpvarcov %*% tmpinv %*% tmpvarcov)
      out[i] <- rtruncnorm(1, a = - abs(out[i - 1]), b = abs(out[i - 1]), mean = mu, sd = sigma)
    } else {
      mu <- mean[i]
      sigma <- sqrt(varcov[i, i])
      out[i] <- rnorm(1, mu, sigma)
    }
    
  }
  return(out)
} 


gibbsBLasso <- function(Y, V, X, 
                       lambda = NULL,
                       updateLambda = TRUE,
                       r = 1, 
                       delta = 0.1,
                       nsamp = 1000,
                       burnin = 100,
                       step = 5#,
                       #max.steps = 100000, 
                       #intercept = TRUE
                       ) {
  n <- length(Y)
  q <- ncol(V)
  p <- ncol(X)
  m <- p + q
  
  #if (intercept == TRUE) {
  #  intercept <- mean(Y)
  #  Y <- Y - intercept
  #} else {
  #  intercept <- 0
  #}
  
  XtX <- t(X) %*% X
  VtV <- t(V) %*% V
  
  VX <- cbind(V, X)
  VXtVX <- t(VX) %*% VX
  VXY <- t(VX) %*% Y
  
  beta <- drop(backsolve(VXtVX + diag(nrow=m), VXY))
  beta1 <- beta[1:q]
  beta2 <- beta[(q + 1):m]
  residue <- drop(Y - VX %*% beta)
  sigma2 <- drop((t(residue) %*% residue) / n)
  invTau2 <- 1 / (beta * beta)
  invTau21 <- invTau2[1:q]
  invTau22 <- invTau2[(q + 1):m]
  
  if (is.null(lambda)) {
    lambda <- m * sqrt(sigma2) / sum(abs(beta))
  }
  
  totSim <- burnin + nsamp * step
  
  
  beta1Samples <- matrix(0, totSim, q)
  beta2Samples <- matrix(0, totSim, p)
  sigma2Samples <- rep(0, totSim)
  invTau21Samples <- matrix(0, totSim, q)
  invTau22Samples <- matrix(0, totSim, p)
  lambdaSamples <- rep(0, totSim)
  
 
  
  k <- 0
  while (k < totSim) {
    k <- k + 1
    
    #if (k %% 100 == 0) {
    #  cat('Iteration:', k, "\r")
    #}
    
    # sample beta1
    Ytilda <- drop(Y - X %*% beta2)
    VY <- t(V) %*% Ytilda
    
    invD1 <- diag(invTau21)
    invA1 <- solve(VtV + invD1)
    mean1 <- invA1 %*% VY
    varcov1 <- sigma2 * invA1
    #beta <- drop(rmnorm(1, mean, varcov))
    #beta2 <- drop(rmvnorm(1, mean2, varcov2)) # this one needs to be constrainted
    beta1 <- rmvnormMono(mean1, varcov1)
    beta1Samples[k,] <- beta1
    
    # sample beta2
    Ytilda <- drop(Y - V %*% beta1)
    XY <- t(X) %*% Ytilda
    
    invD2 <- diag(invTau22)
    invA2 <- solve(XtX + invD2)
    mean2 <- invA2 %*% XY
    varcov2 <- sigma2 * invA2
    #beta <- drop(rmnorm(1, mean, varcov))
    beta2 <- drop(rmvnorm(1, mean2, varcov2))
    beta2Samples[k,] <- beta2
    
    beta <- c(beta1, beta2)
    
    
    # sample sigma2
    shape <- (n+m-1)/2
    residue <- drop(Y - VX %*% beta)
    scale <- (t(residue) %*% residue + t(beta) %*% diag(c(invTau2)) %*% beta)/2
    sigma2 <- 1/rgamma(1, shape, 1/scale)
    sigma2Samples[k] <- sigma2
    
    # sample tau2
    muPrime <- sqrt(lambda^2 * sigma2 / beta^2)
    lambdaPrime <- lambda^2
    invTau2 <- rep(0, m)
    for (i in seq(m)) {
      invTau2[i] <- rinv.gaussian(1, muPrime[i], lambdaPrime)
    }
    invTau21 <- invTau2[1:q]
    invTau22 <- invTau2[(q + 1):m]
    
    invTau21Samples[k, ] <- invTau21
    invTau22Samples[k, ] <- invTau22
    
    # update lambda
    if (updateLambda == TRUE) {
      shape = r + m/2
      scale = delta + sum(1/invTau2)/2
      lambda <- rgamma(1, shape, 1/scale)
    }
    # if (k %% 10 == 0) {
    # low <- k - 9
    # high <- k
    # lambda <- sqrt( 2*m / sum(colMeans(invTau2Samples[low:high, ])) )
    # }
    lambdaSamples[k] <- lambda
  }
  
  #colMedians(betaSamples[seq(max.steps/2, max.steps, 5), ])
  
  select <- seq(burnin + 1, totSim, step)
  
  out <- list(
    beta1 = beta1Samples[select, ],
    beta2 = beta2Samples[select, ],
    #invTau21 = invTau21Samples[select, ],
    #invTau22 = invTau22Samples[select, ],
    sigma2 = sigma2Samples[select],
    lambda = lambdaSamples[select]
  )
  
  return(out)
}

######################################

library(glmnet)
library(BayesianLassoMonitoring)

seed <- 12345

controlStat <- 0.2

q <- 5
T <- 300

#######################################
set.seed(seed)



#######################################

# get a simulated process
Y <- arima.sim(list(ar = 0.5), n = T)

#####
#Y[(round(T/2) + 1):T] <- Y[(round(T/2) + 1):T] + 2 * sqrt(1 / (1 - 0.5 ^ 2))

intercept <- mean(Y)
Y <- Y - intercept

# get lagged variables
V <- getV(Y, q)

# get shifts
IS <- IsolatedShift(T)
SS <- SustainedShift(T)

X <- cbind(IS, SS)

aa1 <- gibbsBLasso(Y, V, X, 
                       lambda = NULL,
                       r = 1, 
                       delta = 0.1,
                       nsamp = 10000,
                       burnin = 100,
                       step = 1) 

bb1 <- getPvalue(aa1$beta2, "two-sided")

cc1 <- BenjaminiHochberg(0.2, aa1$beta2, "two-sided", 
                  1, 0.0001)

dd1 <- BonferroniCorrection(0.2, aa1$beta2, "two-sided", 
                        1, 0.0001)

#######################################

# get a simulated process
Y <- arima.sim(list(ar = 0.5), n = T)

#####
#Y[(round(T/2) + 1):T] <- Y[(round(T/2) + 1):T] + 2 * sqrt(1 / (1 - 0.5 ^ 2))

intercept <- mean(Y)
Y <- Y - intercept

# get lagged variables
V <- getV(Y, q)

# get shifts
IS <- IsolatedShift(T)
SS <- SustainedShift(T)

X <- cbind(IS, SS)

aa2 <- gibbsBLasso(Y, V, X, 
                  lambda = NULL,
                  r = 1, 
                  delta = 0.1,
                  nsamp = 10000,
                  burnin = 100,
                  step = 1) 

bb2 <- getPvalue(aa2$beta2, "two-sided")

cc2 <- BenjaminiHochberg(0.2, aa2$beta2, "two-sided", 
                        1, 0.0001)

dd2 <- BonferroniCorrection(0.2, aa2$beta2, "two-sided", 
                           1, 0.0001)

#######################################

# get a simulated process
Y <- arima.sim(list(ar = 0.5), n = T)

#####
#Y[(round(T/2) + 1):T] <- Y[(round(T/2) + 1):T] + 2 * sqrt(1 / (1 - 0.5 ^ 2))

intercept <- mean(Y)
Y <- Y - intercept

# get lagged variables
V <- getV(Y, q)

# get shifts
IS <- IsolatedShift(T)
SS <- SustainedShift(T)

X <- cbind(IS, SS)

aa3 <- gibbsBLasso(Y, V, X, 
                  lambda = NULL,
                  r = 1, 
                  delta = 0.1,
                  nsamp = 10000,
                  burnin = 100,
                  step = 1) 

bb3 <- getPvalue(aa3$beta2, "two-sided")

cc3 <- BenjaminiHochberg(0.2, aa3$beta2, "two-sided", 
                        1, 0.0001)

dd3 <- BonferroniCorrection(0.2, aa3$beta2, "two-sided", 
                           1, 0.0001)
