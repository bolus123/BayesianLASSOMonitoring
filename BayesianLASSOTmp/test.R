# library(pscl) # rigamma
#library(mnormt) # rmnorm
library(VGAM) # rinv.gaussian
#library(miscTools) # colMeans
library(mvtnorm)
library(truncnorm)
library(pracma)
library(BayesianLassoMonitoring)
library(parallel)

getV <- function(Y, q) {
  n <- length(Y)
  out <- matrix(0, nrow = n, ncol = q)
  for (i in 2:n) {
    if (i <= q) {
      out[i, 1:(i - 1)] <- Y[(i -1):1] 
    } else {
      out[i, ] <- Y[(i - 1):(i - q)]
    }
  }
  out
}


IsolatedShift <- function(n) {
  eye(n)
}

SustainedShift <- function(n) {
  out <- matrix(1, n, n)
  out <- tril(out)
  out[, c(-1, -n)]
}


loglik <- function(Y, V, X, beta1, beta2, sigma2) {
  n <- length(Y)
  q <- dim(V)[2]
  residuals <- rep(NA, n)
  
  for (i in (q+1):n) {
    Xq <- X[(i-q):(i-1), ]
    residuals[i] <- Y[i] - (X[i, ] %*% beta2 + V[i, ] %*% beta1 - t(beta2) %*% t(Xq) %*% beta1)
  }
  
  loglik <- - n / 2 * log(2 * pi * sigma2) - 1 / 2 / sigma2 * sum(residuals ^ 2, na.rm = TRUE)
  loglik
}

diflogfbetai <- function(newbetai, oldbetai, invTaui, sigma2) {
  -1 / 2 / sigma2 * invTaui * (newbetai ^ 2 - oldbetai ^ 2)
}

updatebeta2mcmc <- function(Y, V, X, beta1, beta2, sigma2, invTau22, burnin = 50, ntry = 1) {
  p <- length(beta2)
  oldbeta2vec <- beta2
  
  for (i in 1:p) {
    #cat("i:", i, "\n")
    oldbeta2 <- oldbeta2vec[i]
    invTau22i <- invTau22[i]
    oldll <- loglik(Y, V, X, beta1, oldbeta2vec, sigma2)
    newbeta2vec <- oldbeta2vec
    for (j in 1:(burnin + ntry)) {
      newbeta2 <- rnorm(1, mean = oldbeta2, sd = 1)
      newbeta2vec[i] <- newbeta2
      #cat("newbeta2:", newbeta2, ", invTau22i:", invTau22i, "\n")
      newll <- loglik(Y, V, X, beta1, newbeta2vec, sigma2)
      difb2i <- diflogfbetai(newbeta2, oldbeta2, invTau22i, sigma2)
      difll <- newll - oldll + difb2i
      ratio <- exp(difll)
      #cat("ratio:", ratio, "\n")
      alpha <- min(c(1, ratio))
      #cat("alpha:", alpha, "\n")
      coin <- rbinom(1, 1, alpha)
      if (coin == 1) {
        oldbeta2vec <- newbeta2vec
        oldbeta2 <- newbeta2
        if (j > burnin) {
          break 
        }
      }
      
    }
  }
  oldbeta2vec
}

rtwosegnorm <- function(n, a, b, mean = 0, sd = 1) {
  out <- rep(NA, n)
  U <- runif(n)
  
  a1 <- (-a - mean) / sd
  b1 <- (-b - mean) / sd
  a2 <- (a - mean) / sd
  b2 <- (b - mean) / sd
  
  p1 <- pnorm(a1) - pnorm(b1)
  p2 <- pnorm(b2) - pnorm(a2)
  p <- p1 + p2
  
  for (i in 1:n) {
    if (U[i] <= p1 / p) {
      out[i] <- qnorm(U[i] * p + pnorm(b1)) * sd + mean
    } else {
      out[i] <- qnorm(U[i] * p + pnorm(a2) - p1) * sd + mean
    }
  }
  out
}


updatebeta1mcmc <- function(Y, V, X, beta1, beta2, sigma2, invTau21, burnin = 50, ntry = 1, 
                            beta10 = Inf, beta1qplus1 = 0) {
  q <- length(beta1)
  oldbeta1vec <- beta1
  
  for (i in 1:q) {
    
    #cat("i:", i, "\n")
    oldbeta1 <- oldbeta1vec[i]
    invTau21i <- invTau21[i]
    oldll <- loglik(Y, V, X, oldbeta1vec, beta2, sigma2)
    newbeta1vec <- oldbeta1vec
    
    if (i > 1) {
      previousbeta <- abs(oldbeta1vec[i - 1])
    } else {
      previousbeta <- abs(beta10)
    }
    
    for (j in 1:(burnin + ntry)) {
      newbeta1 <- rtruncnorm(1, a = -previousbeta, b = previousbeta, 
                             mean = oldbeta1, sd = 1)
      newbeta1vec[i] <- newbeta1
      #cat("newbeta2:", newbeta2, ", invTau22i:", invTau22i, "\n")
      newll <- loglik(Y, V, X, newbeta1vec, beta2, sigma2)
      difb1i <- diflogfbetai(newbeta1, oldbeta1, invTau21i, sigma2)
      difll <- newll - oldll + difb1i
      ratio <- exp(difll) * 
        exp(log(dtruncnorm(oldbeta1, a = -previousbeta, b = previousbeta, mean = newbeta1, sd = 1)) - 
              log(dtruncnorm(newbeta1, a = -previousbeta, b = previousbeta, mean = oldbeta1, sd = 1)))
      #cat("ratio:", ratio, "\n")
      alpha <- min(c(1, ratio))
      #cat("alpha:", alpha, "\n")
      coin <- rbinom(1, 1, alpha)
      if (coin == 1) {
        oldbeta1vec <- newbeta1vec
        oldbeta1 <- newbeta1
        if (j > burnin) {
          break 
        }
      }
      
    }
  }
  oldbeta1vec
}


 
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
    beta2 <- drop(rmvnorm(1, mean2, varcov2, checkSymmetry = FALSE))
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
    invTau21 = invTau21Samples[select, ],
    invTau22 = invTau22Samples[select, ],
    sigma2 = sigma2Samples[select],
    lambda = lambdaSamples[select]
  )
  
  return(out)
}


gibbsBLassomcmc <- function(Y, V, X, 
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
  XY <- t(X) %*% Y
  beta2 <- drop(backsolve(XtX + diag(nrow=p), XY))
  
  rX <- drop(Y - X %*% beta2)
  
  VtV <- t(V) %*% V
  VrX <- t(V) %*% Y
  beta1 <- drop(backsolve(VtV + diag(nrow=q), VrX)) 
  
  beta <- c(beta1, beta2)
  
  residue <- drop(rX - V %*% beta1)
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
    beta1 <- updatebeta1mcmc(Y, V, X, beta1, beta2, sigma2, invTau21, burnin = 50, ntry = 1)
    
    beta1Samples[k,] <- beta1
    
    # sample beta2
    beta2 <- updatebeta2mcmc(Y, V, X, beta1, beta2, sigma2, invTau22, burnin = 50, ntry = 1)

    beta2Samples[k,] <- beta2
    
    beta <- c(beta1, beta2)
    
    
    # sample sigma2
    shape <- (n+m-1)/2
    rX <- drop(Y - X %*% beta2)
    residue <- drop(rX - V %*% beta1)
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
    invTau21 = invTau21Samples[select, ],
    invTau22 = invTau22Samples[select, ],
    sigma2 = sigma2Samples[select],
    lambda = lambdaSamples[select]
  )
  
  return(out)
}


######################################

getSim <- function(X, pars, q, seed, nsamp = 1000, burnin = 50, step = 1) {
  
  set.seed(seed + X)
  
  fap0 <- pars[X, 1]
  n <- pars[X, 2]
  phi <- pars[X, 3]
  delta <- pars[X, 4]
  
  ############################
  
  Y <- arima.sim(list(ar = phi), n = n)
  Y[(round(n/2) + 1):n] <- Y[(round(n/2) + 1):n] + delta * sqrt(1 / (1 - phi ^ 2))
  
  ############################
  
  intercept <- mean(Y)
  Y1 <- Y - intercept
  
  ############################
  
  V <- getV(Y1, q)
  
  Y1 <- Y1[-c(1:q)]
  V <- V[-c(1:q), ]
  
  
  IS <- IsolatedShift(n - q)
  SS <- SustainedShift(n - q)
  
  XX <- cbind(IS, SS)
  
  aa1 <- gibbsBLasso(Y1, V, XX, 
                     lambda = NULL,
                     r = 1, 
                     delta = 0.1,
                     nsamp = nsamp,
                     burnin = burnin,
                     step = step) 
  
  #bb1 <- getPvalue(aa1$beta2, "two-sided")
  
  cc1 <- BenjaminiHochberg(fap0, aa1$beta2, "two-sided", 
                           1, 0.0001)
  
  dd1 <- BonferroniCorrection(fap0, aa1$beta2, "two-sided", 
                              1, 0.0001)
  
  
  ############################
  
  c(sum(cc1[, 4]) == 0, sum(dd1[, 4]) == 0)
  
  
}


######################################

nsim <- 1000

seed <- 12345
q <- 5

fap0vec <- c(0.2)
phiVec <- c(0.5)
#tVec <- c(100, 200, 300)
tVec <- c(100)
deltaVec <- c(0, 0.1)

pars <- expand.grid(fap0vec, tVec, phiVec, deltaVec)

######################################

cl <- parallel::makeCluster(2)
parallel::clusterEvalQ(cl, require(VGAM))
parallel::clusterEvalQ(cl, require(mvtnorm))
parallel::clusterEvalQ(cl, require(truncnorm))
parallel::clusterEvalQ(cl, require(pracma))
parallel::clusterEvalQ(cl, require(BayesianLassoMonitoring))
parallel::clusterExport(cl, c("rmvnormMono", "gibbsBLasso", "getSim"))

parallel::parLapplyLB(cl = cl, X = 1:nsim, fun = getSim, pars = pars, q = q, 
                      seed = seed, nsamp = 1000, burnin = 50, step = 1)

