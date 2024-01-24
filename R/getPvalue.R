#' obtain the residual statistics 
#'
#' @param Y.hat is the transformed vector 
#' @param sigma2.hat is the transformed vector
#' @param Y.sim is the transformed vector
#' @param FAP0 is the matrix of laggy coefficients
#' @param side is the side
#' @param tol is the tolerance
#' @export
#' 
cc.ph1 <- function(Y.hat, sigma2.hat, Y.sim, FAP0 = 0.3, side = "two-sided", tol = 1e-6) {
  root.finding <- function(cc, FAP0, resi, nsim, TT, side) {
  
    sig <- rep(NA, nsim)
    for (i in 1:nsim) {
      if (side == "two-sided") {
        sig[i] <- sum((-cc <= resi[, i]) & (resi[, i] <= cc)) == TT
      } else if (side == "right-sided") {
        sig[i] <- sum((resi[, i] <= cc)) == TT
      } else if (side == "left-sided") {
        sig[i] <- sum((-cc <= resi[, i])) == TT
      }
      
    }
    
    tmpFAP0 <- 1 - mean(sig)
    d <- FAP0 - tmpFAP0
    
    #cat("cc:", cc, ", tmpFAP0:", tmpFAP0, "and d:", d, "\n")
    
    d
  }
  
  nsim <- dim(Y.sim)[2]
  TT <- dim(Y.sim)[1]
  
  resi <- matrix(NA, nrow = TT, ncol = nsim)
  for (i in 1:nsim) {
    resi[, i] <- (Y.sim[, i] - Y.hat)
    resi[, i] <- resi[, i] / sqrt(sigma2.hat)
  }
  
  interval <- c(0, max(abs(resi)))
  
  list("cc" = uniroot(root.finding, interval, FAP0 = FAP0, resi = resi, nsim = nsim, TT = TT, 
          side = side, tol = tol)$root, 
       "resi" = resi)
  
}  

#' obtain the residual statistics 
#'
#' @param Y.hat is the transformed vector 
#' @param sigma2.hat is the transformed vector
#' @param Y.sim is the transformed vector
#' @param FAP0 is the matrix of laggy coefficients
#' @param side is the side
#' @param tol is the tolerance
#' @export
#' 
adjalpha.ph1 <- function(Y.hat, sigma2.hat, Y.sim, FAP0 = 0.3, side = "two-sided", tol = 1e-6) {
  root.finding <- function(adjalpha, FAP0, resi, nsim, TT, side) {
    tmplower <- rep(NA, TT)
    tmpupper <- rep(NA, TT)
    
    for (i in 1:TT) {
      if (side == "two-sided") {
        tmplower[i] <- quantile(resi[i, ], adjalpha/2)
        tmpupper[i] <- quantile(resi[i, ], 1 - adjalpha/2)
      } else if (side == "right-sided") {
        tmplower[i] <- -Inf
        tmpupper[i] <- quantile(resi[i, ], 1 - adjalpha)
      } else if (side == "left-sided") {
        tmplower[i] <- quantile(resi[i, ], adjalpha)
        tmpupper[i] <- Inf
      }
      
    }
    
    sig <- rep(NA, nsim)
    for (i in 1:nsim) {
      sig[i] <- sum((tmplower <= resi[, i]) & (resi[, i] <= tmpupper)) == TT
    }
    
    tmpFAP0 <- 1 - mean(sig)
    d <- FAP0 - tmpFAP0
    
    #cat("cc:", cc, ", tmpFAP0:", tmpFAP0, "and d:", d, "\n")
    
    d
  }
  
  nsim <- dim(Y.sim)[2]
  TT <- dim(Y.sim)[1]
  
  resi <- matrix(NA, nrow = TT, ncol = nsim)
  for (i in 1:nsim) {
    resi[, i] <- (Y.sim[, i] - Y.hat)
    resi[, i] <- resi[, i] / sqrt(sigma2.hat)
  }
  
  interval <- c(0, 0.5)
  
  list("adjalpha" = uniroot(root.finding, interval, FAP0 = FAP0, resi = resi, nsim = nsim, TT = TT, 
                      side = side, tol = tol)$root, 
       "resi" = resi)
  
} 


#' obtain the residual statistics 
#'
#' @param Y.hat is the transformed vector 
#' @param sigma2.hat is the transformed vector
#' @param Y.sim is the transformed vector
#' @param ARL0 is the matrix of laggy coefficients
#' @param side is the side
#' @param tol is the tolerance
#' @export
#' 
cc.ph2 <- function(Y.hat, sigma2.hat, Y.sim, 
                   cs.mean, cs.sd, lambda = 0.05, ARL0 = 360, side = "two-sided", tol = 1e-6) {
  root.finding <- function(cc, ARL0, ewma, nsim, TT, side) {
 
    RL <- rep(NA, nsim)
    for (i in 1:nsim) {
      if (side == "two-sided") {
        tmp <- (-cc <= ewma[, i]) & (ewma[, i] <= cc)
      } else if (side == "right-sided") {
        tmp <- (ewma[, i] <= cc)
      } else if (side == "left-sided") {
        tmp <- (-cc <= ewma[, i])
      }
      tmpsig <- which(tmp == FALSE)
      if (length(tmpsig) > 0) {
        RL[i] <- min(tmpsig)
      } else {
        RL[i] <- TT
      }
    }
    
    tmpARL <- mean(RL)
    d <- ARL0 - tmpARL
    
    #cat("cc:", cc, ", tmpARL:", tmpARL, "and d:", d, "\n")
    
    d
  }
  
  nsim <- dim(Y.sim)[2]
  TT <- dim(Y.sim)[1]
  
  resi <- matrix(NA, nrow = TT, ncol = nsim)
  for (i in 1:nsim) {
    resi[, i] <- (Y.sim[, i] - Y.hat) / sqrt(sigma2.hat)
    resi[, i] <- (resi[, i] - cs.mean) / cs.sd
  }
  
  ewma.sim <- resi
  
  for (i in 2:TT) {
    ewma.sim[i, ] <- lambda * ewma.sim[i, ] + (1 - lambda) * ewma.sim[i - 1, ]
  }
  
  interval <- c(0, max(abs(resi)))
  
  list("cc" = uniroot(root.finding, interval, ARL0 = ARL0, ewma = ewma.sim, nsim = nsim, TT = TT, 
                      side = side, tol = tol)$root, 
       "ewma" = ewma.sim)
  
}

#' obtain the residual statistics 
#'
#' @param Y.hat is the transformed vector 
#' @param sigma2.hat is the transformed vector
#' @param Y.sim is the transformed vector
#' @param ARL0 is the matrix of laggy coefficients
#' @param side is the side
#' @param tol is the tolerance
#' @export
#' 
adjalpha.ph2 <- function(Y.hat, sigma2.hat, Y.sim, cs.mean, cs.sd, 
                         ARL0 = 360, side = "two-sided", tol = 1e-6) {
  root.finding <- function(adjalpha, ARL0, ewma, nsim, TT, side) {
    tmplower <- rep(NA, TT)
    tmpupper <- rep(NA, TT)
    
    tmplower <- rep(NA, TT)
    tmpupper <- rep(NA, TT)
    
    for (i in 1:TT) {
      if (side == "two-sided") {
        tmplower[i] <- quantile(ewma[i, ], adjalpha/2)
        tmpupper[i] <- quantile(ewma[i, ], 1 - adjalpha/2)
      } else if (side == "right-sided") {
        tmplower[i] <- -Inf
        tmpupper[i] <- quantile(ewma[i, ], 1 - adjalpha)
      } else if (side == "left-sided") {
        tmplower[i] <- quantile(ewma[i, ], adjalpha)
        tmpupper[i] <- Inf
      }
      
    }
    
    RL <- rep(NA, nsim)
    for (i in 1:nsim) {
      tmp <- (tmplower <= ewma[, i]) & (ewma[, i] <= tmpupper)
      tmpsig <- which(tmp == FALSE)
      if (length(tmpsig) > 0) {
        RL[i] <- min(tmpsig)
      } else {
        RL[i] <- TT
      }
    }
    
    tmpARL <- mean(RL)
    d <- ARL0 - tmpARL
    
    #cat("cc:", cc, ", tmpARL:", tmpARL, "and d:", d, "\n")
    
    d
  }
  
  nsim <- dim(Y.sim)[2]
  TT <- dim(Y.sim)[1]
  
  resi <- matrix(NA, nrow = TT, ncol = nsim)
  for (i in 1:nsim) {
    resi[, i] <- (Y.sim[, i] - Y.hat) / sqrt(sigma2.hat)
    resi[, i] <- (resi[, i] - cs.mean) / cs.sd
  }
  
  ewma.sim <- resi
  
  for (i in 2:TT) {
    ewma.sim[i, ] <- lambda * ewma.sim[i, ] + (1 - lambda) * ewma.sim[i - 1, ]
  }
  
  interval <- c(0, 0.5)
  
  list("adjalpha" = uniroot(root.finding, interval, ARL0 = ARL0, ewma = ewma.sim, nsim = nsim, TT = TT, 
                      side = side, tol = tol)$root, 
       "ewma" = ewma.sim)
  
}



#' obtain the residual statistics 
#'
#' @param Y.hat is the transformed vector 
#' @param sigma2.hat is the transformed vector
#' @param Y.sim is the transformed vector
#' @param FAP0 is the matrix of laggy coefficients
#' @param side is the side
#' @param tol is the tolerance
#' @export
#' 
lim.ph1 <- function(Y.hat, Y1.hat, sigma2.hat, Y.sim, FAP0 = 0.3, side = "two-sided") {
  
  nsim <- dim(Y.sim)[2]
  TT <- dim(Y.sim)[1]
  
  llr <- matrix(NA, nrow = TT, ncol = nsim)
  llr.max <- rep(NA, nsim)
  for (i in 1:nsim) {
    llr[, i] <- 2 * (dnorm(Y.sim[, i], Y1.hat, sqrt(sigma2.hat), log = TRUE) - 
           dnorm(Y.sim[, i], Y.hat, sqrt(sigma2.hat), log = TRUE))
    if (side == "right-sided") {
      llr[, i] <- llr[, i] * (Y1.hat > Y.hat)
    } else if (side == "left-sided") {
      llr[, i] <- llr[, i] * (Y1.hat < Y.hat)
    }
    llr.max[i] <- max(llr[, i])
  }
  
  lim <- quantile(llr.max, 1 - FAP0)
  
  lim
  
} 
