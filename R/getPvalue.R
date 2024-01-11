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
    tmplower <- rep(NA, TT)
    tmpupper <- rep(NA, TT)
    
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
    resi[, i] <- (Y.sim[, i] - Y.hat) / sqrt(sigma2.hat)
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
#' @param ARL0 is the matrix of laggy coefficients
#' @param side is the side
#' @param tol is the tolerance
#' @export
#' 
cc.ph2 <- function(Y.hat, sigma2.hat, Y.sim, ARL0 = 360, side = "two-sided", tol = 1e-6) {
  root.finding <- function(cc, ARL0, resi, nsim, TT, side) {
    tmplower <- rep(NA, TT)
    tmpupper <- rep(NA, TT)
    
    RL <- rep(NA, nsim)
    for (i in 1:nsim) {
      if (side == "two-sided") {
        tmp <- (-cc <= resi[, i]) & (resi[, i] <= cc)
      } else if (side == "right-sided") {
        tmp <- (resi[, i] <= cc)
      } else if (side == "left-sided") {
        tmp <- (-cc <= resi[, i])
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
  }
  
  interval <- c(0, max(abs(resi)))
  
  list("cc" = uniroot(root.finding, interval, ARL0 = ARL0, resi = resi, nsim = nsim, TT = TT, 
                      side = side, tol = tol)$root, 
       "resi" = resi)
  
}  