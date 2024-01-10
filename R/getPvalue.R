#' obtain the residual statistics 
#' 
#' @param Y.sim is the transformed vector
#' @param FAP0 is the matrix of laggy coefficients
#' @param interval is the interval
#' @param tol is the tolerance
#' @export
#' 
adjalp <- function(Y.sim, FAP0 = 0.3, side = "two-sided", interval = c(0.00001, 0.4), tol = 1e-6) {
  root.finding <- function(adjalp, FAP0, Y.sim, nsim, TT, side) {
    tmplower <- rep(NA, TT)
    tmpupper <- rep(NA, TT)
    
    for (i in 1:TT) {
      if (side == "two-sided") {
        tmplower[i] <- quantile(Y.sim[i, ], adjalp / 2)
        tmpupper[i] <- quantile(Y.sim[i, ], 1 - adjalp / 2)
      } else if (side == "right-sided") {
        tmplower[i] <- -Inf #quantile(Y.sim[i, ], adjalp)
        tmpupper[i] <- quantile(Y.sim[i, ], 1 - adjalp)
      } else if (side == "left-sided") {
        tmplower[i] <- quantile(Y.sim[i, ], adjalp)
        tmpupper[i] <- Inf #quantile(Y.sim[i, ], 1 - adjalp)
      }
      
    }
    
    tmp <- rep(NA, TT)
    for (i in 1:nsim) {
      tmp[i] <- sum((tmplower <= Y.sim[, i]) & (Y.sim[, i] <= tmpupper)) == TT
    }
    
    tmpFAP0 <- 1 - mean(tmp)
    d <- FAP0 - tmpFAP0
    
    #cat("adjalp:", adjalp, ", tmpFAP0:", tmpFAP0, "and d:", d, "\n")
    
    d
  }
  
  nsim <- dim(Y.sim)[2]
  TT <- dim(Y.sim)[1]
  
  uniroot(root.finding, interval, FAP0 = FAP0, Y.sim = Y.sim, nsim = nsim, TT = TT, 
          side = side, tol = tol)$root
  
}  