#' get the P value for RFLSM
#' 
#' @param TauGamma is the distributions of shifts.
#' @param tail is type of the test.
#' @examples
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, q, diag(nrow = q), 0.1, 0.1, 0.1, 0.1, 
#' 1, 1, 0.1, "MonoALASSO", 0, Inf, 1000, 1, 100, 1e-10, H)
#'
#' getPvalueRFLSM(result$Tau * result$Gamma)
getPvalueRFLSM <- function(TauGamma, tail = "2-sided") {
  tmp <- TauGamma
  
  nn <- dim(tmp)[1]
  pvalue <- rep(0, nn)
  
  for (i in 1:nn) {
    
    if (tail == "2-sided") {
      pvalue[i] <- 2 * min(1 - mean(tmp[i, ] <= 0), mean(tmp[i, ] <= 0))
    } else if (tail == "left-sided") {
      pvalue[i] <- mean(tmp[i, ] <= 0)
    } else if (tail == "right-sided") {
      pvalue[i] <- 1 - mean(tmp[i, ] <= 0)
    }
    
  }
  return(pvalue)
}