#' get the P value for RFLSM with kernel smoothing
#' 
#' @param TauGamma is the distributions of shifts.
#' @param tail is type of the test.
#' @export
#' @examples
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, q, diag(nrow = q), 0.1, 0.1, 0.1, 0.1, 
#' 1, 1, 0.1, "MonoALASSO", Inf, 0, 1000, 1, 100, 1e-10, H)
#'
#' getPvalueKSRFLSM(result$Tau * result$Gamma)
getPvalueKSRFLSM <- function(TauGamma, tail = "2-sided") {
  tmp <- TauGamma
  
  nn <- dim(tmp)[1]
  pvalue <- rep(0, nn)
  
  for (i in 1:nn) {
    
    tmpkde <- density(tmp[i, ])
    rtmpdens <- spatstat.explore::CDF(tmpkde)
    tmpp <- rtmpdens(0)
    
    if (tail == "2-sided") {
      pvalue[i] <- 2 * min(1 - tmpp, tmpp)
    } else if (tail == "left-sided") {
      pvalue[i] <- tmpp
    } else if (tail == "right-sided") {
      pvalue[i] <- 1 - tmpp
    }
    
  }
  return(pvalue)
}