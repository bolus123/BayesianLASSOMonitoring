#' do the multiple testing for all shifts using Benjamini-Hochberg procedure
#' 
#' @param TauGamma is the distributions of shifts.
#' @param FDR0 is a given constant for FDR.
#' @param tail is type of the test.
#' @param method do the procedure with or without kernel smoothing.
#' @examples
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' result <- GibbsRFLSM(Y, q, diag(nrow = q), 0.1, 0.1, 0.1, 0.1, 1, 1, 0.1, "MonoALASSO", 0, Inf, 1000, 1, 100, 1e-10, H)
#'
#' BHProcRFLSM(result$Tau * result$Gamma)
BHProcRFLSM <- function(TauGamma, FDR0 = 0.1, tail = "2-sided", method = "raw") {
  
  if (method == "raw") {
    pvalue <- getPvalueRFLSM(TauGamma, tail)
  } else if (method == "ks") {
    pvalue <- getPvalueKSRFLSM(TauGamma, tail) 
  }
  
  nn <- length(pvalue)
  rr <- rank(pvalue)
  bound <- rr / nn * FDR0
  sig <- pvalue <= bound
  if (sum(sig == 1) > 0) {
    sig <- rr <= max(rr[sig == 1])
    out <- cbind(pvalue, rr, bound, sig)
  } else {
    out <- cbind(pvalue, rr, bound, 0)
  }
  return(out)
}