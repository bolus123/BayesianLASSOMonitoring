#' get a true proportion of H0
#' 
#' @param TauGamma is the distributions of shifts.
#' @param tail is type of the test.
#' @param method get p values with or without kernel smoothing.
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
#' getH0(result$Tau * result$Gamma)
getH0 <- function(TauGamma, tail = "2-sided", FDR.method = "meinshausen", pvalue.method = "raw") {
  
  pvalue <- getPvalue(TauGamma, tail, pvalue.method)
  
  get.pi0(pvalue, estim.method = FDR.method)
  
}