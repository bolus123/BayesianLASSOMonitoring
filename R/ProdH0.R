#' obtain the posterior probablity under H0
#' 
#' @param Tau is the distributions of Tau
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
#' ProdH0(result$Tau)
ProdH0 <- function(Tau) {
  
  mean(colSums(Tau) == 0)
  
}
