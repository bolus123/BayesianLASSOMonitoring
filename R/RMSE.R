#' obtain the root squared error
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
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
#' RMSE(Y, result$Phi, result$Mu)
RMSE <- function(Y, Phi, Mu) {
  
  T <- length(Y)
  q <- dim(Phi)[1]
  nsim <- dim(Phi)[2]
  ee <- matrix(NA, nrow = T - q, ncol = nsim)
  rmse <- rep(0, nsim)
  
  for (ii in seq(nsim)) {
    V <- matrix(Y, ncol = 1) - Mu[, ii]
    Vas <- getV(V, q)
    V <- V[-c(1:q)]
    Vas <- Vas[-c(1:q), ]
    ee[, ii] <- V - Vas %*% Phi[, ii]
    rmse[ii] <- sqrt(sum(ee[, ii] ^ 2))
  }
  rmse
}