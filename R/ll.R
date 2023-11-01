#' obtain the log likelihood
#' 
#' @param Y is a vector
#' @param Phi is the coefficient
#' @param Mu is the mean
#' @param sigma2 is the variance
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
#' ll(Y, result$Phi, result$Mu, result$sigma2)
ll <- function(Y, Phi, Mu, sigma2) {
  
  T <- length(Y)
  q <- dim(Phi)[1]
  nsim <- length(sigma2)
  ll <- matrix(NA, nrow = T - q, ncol = nsim)
  
  for (ii in seq(nsim)) {
    V <- matrix(Y, ncol = 1) - Mu[, ii]
    Vas <- getV(V, q)
    V <- V[-c(1:q)]
    Vas <- Vas[-c(1:q), ]
    resi <- V - Vas %*% Phi[, ii]
    ll[, ii] <- dnorm(resi, mean = 0, sd = sqrt(sigma2[ii]), log = TRUE)
  }
  
  ll
  
}
