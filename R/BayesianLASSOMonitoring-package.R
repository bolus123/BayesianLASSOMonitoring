#' @keywords package
#' @aliases BayesianLASSOMonitoring-package
"_PACKAGE"

#' Package BayesianLASSOMonitoring
#'
#' Package description.
#'
#' @name BayesianLASSOMonitoring-package
#' @references McCulloch, R. E., & Tsay, R. S. (1993). Bayesian inference and prediction for mean and variance shifts in autoregressive time series. Journal of the american Statistical association, 88(423), 968-978.
#' @import Rcpp VGAM ZINARp forecast
#' @importFrom stats dnorm
#' @importFrom stats rnorm
#' @importFrom stats arima
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom Rcpp evalCpp
#' @importFrom pracma sqrtm
#' @importFrom pracma Diag
#' @useDynLib BayesianLASSOMonitoring, .registration = TRUE
NULL

