% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/others.R
\name{RMSE.ph1}
\alias{RMSE.ph1}
\title{simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift}
\usage{
RMSE.ph1(
  Ph1BayesianLASSO.model,
  log = TRUE,
  const = 1,
  sta = TRUE,
  lowerbound = 0
)
}
\arguments{
\item{Ph1BayesianLASSO.model}{is the length}

\item{log}{is the log}

\item{const}{is the constant}

\item{sta}{is the sta}
}
\description{
simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
}
\examples{
nsim <- 100
burnin <- 100
T <- 100
q <- 5
H <- getHMatMT(T, q)
Y <- arima.sim(list(ar = 0.5), n = T)

alpha <- c(0.03083069, 0.06242601, 0.09120189)
lambda <- 0.239385
pi <- 0.1453097

TT <- 183
w <- 28
Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)

}
