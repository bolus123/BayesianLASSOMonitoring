% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/others.R
\name{rarma}
\alias{rarma}
\title{simulate realizations using ARMA(p, q) and one sustained shift}
\usage{
rarma(
  object,
  n,
  h,
  delta,
  xreg = NULL,
  nsim = 100,
  burnin = 50,
  lowerbound = 0
)
}
\arguments{
\item{n}{is the length}

\item{h}{is the proportion of zeros}

\item{delta}{is the start point of shift}

\item{burnin}{is the length of the burn-in period}

\item{phi}{is the alpha}

\item{theta}{is the mean of poisson mixture}

\item{sigma2}{is the mean of poisson mixture}
}
\description{
simulate realizations using ARMA(p, q) and one sustained shift
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
