% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/others.R
\name{rzinpoisinar3}
\alias{rzinpoisinar3}
\title{simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift}
\usage{
rzinpoisinar3(n, alpha, lambda, pi, h, delta, burnin = 100)
}
\arguments{
\item{n}{is the length}

\item{alpha}{is the alpha}

\item{lambda}{is the mean of poisson mixture}

\item{pi}{is the proportion of zeros}

\item{h}{is the start point of shift}

\item{delta}{is the value of the standardized shift}

\item{burnin}{is the length of the burn-in period}
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
