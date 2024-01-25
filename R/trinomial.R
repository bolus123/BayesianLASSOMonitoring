#' Get a pmf of trinomial distribution
#' 
#' gets the simulated data
#' @param x is x
#' @param n is trials
#' @param p0 is ties
#' 
#' @export
dtrinomial <- function(x, n, p0) {
  nd <- x
  nd <- abs(nd)
  kk <- 0:((n - nd) / 2)
  tmp <- rep(length(kk))
  for (k in kk) {
    x <- c(nd + k, k, n - nd - 2 * k)
    size <- n
    prob <- c((1 - p0) / 2, (1 - p0) / 2, p0)
    tmp[k + 1] <- dmultinom(x, size, prob)
  }
  sum(tmp)
}

dtrinomial <- Vectorize(dtrinomial, "x")

#' Get a cdf of trinomial distribution
#' 
#' gets the simulated data
#' @param x is x
#' @param n is trials
#' @param p0 is ties
#' 
#' @export
ptrinomial <- function(x, n, p0) {
  sum(dtrinomial((-n):x, n, p0))
}

#' Get a quantile of trinomial distribution
#' 
#' gets the simulated data
#' @param p is p
#' @param n is trials
#' @param p0 is ties
#' 
#' @export
qtrinomial <- function(p, n, p0) {
  
  rootfinding <- function(x, n, p0, p) {
    nd <- round(x)
    dif <- ptrinomial(nd, n, p0) - p
    #cat("nd:", nd, "and dif:", dif, "\n")
    dif
  }
  
  if (p > 0.5) {
    interval <- c(0, n)
  } else {
    interval <- c(-n, 0)
  }

  out <- round(uniroot(rootfinding, interval, n = n, p0 = p0, p = p)$root)
  if (ptrinomial(out, n, p0) > p) {
    out <- out - 1
    if (out < (-n)) {
      out <- -n
    }
  }
  out
}