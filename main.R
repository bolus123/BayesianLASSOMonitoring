
set.seed(1234)
n <- 100
p <- 5

Y <- arima.sim(list(ar = 0.5), n = n)
Y[1:50] <- Y[1:50] + 2

aa <- getPosteriorLayer2(Y, 5, 5, 1, 1, 1, 1, 
                         rep(1, 5) / 5,
                         100, 50, 1000, simlambda = TRUE, lambda = 5);

plot(Y)
points(aa$fit0[100, ], col = 'red')
points(aa$fit1[100, ], col = 'blue')

points(aa$beta00[1] + getSUMat(aa$U[, , 1]) %*% aa$Delta[1, ], 
       col = 'green')


aa <- getPosteriorLayer2(Y, 5, 5, 1, 1, 1, 1, 
                         rep(1, 5) / 5,
                         100, 50, 1000, simlambda = FALSE, lambda = 5);

bb <- getPosteriorLayer2(Y, 1, 0, 1, 1, 1, 1, 
                         rep(1, 5) / 5,
                         100, 50, 200, simlambda = FALSE, lambda = 1); 




ru <- aa$RU
bb <-getU(ru, aa$Delta, aa$sigma2, rep(1, 5) / 5)


X <- NULL

XX <- matrix(c(1, 2, 3, 4), ncol = 2)
getPosterior(Y, 1, 2, 3, 4, XX)


U <- matrix(c(1, 0, 0, 0, 1, 0, 1, 0, 0), ncol = 3, byrow = TRUE)

lambda2 <- 1000

f <- lm.fit(cbind(1, TT), Y)

oldbeta00 <- f$coefficients[1]
oldBeta00 <- getBeta00(beta00, n)
oldBeta20 <- f$coefficients[2:6]

oldsigma2 <- var(f$residuals)

oldtau200 <- getTau2(beta00, sigma2, lambda2)
oldTau220 <- getTau2Vec(Beta20, sigma2, lambda2)


newBeta00 <- getBeta00(0, n)
R00 <- getR(Y, 0, p, 0, 
            oldBeta00, 0, oldBeta20, 
            0, matrix(0, 1, 1), TT)
newbeta00 <- getbeta00(Y, R00, oldtau200, oldsigma2) 
newBeta00 <- getBeta00(newbeta00, n)


R20 <- getR(Y, 0, 0, 0, 
           newBeta00, 0, oldBeta20, 
           0, matrix(0, 1, 1), TT)


newBeta20 <- getBeta20Mono(Y, R20, TT, oldTau220, oldsigma2)

getR(arma::vec V, int q, int p, int K,
     arma::vec Beta00, arma::vec Beta10, arma::vec Beta20, 
     arma::vec UDelta, arma::mat X, arma::mat T) 