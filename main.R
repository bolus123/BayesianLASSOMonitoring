nn <- 60

Y <- arima.sim(list(ar = 0.5), n = nn)
Ytilde <- Y - mean(Y)
X <- getX(Ytilde, nn, 5, 1, 1)

Beta <- getInitialBetaLS(Ytilde, nn, 5, 1, 1)
result <- optim(Beta, fn = getSSE, method = "BFGS", Y = Ytilde, X = X)
InitialBeta <- result$par
InitialSigma2 <- result$value / nn
InitialLambda <- sqrt(rgamma(1, shape = 0.1, rate = 0.1) * InitialSigma2)
aa <- getOmega(InitialBeta, InitialLambda)
#bb <- getBetaLASSO(Ytilde, X, aa)
#getLambdaEM(bb, InitialLambda, InitialSigma2, 
#            1000, 1e-6)

ddd <- getPosterior(Ytilde, X, 
                        aa,
                        InitialBeta, 
                        InitialLambda, 
                        InitialSigma2, 
                        1000, 
                        1000, 100, 1e-6) 

BetaNull <- getInitialBetaLSNull(Ytilde, nn, 5)
XNull <- getXNull(Ytilde, nn, 5)
resultNull <- optim(BetaNull, fn = getSSE, method = "BFGS", Y = Ytilde, 
                 X = XNull)
InitialBetaNull <- resultNull$par
InitialSigma2Null <- resultNull$value / nn
InitialLambdaNull <- sqrt(rgamma(1, shape = 0.1, rate = 0.1) * InitialSigma2Null)
InitialOmegaNull <- getOmega(InitialBetaNull, InitialLambdaNull)


dddNull <- getPosterior(Ytilde, XNull, 
                        InitialOmegaNull,
                        InitialBetaNull, 
                        InitialLambdaNull, 
                        InitialSigma2Null, 
                    1000, 
                    1000, 100, 1e-6) 

plot(Ytilde, type= 'l')
points(XNull %*% dddNull$Beta[2,], type = 'l', lty = 2)
ee <- forecast::auto.arima(Ytilde)
points(ee$fitted)