source(file = 'C:/Users/yyao17/Documents/GitHub/BayesianLassoMonitoring/head.R')



tol <- 1e-6

nsim <- 1000
burnin <- 100

lambda2 <- 100
K <- 20
p <- 10

n <- 730

spl <- 0.7

Y <- arima.sim(list(ar = 0.5), n = n)
Y[100:190] <- Y[100:190] + 3
#Y[1:50] <- Y[1:50] + 1
#Y[51:100] <- Y[51:100] + 5
#Y[101:150] <- Y[101:150] + 10
#Y[151:200] <- Y[151:200] + 20
#Y[201:250] <- Y[201:250] + 30


aa <- posteriorGaussian(Y, X = NULL, W = NULL, p = p, K = K, lambda = 10,
                            updatelambda = TRUE, 
                            clustering = "kmeans", spl = spl, 
                            changepoint = "idetect", 
                            burnin = burnin, nsim = nsim, tol = 1e-6)



plot(Y)

for (i in 1:nsim) {
  points(aa$Y1[i, ], col = 'blue')
  points(aa$Y0[i, ], col = 'red')
}
