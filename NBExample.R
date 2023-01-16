source(file = 'C:/Users/yyao17/Documents/GitHub/BayesianLassoMonitoring/head.R')



tol <- 1e-6

nsim <- 1000
burnin1 <- 50
burnin2 <- 30

lambda2 <- 100
K <- 20
p <- 10

n <- 730

spl <- 0.7

V <- arima.sim(list(ar = 0.5), n = n)
V[100:190] <- V[100:190] + 3
#Y[1:50] <- Y[1:50] + 1
#Y[51:100] <- Y[51:100] + 5
#Y[101:150] <- Y[101:150] + 10
#Y[151:200] <- Y[151:200] + 20
#Y[201:250] <- Y[201:250] + 30

psi <- 1

Y <- rnbinom(n, psi, mu = exp(V))

#####################################################################


#debug(posteriorNB)
aa <- posteriorNB(Y = Y, X = NULL, W = NULL, p = p, K = K, 
                  lambda = 10, c1 = 1, c2 = 1, 
                  updatelambda = TRUE, updatepsi = TRUE, 
                  clustering = "kmeans", spl = spl, 
                  changepoint = "idetect", 
                  burnin1 = burnin1, burnin2 = burnin2, nsim = nsim, tol = 1e-6)




plot(Y)

for (i in 1:nsim) {
  points(aa$Y1[i, ], col = 'blue')
  points(aa$Y0[i, ], col = 'red')
}


plot(c(1, n), c(min(V, aa$Shift$V, aa$NoShift$V), max(V, aa$Shift$V, aa$NoShift$V)), type = 'n')

for (i in 1:nsim) {
  points(aa$Shift$V[i, ], col = 'blue')
  points(aa$NoShift$V[i, ], col = 'red')
}

points(V, type = 'l')
