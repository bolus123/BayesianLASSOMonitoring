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

V1 <- arima.sim(list(ar = 0.5), n = n)
V1[100:190] <- V1[100:190] - 3
#Y[1:50] <- Y[1:50] + 1
#Y[51:100] <- Y[51:100] + 5
#Y[101:150] <- Y[101:150] + 10
#Y[151:200] <- Y[151:200] + 20
#Y[201:250] <- Y[201:250] + 30

V2 <- arima.sim(list(ar = 0.5), n = n)
V2[200:290] <- V2[200:290] - 3
#Y[1:50] <- Y[1:50] + 1
#Y[51:100] <- Y[51:100] + 5
#Y[101:150] <- Y[101:150] + 10
#Y[151:200] <- Y[151:200] + 20
#Y[201:250] <- Y[201:250] + 30


psi <- 1
prob <- 1 - 1 / (1 + exp(-V2))

V2tmp <- rbinom(n, 1, prob)
Y <- V2tmp * rnbinom(n, psi, mu = exp(V1))



#####################################################################


#debug(posteriorNB)
aa <- posteriorZinfNB(Y = Y, X = NULL, W = NULL, p = p, K = K, 
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
