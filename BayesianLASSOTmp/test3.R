q <- 5

IS <- IsolatedShift(n - q)
SS <- SustainedShift(n - q)

XX <- cbind(IS, SS)

n <- length(Y)

fap0 <- pars[X, 1]
n <- pars[X, 2]
phi <- pars[X, 3]
delta <- pars[X, 4]

############################

Y <- arima.sim(list(ar = phi), n = n)
Y[(round(n/2) + 1):n] <- Y[(round(n/2) + 1):n] + delta * sqrt(1 / (1 - phi ^ 2))

############################

intercept <- mean(Y)
Y1 <- Y - intercept

############################

V <- getV(Y, q)

YY <- Y1[-c(1:q)]
V <- V[-c(1:q), ]




#######################################################

set.seed(12345)
start <- Sys.time()
test1 <- gibbsBLasso(YY, V, XX, 
                lambda = 100,
                updateLambda = FALSE,
                r = 1, 
                delta = 0.1,
                nsamp = 1000,
                burnin = 100,
                step = 5#,
                #max.steps = 100000, 
                #intercept = TRUE
)
end <- Sys.time()
end - start

set.seed(12345)
start <- Sys.time()
test2 <- gibbsBLasso(YY, V, XX, 
                     lambda = 200,
                     updateLambda = FALSE,
                     r = 1, 
                     delta = 0.1,
                     nsamp = 1000,
                     burnin = 100,
                     step = 5#,
                     #max.steps = 100000, 
                     #intercept = TRUE
)
end <- Sys.time()
end - start

set.seed(12345)
start <- Sys.time()
test3 <- gibbsBLasso(YY, V, XX, 
                     lambda = 400,
                     updateLambda = FALSE,
                     r = 1, 
                     delta = 0.1,
                     nsamp = 1000,
                     burnin = 100,
                     step = 5#,
                     #max.steps = 100000, 
                     #intercept = TRUE
)
end <- Sys.time()
end - start

set.seed(12345)
start <- Sys.time()
test4 <- gibbsBLasso(YY, V, XX,  
                     lambda = 800,
                     updateLambda = FALSE,
                     r = 1, 
                     delta = 0.1,
                     nsamp = 1000,
                     burnin = 100,
                     step = 5#,
                     #max.steps = 100000, 
                     #intercept = TRUE
)
end <- Sys.time()
end - start

set.seed(12345)
start <- Sys.time()
test5 <- gibbsBLasso(YY, V, XX, 
                     lambda = 1600,
                     updateLambda = FALSE,
                     r = 1, 
                     delta = 0.1,
                     nsamp = 1000,
                     burnin = 100,
                     step = 5#,
                     #max.steps = 100000, 
                     #intercept = TRUE
)
end <- Sys.time()
end - start



aa1 <- BayesianLassoMonitoring::BenjaminiHochberg(0.2, test5$beta2, "two-sided", 1, 0.0001
)
