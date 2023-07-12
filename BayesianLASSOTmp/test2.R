q <- 5

IS <- IsolatedShift(n - q)
SS <- SustainedShift(n - q)

XX <- cbind(IS, SS)

n <- length(Y)

YY <- Y[(q + 1):n]

V <- getV(YY, q)

loglik(YY, V, XX, aa1$beta1[1, ], aa1$beta2[1, ], aa1$sigma2[1])

#debug(updatebeta2mcmc)
test2 <- updatebeta2mcmc(YY, V, XX, aa1$beta1[1, ], aa1$beta2[1, ], aa1$sigma2[1], aa1$invTau22[1, ], burnin = 50, ntry = 1)

debug(updatebeta1mcmc)
set.seed(12345)
test1 <- updatebeta1mcmc(YY, V, XX, aa1$beta1[1, ], aa1$beta2[1, ], aa1$sigma2[1], aa1$invTau21[1, ], burnin = 50, ntry = 1)

#debug(gibbsBLassomcmc)

debug(updatebeta1mcmc)
set.seed(12345)
start <- Sys.time()
test <- gibbsBLassomcmc(YY, V, XX, 
                lambda = NULL,
                updateLambda = TRUE,
                r = 1, 
                delta = 0.1,
                nsamp = 1000,
                burnin = 50,
                step = 1
)
end <- Sys.time()
end - start