library(GIGrvg)
library(breakfast)
library(pscl)


dat <-read.csv(file = "/Users/yuihuiyao/Library/CloudStorage/Box-Box/Yuhui R21/Walker County De-Identified 2016-2021 Opioid ER Visits.csv")

DateTime <- strptime(dat$Admit.Date.Time, format = "%m/%d/%Y %I:%M:%S %p", tz = "America/Chicago")

DateTime1 <- sort(DateTime)

dif <- diff(DateTime1) / 3600

#############

#which(dif == 0)

#218, 531, 555, 561#

#############

target <- 218

tmp <- dat[which(DateTime == DateTime1[target]), ]

if (tmp$Birth.Date.Time[1] == tmp$Birth.Date.Time[2]) {
  if (tmp$CRace.CEth.Combined.Broad[1] == tmp$CRace.CEth.Combined.Broad[2]) {
    if (tmp$Patient.Zip[1] == tmp$Patient.Zip[2]) {
      if (tmp$Hospital[1] == tmp$Hospital[2]) {
        dat <- dat[-which(DateTime == DateTime1[target]), ]
        dat <- rbind(dat, tmp[1, ])
      }
    }
  }
}

#############

target <- 531

tmp <- dat[which(DateTime == DateTime1[target]), ]

if (tmp$Birth.Date.Time[1] == tmp$Birth.Date.Time[2]) {
  if (tmp$CRace.CEth.Combined.Broad[1] == tmp$CRace.CEth.Combined.Broad[2]) {
    if (tmp$Patient.Zip[1] == tmp$Patient.Zip[2]) {
      if (tmp$Hospital[1] == tmp$Hospital[2]) {
        dat <- dat[-which(DateTime == DateTime1[target]), ]
        dat <- rbind(dat, tmp[1, ])
      }
    }
  }
}

#############

target <- 555

tmp <- dat[which(DateTime == DateTime1[target]), ]

if (tmp$Birth.Date.Time[1] == tmp$Birth.Date.Time[2]) {
  if (tmp$CRace.CEth.Combined.Broad[1] == tmp$CRace.CEth.Combined.Broad[2]) {
    if (tmp$Patient.Zip[1] == tmp$Patient.Zip[2]) {
      if (tmp$Hospital[1] == tmp$Hospital[2]) {
        dat <- dat[-which(DateTime == DateTime1[target]), ]
        dat <- rbind(dat, tmp[1, ])
      }
    }
  }
}

#############

target <- 561

tmp <- dat[which(DateTime == DateTime1[target]), ]

if (tmp$Birth.Date.Time[1] == tmp$Birth.Date.Time[2]) {
  if (tmp$CRace.CEth.Combined.Broad[1] == tmp$CRace.CEth.Combined.Broad[2]) {
    if (tmp$Patient.Zip[1] == tmp$Patient.Zip[2]) {
      if (tmp$Hospital[1] == tmp$Hospital[2]) {
          dat <- dat[-which(DateTime == DateTime1[target]), ]
          dat <- rbind(dat, tmp[1, ])
      }
    }
  }
}

#############

DateTime <- strptime(dat$Admit.Date.Time, format = "%m/%d/%Y %I:%M:%S %p", tz = "America/Chicago")

DateTime1 <- sort(DateTime)

dif <- diff(DateTime1) / 3600

ts <- data.frame(
  "DateTime" = DateTime1[-1], 
  "Waiting" = as.numeric(dif))

#############

tcubicroot <- ts$Waiting ^ (1/3)
tlog <- log(ts$Waiting + 0.5)

#############

w <- 5

tma <- rep(NA, dim(ts)[1])

for (i in 1:dim(ts)[1]) {
  
  if (i < w) {
    tma[i] <- mean(ts[1:i, 2])
  } else {
    tma[i] <- mean(ts[(i - w + 1):i, 2])
  }

}

#############

check <- which("2017-01-01" <= ts$DateTime & ts$DateTime <= "2020-12-31")

tcubicroot <- tcubicroot[check]
tlog <- tlog[check]
tma <- tma[check]

#############

T <- length(tma)

p <- 5

lambda20 <- 5.0
lambda21 <- 5.0

burnin <- 50
nsim <- 100

##########################

X1 <- IsolatedShift(length(tma))
X2 <- SustainedShift(length(tma))
#X3 <- GradualShift(length(tma))
#XShift <- cbind(X1, X2, X3)
XShift <- cbind(X1, X2)
                
q <- dim(XShift)[2]

betamat1 <- matrix(NA, nrow = nsim, ncol = 1 + q + p)
taumat1 <- matrix(NA, nrow = nsim, ncol = 1 + q + p)
sigmamat1 <- matrix(NA, nrow = nsim, ncol = 1)
psimat1 <- matrix(NA, nrow = nsim, ncol = 1)
fit0mat1 <- matrix(NA, nrow = nsim, ncol = T)
fit1mat1 <- matrix(NA, nrow = nsim, ncol = T)
lambda2mat1 <- matrix(NA, nrow = nsim, ncol = 1)

betamat2 <- matrix(NA, nrow = nsim, ncol = 1 + p)
taumat2 <- matrix(NA, nrow = nsim, ncol = 1 + p)
sigmamat2 <- matrix(NA, nrow = nsim, ncol = 1)
psimat2 <- matrix(NA, nrow = nsim, ncol = 1)
fit0mat2 <- matrix(NA, nrow = nsim, ncol = T)
fit1mat2 <- matrix(NA, nrow = nsim, ncol = T)
lambda2mat2 <- matrix(NA, nrow = nsim, ncol = 1)

##########################

Y <- tma;
V <- getT(tma, p)

X <- cbind(V, XShift)

m0 <- glmnet::glmnet(X, Y, lambda = lambda20, intercept = TRUE)

beta00 <- unlist(as.numeric(m0$a0[1]))
beta02 <- unlist(m0$beta[1:(p)])
beta01 <- unlist(m0$beta[(p + 1):dim(X)[2]])

fit1 <- as.vector(beta00 + X %*% m0$beta)

sigma02 <- var(Y - fit1)

beta <- c(beta00, beta02, beta01)

tau0 <- getTau2(beta, sigma02, lambda20)

taubeta00 <- tau0[1]
taubeta02 <- tau0[2:(p + 1)]
taubeta01 <- tau0[(p + 2):(dim(X)[2] + 1)]

#############

kk <- 0
for (i in 1:(nsim + burnin)) {
  tmpmodel <- getGaussianPosteriorCM(matrix(Y, ncol = 1), beta0 = beta00, tau2beta0 = taubeta00, sigma2 = sigma02, 
                               lambda2 = lambda20, 
                               X = XShift, V = V, beta1 = beta01, beta2 = beta02, 
                               tau2beta1 = taubeta01, tau2beta2 = taubeta02, q = q, p = p)
  sigma02 <- var(Y - tmpmodel$fit1)
  beta00 <- tmpmodel$betadelta[1]
  taubeta00 <- tmpmodel$tau2all[1]
  
  beta01 <- tmpmodel$betadelta[2:(dim(XShift)[2] + 1)]
  taubeta01 <- tmpmodel$tau2all[2:(dim(XShift)[2] + 1)]
    
  beta02 <- tmpmodel$betadelta[(dim(XShift)[2] + 2):(dim(XShift)[2] + 1 + p)]
  taubeta02 <- tmpmodel$tau2all[(dim(XShift)[2] + 2):(dim(XShift)[2] + 1 + p)]
  
  lambda20 <- getLambda2EM(tmpmodel$expectedtau2all)
  
  if (i > burnin) {
    kk <- kk + 1
    
    betamat1[kk, ] <- tmpmodel$betadelta
    taumat1[kk, ] <- tmpmodel$tau2all
    sigmamat1[kk, ] <- sigma02
    fit0mat1[kk, ] <- tmpmodel$fit0
    fit1mat1[kk, ] <- tmpmodel$fit1
 
    lambda2mat1[kk, ] <- lambda20
    
  }
}

##########################

alpha <- 0.5
CI <- matrix(NA, nrow = 1 + q + p, ncol = 3)

for (i in 1:(1 + q + p)) {
  
  CI[i, 1:2] <- quantile(betamat1[, i], c(alpha / 2, 1 - alpha / 2)) 
  CI[i, 3] <- 1 - (CI[i, 1] <= 0 & 0 <= CI[i, 2])
  
}

which(CI[, 3] == 1)

truebetamat1 <- colMeans(betamat1)
truebetamat1[which(CI[, 3] == 0)] <- 0

plot(tma)
points(cbind(1, XShift) %*% colMeans(betamat1[, 1:1167]), type = 'l', col = 'red')

plot(truebetamat1[2:(length(tma)+1)], col = 'red', type = 'l')
plot(X2 %*% (truebetamat1[(length(tma) + 2):(2 * length(tma)-1)]), col = 'red', type = 'l')
#plot(X3 %*% truebetamat1[(2 * length(tma)):(3 * length(tma)-3)], col = 'red', type = 'l')
plot(X2 %*% (truebetamat1[(length(tma) + 2):(2 * length(tma)-1)]) + X3 %*% truebetamat1[(2 * length(tma)):(3 * length(tma)-3)], 
     col = 'red', type = 'l')


plot(XShift[, 1:length(tma)] %*% colMeans(betamat1[, 1:length(tma)])


##########################

q <- 0

Y <- tma;
V <- getT(tma, p)

X <- cbind(V)

m0 <- glmnet::glmnet(X, Y, lambda = lambda21, intercept = TRUE)

beta00 <- as.numeric(m0$a0[1])
beta02 <- m0$beta[1:(p)]

fit1 <- as.vector(beta00 + X %*% m0$beta)

sigma02 <- var(Y - fit1)

beta <- matrix(c(beta00, beta02), ncol = 1)

tau0 <- getTau2(beta, sigma02, lambda21)

taubeta00 <- tau0[1]
taubeta02 <- tau0[2:(p + 1)]

kk <- 0
for (i in 1:(nsim + burnin)) {
  tmpmodel <- getGaussianPosteriorCM(matrix(Y, ncol = 1), beta0 = beta00, tau2beta0 = taubeta00, sigma2 = sigma02, 
                                     lambda2 = lambda21, 
                                     X = XShift, V = V, beta1 = beta01, beta2 = beta02, 
                                     tau2beta1 = taubeta01, tau2beta2 = taubeta02, q = q, p = p)
  sigma02 <- var(Y - tmpmodel$fit1)
  beta00 <- tmpmodel$betadelta[1]
  taubeta00 <- tmpmodel$tau2all[1]
  
  beta02 <- tmpmodel$betadelta[(2):(1 + p)]
  taubeta02 <- tmpmodel$tau2all[(2):(1 + p)]
  
  lambda21 <- getLambda2EM(tmpmodel$expectedtau2all)
  
  if (i > burnin) {
    kk <- kk + 1
    
    betamat2[kk, ] <- tmpmodel$betadelta
    taumat2[kk, ] <- tmpmodel$tau2all
    sigmamat2[kk, ] <- sigma02
    fit0mat2[kk, ] <- tmpmodel$fit0
    fit1mat2[kk, ] <- tmpmodel$fit1
    
    lambda2mat2[kk, ] <- lambda21
    
  }
}

#############



#############

nsim1 <- 10000

statmat <- matrix(NA, nrow = 10000, ncol = 1)

for (i in 1:nsim1) {

  target1 <- sample(1:nsim, 1)
  fit11 <- fit1mat1[target1, ]
  sigma21 <- sigmamat1[target1]

  target2 <- sample(1:nsim, 1)
  fit21 <- fit1mat2[target2, ]
  sigma22 <- sigmamat2[target2]
  
  statmat[i, 1] <- ((sum((Y - fit21)^2)) - (sum((Y - fit11)^2) )/ sigma21)
  
}