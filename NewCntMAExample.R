library(GIGrvg)
library(breakfast)
library(pscl)
library(glmnet)

dat <-read.csv(file = "/Users/yuihuiyao/Library/CloudStorage/Box-Box/Yuhui R21/Walker County De-Identified 2016-2021 Opioid ER Visits.csv")

#dat <- read.csv(file = "C:/Users/Yuhui/Box/Yuhui R21/Walker County De-Identified 2016-2021 Opioid ER Visits.csv")

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

#check <- which("2017-01-01" <= ts$DateTime & ts$DateTime <= "2020-12-31")

aa <- vector()

ii <- as.Date("2016-01-01", tz = "America/Chicago")

while(ii <= as.Date("2020-12-31")) {
  aa <- c(aa, as.character(ii))
  ii <- ii + 1
}

dd <- as.Date(ts$DateTime, tz = "America/Chicago")
cnt <- table(dd)

n <- length(aa)

bb <- matrix(0, nrow = n, ncol = 2)
bb[, 1] <- aa

for (i in 1:n) {
  target <- which(names(cnt) %in% bb[i, 1])
  if (length(target) != 0) {
    bb[i, 2] <- cnt[target]
  }
}

dat1 <- bb

#############

w <- 30

cntma <- rep(NA, dim(dat1)[1])

for (i in 1:dim(dat1)[1]) {
  
  if (i < w) {
    cntma[i] <- mean(as.numeric(dat1[1:i, 2]))
  } else {
    cntma[i] <- mean(as.numeric(dat1[(i - w + 1):i, 2]))
  }
  
}

#############
q <- 5
V <- getT(cntma, q)

check <- which(as.Date("2018-01-01") <= as.Date(dat1[, 1]) & 
                 as.Date(dat1[, 1]) <= as.Date("2018-12-31"))

cntma <- cntma[check]
V <- V[check, ]

ddate <- dat1[check, 1]

#############

T <- length(cntma)



lambda20 <- 5.0
lambda21 <- 5.0

burnin <- 50
nsim <- 100


X1 <- IsolatedShift(length(cntma))
X2 <- SustainedShift(length(cntma))
X3 <- GradualShift(length(cntma))
#XShift <- cbind(X1, X2, X3)
XShift <- cbind(X1, X2)
                
p <- dim(XShift)[2]

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

########################

########################

Y <- cntma;

X <- cbind(V, XShift)

m0 <- glmnet_cpp(X, Y, 5)

aaa(m0)

m0 <- glmnet::glmnet(X, Y, lambda = lambda20, intercept = TRUE)

beta0 <- m0$a0
beta <- as.vector(m0$beta)
beta1 <- beta[1:q]
beta2 <- beta[(q + 1):(q + p)]

ee <- getPosterior(Y, V, XShift, lambda20, 
             beta0, beta1, beta2, 
             burnin = 50, nsim = 100) 

fit0 <- cbind(1, V) %*% colMeans(ee$beta[, 1:(q + 1)])
fit1 <- cbind(1, V, XShift) %*% colMeans(ee$beta)

plot(Y)
points(fit0, type = 'l', col = 'blue')
points(fit1, type = 'l', col = 'red')

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

##############################


kk <- 0
for (i in 1:(nsim + burnin)) {
  
  cat("i:", i, "\n")
  
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

#save(betamat1, file = "/Users/yuihuiyao/Library/CloudStorage/Box-Box/Yuhui R21/betamat.csv")
load(file = "C:/Users/Yuhui/Box/Yuhui R21/betamat.csv")
load(file = '/Users/yuihuiyao/Library/CloudStorage/Box-Box/Yuhui R21/betamat.csv')
##########################

alpha <- 0.05
CI <- matrix(NA, nrow = 1 + q + p, ncol = 3)

for (i in 1:(1 + q + p)) {
  
  CI[i, 1:2] <- quantile(betamat1[, i], c(alpha / 2, 1 - alpha / 2)) 
  CI[i, 3] <- 1 - (CI[i, 1] <= 0 & 0 <= CI[i, 2])
  
}

which(CI[, 3] == 1)

truebetamat1 <- colMeans(betamat1)
truebetamat1[which(CI[, 3] == 0)] <- 0

plot(cntma)
points(cbind(1, V) %*% colMeans(betamat1[, c(1, (dim(betamat1)[2] - p + 1):dim(betamat1)[2])]), type = 'l', col = 'blue')
points(cbind(1, XShift, V) %*% colMeans(betamat1), type = 'l', col = 'red')



check1 <- which(as.Date("2017-01-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2017-03-31"))
check2 <- which(as.Date("2017-04-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2017-06-30"))
check3 <- which(as.Date("2017-07-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2017-09-30"))
check4 <- which(as.Date("2017-10-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2017-12-31"))


check5 <- which(as.Date("2018-01-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2018-03-31"))
check6 <- which(as.Date("2018-04-01") <= as.Date(ddate) & as.Date(ddate) <= 
                 as.Date("2018-06-30"))
check7 <- which(as.Date("2018-07-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2018-09-30"))
check8 <- which(as.Date("2018-10-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2018-12-31"))



check9 <- which(as.Date("2019-01-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2019-03-31"))
check10 <- which(as.Date("2019-04-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2019-06-30"))
check11 <- which(as.Date("2019-07-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2019-09-30"))
check12 <- which(as.Date("2019-10-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2019-12-31"))


check13 <- which(as.Date("2020-01-01") <= as.Date(ddate) & as.Date(ddate) <= 
                  as.Date("2020-03-31"))
check14 <- which(as.Date("2020-04-01") <= as.Date(ddate) & as.Date(ddate) <= 
                   as.Date("2020-06-30"))
check15 <- which(as.Date("2020-07-01") <= as.Date(ddate) & as.Date(ddate) <= 
                   as.Date("2020-09-30"))
check16 <- which(as.Date("2020-10-01") <= as.Date(ddate) & as.Date(ddate) <= 
                   as.Date("2020-12-31"))


plot(truebetamat1[2:(length(cntma)+1)], col = 'red', type = 'l')
plot(X2 %*% (truebetamat1[(length(cntma) + 2):(2 * length(cntma)-1)]), col = 'red', type = 'l')
#plot(X2 %*% (truebetamat1[(length(cntma) + 2):(2 * length(cntma)-1)]) + X3 %*% truebetamat1[(2 * length(cntma)):(3 * length(cntma)-3)], 
#     col = 'red', type = 'l')

##########################



logpdf <- function(Y, fit) {
  sigma2 <- var(Y - fit)
  log(dnorm(Y, fit, sqrt(sigma2)))
}


getKS <- function(Y, fit0, fit1, nsim, nnsim) {
  nn <- length(Y)
  resi0 <- Y - fit0
  resi1 <- Y - fit1
  out <- matrix(NA, nrow = nnsim, ncol = nn)
  for (j in 1:nnsim) {
    for (i in 1:nn) {
      #cat("j:", j, "i:", i, '\n');
      tmpx <- sample(resi0[, i], nsim, replace = TRUE)
      tmpy <- sample(resi1[, i], nsim, replace = TRUE)
      out[j, i] <- ks.test(tmpx, tmpy)$statistic
    }
  }
  out
}

getKSlogpdf <- function(Y, fit0, fit1, nsim, nnsim) {
  nn <- length(Y)
  
  logpdf0 <- matrix(NA, nrow = nsim, ncol = nn)
  logpdf1 <- logpdf0
  
  for (i in 1:nsim) {
    logpdf0[i, ] <- logpdf(Y, fit0[i, ])
    logpdf1[i, ] <- logpdf(Y, fit1[i, ])
  }
  
  out <- matrix(NA, nrow = nnsim, ncol = nn)
  for (j in 1:nnsim) {
    for (i in 1:nn) {
      #cat("j:", j, "i:", i, '\n');
      tmpx <- sample(logpdf0[, i], nsim, replace = TRUE)
      tmpy <- sample(logpdf1[, i], nsim, replace = TRUE)
      out[j, i] <- ks.test(tmpx, tmpy)$statistic
    }
  }
  out
}


findroot <- function(FAP0, Y, fit0, nsim, nnsim, interval = c(0, 5), tol = 1e-4) {
  
  rootfinding <- function(cc, FAP0, ref) {
    check <- ref <= cc
    ProbNSE <- mean(rowMeans(check))
    diff <- 1 - FAP0 - ProbNSE
    cat("cc:", cc, "and diff:", diff, '\n')
    return(diff)
  }
  
  ref <- getKS(Y, fit0, fit0, nsim, nnsim)
  
  uniroot(rootfinding, interval = interval, FAP0 = FAP0, 
          ref = ref, tol = tol)$root
  
}

findrootlogpdf <- function(FAP0, Y, fit0, nsim, nnsim, interval = c(0, 5), tol = 1e-4) {
  
  rootfinding <- function(cc, FAP0, ref) {
    check <- ref <= cc
    ProbNSE <- mean(rowMeans(check))
    diff <- 1 - FAP0 - ProbNSE
    cat("cc:", cc, "and diff:", diff, '\n')
    return(diff)
  }
  
  ref <- getKSlogpdf(Y, fit0, fit0, nsim, nnsim)
  
  uniroot(rootfinding, interval = interval, FAP0 = FAP0, 
          ref = ref, tol = tol)$root
  
}


getKSCS <- function(Y, fit0, fit1, nsim) {
  nn <- length(Y)
  resi0 <- Y - fit0
  resi1 <- Y - fit1
  out <- rep(NA, ncol = nn)
    for (i in 1:nn) {
      #cat("j:", j, "i:", i, '\n');
      tmpx <- resi0[, i]
      tmpy <- resi1[, i]
      out[i] <- ks.test(tmpx, tmpy)$statistic
    }
  out
}


getKSlogpdfCS <- function(Y, fit0, fit1, nsim) {
  nn <- length(Y)
  
  logpdf0 <- matrix(NA, nrow = nsim, ncol = nn)
  logpdf1 <- logpdf0
  
  for (i in 1:nsim) {
    logpdf0[i, ] <- logpdf(Y, fit0[i, ])
    logpdf1[i, ] <- logpdf(Y, fit1[i, ])
  }
  
  out <- rep(NA, ncol = nn)
    for (i in 1:nn) {
      #cat("j:", j, "i:", i, '\n');
      tmpx <- logpdf0[, i]
      tmpy <- logpdf1[, i]
      out[i] <- ks.test(tmpx, tmpy)$statistic
    }
  out
}


getKSlogpdfCSPPP <- function(Y, fit0, nsim) {
  nn <- length(Y)
  
  logpdf0rep <- matrix(NA, nrow = nsim, ncol = nn)
  logpdf0 <- logpdf0rep
  
  for (i in 1:nsim) {
    Sigma2 <- var(Y - fit0[i, ])
    Yrep <- rnorm(nn, fit0[i, ], sqrt(Sigma2))
    logpdf0rep[i, ] <- logpdf(Yrep, fit0[i, ])
    logpdf0[i, ] <- logpdf(Y, fit0[i, ])
  }
  
  out <- rep(NA, ncol = nn)
  for (i in 1:nn) {
    #cat("j:", j, "i:", i, '\n');
    tmpx <- logpdf0rep[, i]
    tmpy <- logpdf0[, i]
    out[i] <- ks.test(tmpx, tmpy)$statistic
  }
  out
}


getKSResiCSPPP <- function(Y, fit0, nsim) {
  nn <- length(Y)
  
  resi0rep <- matrix(NA, nrow = nsim, ncol = nn)
  resi0 <- resi0rep
  
  for (i in 1:nsim) {
    Sigma2 <- var(Y - fit0[i, ])
    Yrep <- rnorm(nn, fit0[i, ], sqrt(Sigma2))
    resi0rep[i, ] <- Yrep - fit0[i, ]
    resi0[i, ] <- Y - fit0[i, ]
  }
  
  out <- rep(NA, ncol = nn)
  for (i in 1:nn) {
    #cat("j:", j, "i:", i, '\n');
    tmpx <- resi0rep[, i]
    tmpy <- resi0[, i]
    out[i] <- ks.test(tmpx, tmpy)$statistic
  }
  out
}


findrootKSResiCSPPP <- function(FAP0, Y, fit0, nsim, nnsim, interval = c(0, 1), tol = 1e-4) {
  
  rootfinding <- function(cc, FAP0, ref) {
    check <- ref <= cc
    ProbNSE <- mean(rowMeans(check))
    diff <- 1 - FAP0 - ProbNSE
    cat("cc:", cc, "and diff:", diff, '\n')
    return(diff)
  }
  
  ref <- matrix(NA, nrow = nnsim, ncol = length(Y))
  
  for (i in 1:nnsim) {
    ref[i, ] <- getKSResiCSPPP(Y, fit0, nsim) 
  }
  
  uniroot(rootfinding, interval = interval, FAP0 = FAP0, 
          ref = ref, tol = tol)$root
  
}


getKSResiCSPPPAlt <- function(Y, fit0, fit1, nsim) {
  nn <- length(Y)
  
  resi0 <- matrix(NA, nrow = nsim, ncol = nn)
  resi1 <- resi0
  
  for (i in 1:nsim) {
    resi0[i, ] <- Y - fit0[i, ]
    resi1[i, ] <- Y - fit1[i, ]
  }
  
  out <- rep(NA, ncol = nn)
  for (i in 1:nn) {
    #cat("j:", j, "i:", i, '\n');
    tmpx <- resi0[, i]
    tmpy <- resi1[, i]
    out[i] <- ks.test(tmpx, tmpy)$statistic
  }
  out
}


getBinProb <- function(p, q, breaks = 10) {
  
  start <- min(p, q)
  end <- max(p, q)
  step <- (end - start) / breaks
  intervals <- rep(NA, breaks + 1)
  out <- matrix(NA, nrow = 2, ncol = breaks)
  for (i in 1:(breaks + 1)) {
    intervals[i] <- start + (i - 1) * (step)
    if (1 < i & i < (breaks + 1)) {
      out[1, i - 1] <- sum(intervals[i - 1] <= p & p < intervals[i])
      out[2, i - 1] <- sum(intervals[i - 1] <= q & q < intervals[i])
    } else if (i == breaks + 1) {
      out[1, i - 1] <- sum(intervals[i - 1] <= p & p <= intervals[i])
      out[2, i - 1] <- sum(intervals[i - 1] <= q & q <= intervals[i])
    }
  }
  out[1, ] <- out[1, ] / sum(out[1, ])
  out[2, ] <- out[2, ] / sum(out[2, ])
  out
}


SimLinear <- function(Y, beta0, beta1, initial = rep(0, length(beta1))) {
  n <- length(Y)
  Yrep <- rep(NA, n)
  Yrep <- c(initial, Yrep)
  sigma2 <- var(Y - beta0)
  q <- length(beta1)
  for (i in (q + 1):(n + q)) {
    Yrep[i] <- beta0 + Yrep[(i - q):(i - 1)] %*% beta1 + rnorm(1, 0, sqrt(sigma2))
  }
  Yrep[-c(1:q)]
}


getKLResiCSPPP <- function(Y, V, beta0, beta1, nsim, breaks = 10) {
  nn <- length(Y)
  
  resi0rep <- matrix(NA, nrow = nsim, ncol = nn)
  resi0 <- resi0rep
  
  p <- dim(beta1)[2]
  #V <- getT(c(initial, Y), p)[-c(1:p), ]
    
  initial <- V[1, ]
  
  for (i in 1:nsim) {
    Yrep <- SimLinear(Y, beta0[i], beta1[i, ], initial) 
    Vrep <- getT(c(initial, Yrep), p)[-c(1:p), ]
    resi0rep[i, ] <- Yrep -  beta0[i] - Vrep %*% beta1[i, ]
    resi0[i, ] <- Y -  beta0[i] - V %*% beta1[i, ]
  }
  
  out <- rep(NA, ncol = nn)
  for (i in 1:nn) {
    #cat("j:", j, "i:", i, '\n');
    tmpx <- resi0rep[, i]
    tmpy <- resi0[, i]
    tmp <- getBinProb(tmpy, tmpx, breaks = breaks)
    out[i] <- philentropy::KL(tmp, unit = "log")
  }
  out
}

findrootKLResiCSPPP <- function(FAP0, Y, V, beta0, beta1, 
                                nsim, nnsim, breaks = 10, 
                                interval = c(0, 5), tol = 1e-4) {
  
  rootfinding <- function(cc, FAP0, ref, n) {
    check <- ref <= cc
    ProbNSE <- mean(rowSums(check) == n)
    diff <- 1 - FAP0 - ProbNSE
    cat("cc:", cc, "and diff:", diff, '\n')
    return(diff)
  }
  
  ref <- matrix(NA, nrow = nnsim, ncol = length(Y))
  
  for (i in 1:nnsim) {
    ref[i, ] <- getKLResiCSPPP(Y, V, beta0, beta1, nsim, breaks = breaks)
  }
  
  n <- length(Y)
  
  #debug(rootfinding)
  uniroot(rootfinding, interval = interval, FAP0 = FAP0, 
          ref = ref, n = n, tol = tol)$root
  
}

getKLResiCSPPPAlt <- function(Y, V, X, beta0, beta1, beta2, nsim, breaks = 10) {
  nn <- length(Y)
  
  resi0 <- matrix(NA, nrow = nsim, ncol = nn)
  resi1 <- resi0
  
  for (i in 1:nsim) {
    resi1[i, ] <- Y -  beta0[i] - V %*% beta1[i, ] - X %*% beta2[i, ]
    resi0[i, ] <- Y -  beta0[i] - V %*% beta1[i, ]
  }
  
  out <- rep(NA, ncol = nn)
  for (i in 1:nn) {
    #cat("j:", j, "i:", i, '\n');
    tmpx <- resi1[, i]
    tmpy <- resi0[, i]
    tmp <- getBinProb(tmpy, tmpx, breaks = breaks)
    out[i] <- philentropy::KL(tmp, unit = "log")
  }
  out
}

########################

beta0 <- betamat1[, 1]
beta1 <- betamat1[, (dim(betamat1)[2] - p + 1):dim(betamat1)[2]]
beta2 <- betamat1[, (2):(dim(betamat1)[2] - p)]


fit0 <- t(cbind(1, V) %*% t(betamat1[, c(1, (dim(betamat1)[2] - p + 1):dim(betamat1)[2])]))
fit1 <- t(cbind(1, XShift, V) %*% t(betamat1))


#undebug(findrootKLResiCSPPP)
#debug(getKLResiCSPPP)
set.seed(12345)
bb1 <- findrootKLResiCSPPP(0.1, Y, V, beta0, beta1, 
                          nsim, nnsim, breaks = 20, interval = c(0, 20)) 

set.seed(12345)
bb2 <- findrootKLResiCSPPP(0.2, Y, V, beta0, beta1, 
                          nsim, nnsim, breaks = 20, interval = c(0, 20)) 

set.seed(12345)
bb3 <- findrootKLResiCSPPP(0.3, Y, V, beta0, beta1, 
                          nsim, nnsim, breaks = 20, interval = c(0, 20)) 
cc <- getKLResiCSPPPAlt(Y, V, XShift, beta0, beta1, beta2, nsim, breaks = 20)




########################

getResiCS <- function(Y, V, beta01, nsim) {
  
  nn <- length(Y)
  
  X <- cbind(1, V)
  sigma02 <- var(Y - X %*% beta0)
  
  resi0Mat <- matrix(NA, nrow = nsim, ncol = nn)
  resi0repMat <- matrix(NA, nrow = nsim, ncol = nn)
  for (i in 1:nsim) {
    resi0 <- Y - fit0[i, ]
    resi0Mat[i, ] <- resi0
  }
  
  mean0vec <- colMeans(resi0Mat)
  sigma02vec <- rep(NA, nn)
  
  Yrep <- rep(NA, nn)
  
  for (i in 1:nn) {
    sigma02vec[i] <- var(resi0Mat[, i] - mean0vec[i])
  }
  
  for (i in 1:nsim) {
    Yrep <- rnorm(nn, fit0[i, ], sqrt(sigma02vec))
    resi0rep <- Yrep - fit0[i, ]
    resi0repMat[i, ] <- resi0rep
  }
  

  mean0repvec <- colMeans(resi0repMat)
  
  out <- rep(NA, ncol = nn)
  
  out <- nsim * (mean0vec - mean0repvec) ^ 2 / (sigma02vec)
  #out <- log(out)
  out
}











findrootResiCSPPP <- function(FAP0, Y, fit0, nsim, nnsim, interval = c(0, 5), tol = 1e-4) {
  
  rootfinding <- function(cc, FAP0, ref, T) {
    check <- ref <= cc
    NSE <- rowSums(check) == T
    ProbNSE <- mean(NSE)
    diff <- (1 - FAP0) - ProbNSE
    cat("cc:", cc, "and diff:", diff, '\n')
    return(diff)
  }
  
  ref <- matrix(NA, nrow = nnsim, ncol = length(Y))
  
  for (i in 1:nnsim) {
    ref[i, ] <- getResiCS(Y, fit0, nsim) 
  }
  
  T <- length(Y)
  
  #ref <- log(ref)
  #debug(rootfinding)
  uniroot(rootfinding, interval = interval, FAP0 = FAP0, 
          ref = ref, tol = tol, T = T)$root
  
}


getResiCSAlt <- function(Y, fit0, fit1, nsim) {
  nn <- length(Y)
  
  resi0Mat <- matrix(NA, nrow = nsim, ncol = nn)
  resi1Mat <- matrix(NA, nrow = nsim, ncol = nn)
  for (i in 1:nsim) {
    resi0 <- Y - fit0[i, ]
    resi0Mat[i, ] <- resi0
  }
  
  mean0vec <- colMeans(resi0Mat)
  sigma02vec <- rep(NA, nn)
  
  Yrep <- rep(NA, nn)
  
  for (i in 1:nn) {
    sigma02vec[i] <- var(resi0Mat[, i] - mean0vec[i])
  }
  
  for (i in 1:nsim) {
    resi1Mat[i, ] <- rnorm(nn, mean0vec, sqrt(sigma02vec))
  }
  
  
  mean1vec <- colMeans(resi1Mat)
  
  out <- rep(NA, ncol = nn)
  
  out <- nsim * (mean0vec - mean1vec) ^ 2 / (sigma02vec)
  #out <- log(out)
  out
}






getResiCSLinear <- function(Y, V, beta0, beta1, nsim, initial = rep(0, length(beta1))) {
  
  nn <- length(Y)
  
  resi0Mat <- matrix(NA, nrow = nsim, ncol = nn)
  resi0repMat <- matrix(NA, nrow = nsim, ncol = nn)
  
  p <- dim(beta1)[2]
  
  out <- matrix(NA, nrow = nsim, ncol = nn)
  
  for (i in 1:nsim) {
    resi0Mat[i, ] <- Y - beta0[i] - V %*% beta1[i, ]
    sigma2 <- var(resi0Mat[i, ])
    Yrep <- SimLinear(Y, beta0[i], beta1[i, ], initial = initial)
    Vrep <- getT(c(initial, Yrep), p)[-c(1:p), ]
    resi0repMat[i, ] <- Yrep - beta0[i] - Vrep %*% beta1[i, ]
    out[i, ] <- (resi0Mat[i, ] - resi0repMat[i, ]) ^ 2  / sigma2
  }
  
  #mu0 <- colMeans(resi0Mat)
  #mu0rep <- colMeans(resi0repMat)
  #sigma02vec <- rep(NA, nn)
  #
  #for (i in 1:nn) {
  #  sigma02vec[i] <- var(resi0Mat[, i])
  #}
  #
  #out <- (mu0 - mu0rep)^2 / sigma02vec
  #log(out)
  colMeans(out)
}

findrootResiCSPPPLinear <- function(FAP0, Y, V, beta0, beta1, 
                                    nsim, nnsim, initial = rep(0, dim(beta1)[2]), 
                                    interval = c(0.5, 1), tol = 1e-4) {
  
  rootfinding <- function(cc, FAP0, ref, T) {
    
    #sigma2vec <- rep(NA, T)
    #uqvec <- rep(NA, T)
    #nsim <- dim(ref)[1]
    check <- matrix(NA, nrow = nsim, ncol = T)
    
    #for (i in 1:T) {
    #  sigma2vec[i] <- var(ref[, i])
    #}
    
    #for (i in 1:nsim) {
    #  check[i, ] <- ref[i, ] & ref[i, ] <= cc
    #}
    check <- ref <= cc
    NSE <- rowSums(check) == T
    ProbNSE <- mean(NSE)
    diff <- (1 - FAP0) - ProbNSE
    cat("cc:", cc, "and diff:", diff, '\n')
    return(diff)
  }
  
  ref <- matrix(NA, nrow = nnsim, ncol = length(Y))
  
  for (i in 1:nnsim) {
    ref[i, ] <- getResiCSLinear(Y, V, beta0, beta1, 
                                nsim, initial = initial)
  }
  
  T <- length(Y)
  
  #ref <- log(ref)
  #debug(rootfinding)
  uniroot(rootfinding, interval = interval, FAP0 = FAP0, 
          ref = ref, tol = tol, T = T)$root
  
}

getResiCSLinearAlt <- function(Y, V, X, beta0, beta1, beta2, nsim, initial = rep(0, length(beta1))) {
  
  nn <- length(Y)
  
  resi0Mat <- matrix(NA, nrow = nsim, ncol = nn)
  resi1Mat <- matrix(NA, nrow = nsim, ncol = nn)
  
  out <- matrix(NA, nrow = nsim, ncol = nn)
  
  p <- dim(beta1)[2]
  
  for (i in 1:nsim) {
    resi0Mat[i, ] <- Y - beta0[i] - V %*% beta1[i, ]
    sigma2 <- var(resi0Mat[i, ])
    resi1Mat[i, ] <- Y - beta0[i] - V %*% beta1[i, ] - X %*% beta2[i, ]
    out[i, ] <- (resi0Mat[i, ] - resi1Mat[i, ]) ^ 2 / sigma2
  }
  
  colMeans(out)
}



#debug(getResiCSLinear)
ee <- getResiCSLinear(Y, V, beta0, beta1, nsim, initial = V[1, ])

#debug(findrootResiCSPPPLinear)
aa <- findrootResiCSPPPLinear(0.2, Y, V, beta0, beta1, 
                   nsim, nnsim = 100, initial = V[1, ], 
                   interval = c(0.1, 20), tol = 1e-4)

bb <- getResiCSLinearAlt(Y, V, XShift, beta0, beta1, beta2, nsim, initial = V[1, ])


plot(c(1, 365), c(min(bb), max(bb)), 
     ylab = 'Divergence', 
     xlab = "", type = 'n', main = "Walker County, 2018",  xaxt = "n")
points(bb, type = 'o')

axis(side = 1, at = atvec, labels = FALSE)
axis(side = 1, at = atv, tick = FALSE,
     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                "Sep", "Oct", "Nov", "Dec"))


abline(h = aa)
legend("topright", legend = "FAP0 = 0.2", lty = 1, title = 'Control Limit')
#######################

num <- which(bb > aa)

beta21 <- matrix(0, nrow = dim(beta2)[1], ncol = dim(beta2)[2])
beta21[, num] <- beta2[, num]
beta21[, 365 + num - 1] <- beta2[, 365 + num - 1]



fit0 <- beta0 + V %*% t(beta1)
fit0 <- t(fit0)
fit1 <- beta0 + V %*% t(beta1) + XShift %*% t(beta21)
fit1 <- t(fit1)

plot(cntma, type = 'o', main = "Walker County, 2018", 
     ylab = 'Moving-average Count with 30 Days of Sliding Windows', 
     xlab = "", xaxt = "n")

atvec <- rep(NA, 13)
atvec[1] <- 1
atvec[2] <- atvec[1] + 31
atvec[3] <- atvec[2] + 28
atvec[4] <- atvec[3] + 31
atvec[5] <- atvec[4] + 30
atvec[6] <- atvec[5] + 31
atvec[7] <- atvec[6] + 30
atvec[8] <- atvec[7] + 31
atvec[9] <- atvec[8] + 31
atvec[10] <- atvec[9] + 30
atvec[11] <- atvec[10] + 31
atvec[12] <- atvec[11] + 30
atvec[13] <- atvec[12] + 31

atv <- rep(NA, 12)

for (i in 2:13) {
  atv[i - 1] <- (atvec[i] + atvec[i - 1]) / 2
}

axis(side = 1, at = atvec, labels = FALSE)
axis(side = 1, at = atv, tick = FALSE,
     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                "Sep", "Oct", "Nov", "Dec"))

for (i in 1:nsim) {
  points(fit1[i, ], type = 'l', col = 'pink')
  points(fit0[i, ], type = 'l', col = 'skyblue')
}


points(colMeans(fit1), type = 'l', col = 'red')
points(colMeans(fit0), type = 'l', col = 'blue')

legend("topright", legend = c("IC", "OOC"), lty = c(1, 1), col = c("blue", "red"), 
       title = "Fitted Curve")

#######################

eee <- beta0 + t(XShift %*% t(beta21))


plot(c(1, 365), c(min(eee), max(eee)), type = 'n', main = "Walker County, 2018", 
     ylab = 'Trend in the Mean', 
     xlab = "", xaxt = "n")

for (i in 1:nsim) {
  points(as.vector(eee[i, ]), col = 'gray', type = "l")
}


points(colMeans(eee), col = 'black', type = "l")


axis(side = 1, at = atvec, labels = FALSE)
axis(side = 1, at = atv, tick = FALSE,
     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                "Sep", "Oct", "Nov", "Dec"))


#######################

fff <-t(XShift[, 1:365] %*% t(beta21)[1:365, ])


plot(c(1, 365), c(min(fff), max(fff)), type = 'n', main = "Walker County, 2018", 
     ylab = 'Isolated Shift in the Mean', 
     xlab = "", xaxt = "n")

for (i in 1:nsim) {
  points(as.vector(fff[i, ]), col = 'gray', type = "l")
}


points(colMeans(fff), col = 'black', type = "l")


axis(side = 1, at = atvec, labels = FALSE)
axis(side = 1, at = atv, tick = FALSE,
     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                "Sep", "Oct", "Nov", "Dec"))

#######################

fff <-t(XShift[, 366:728] %*% t(beta21)[366:728, ])


plot(c(1, 365), c(min(fff), max(fff)), type = 'n', main = "Walker County, 2018", 
     ylab = 'Sustained Shift in the Mean', 
     xlab = "", xaxt = "n")

for (i in 1:nsim) {
  points(as.vector(fff[i, ]), col = 'gray', type = "l")
}


points(colMeans(fff), col = 'black', type = "l")


axis(side = 1, at = atvec, labels = FALSE)
axis(side = 1, at = atv, tick = FALSE,
     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                "Sep", "Oct", "Nov", "Dec"))

#######################



plot(c(min(beta21[, 366:365]), max(beta21[, 1:365])), type = 'n', main = "Walker County, 2018", 
     ylab = 'Trend in the Mean', 
     xlab = "", xaxt = "n")

for (i in 1:nsim) {
  points(beta21[i, 1:365], type = 'l', col = 'skyblue')
}


points(colMeans(beta21[, 1:365]), type = 'l', col = 'blue')

#######################


#debug(SimLinear)
SimLinear(Y, beta0, beta1, initial = V[1, ])

fit0 <- t(cbind(1, V) %*% t(betamat1[, c(1, (dim(betamat1)[2] - p + 1):dim(betamat1)[2])]))
fit1 <- t(cbind(1, XShift, V) %*% t(betamat1))

#debug(getResiCS)
qq <- getResiCS(Y, fit0, nsim)

#debug(findrootResiCSPPP)
ee1 <- findrootResiCSPPP(0.2, Y, fit0, nsim, nnsim = 1000, interval = c(1e-4, 20), tol = 1e-4)

ee2 <- findrootResiCSPPP(0.5, Y, fit0, nsim, nnsim = 1000, interval = c(1e-4, 1000), tol = 1e-4)
  
dd <- getResiCSAlt(Y, fit0, fit1, nsim)








qq <- getKSlogpdfCSPPP(Y, fit0, nsim)
qq <- getKSResiCSPPP(Y, fit0, nsim)

qq <- getKLResiCSPPP(Y, fit0, nsim)




pp <- getKSResiCSPPPAlt(Y, fit0, fit1, nsim) 

ee <- findrootKLResiCSPPP(0.05, Y, fit0, nsim, nnsim = 100)

nnsim <- 100

cc <- findroot(0.2, Y, fit0, nsim, nnsim, interval = c(0, 5), tol = 1e-4)
 
test <- getKSCS(Y, fit0, fit1, nsim)


cc <- findrootlogpdf(0.1, Y, fit0, nsim, nnsim)

test <- getKSlogpdfCS(Y, fit0, fit1, nsim)

nnsim <- 10000

FAP <- 0.1

findroot(0.1, Y, fit0, nsim)


r <- matrix(NA, nrow = nsim, ncol = length(cntma))

for (i in 1:nsim) {
  
  r[i, ] <- (Y - fit0mat1[i, ]) ^ 2 / var(Y - fit0mat1[i, ])

}

##########################

d <- rep(NA, nsim)
D <- rep(NA, nsim)

for (i in 1:nsim) {
  
  d[i] <- sum((((Y - fit0mat1[i, ]) ^ 2 / var(Y - fit0mat1[i, ]) - (Y - fit1mat1[i, ])) ^ 2 / sigmamat1[i]))
  Yrep <- rnorm(length(cntma), fit1mat1[i, ], sqrt(sigmamat1[i]))
  D[i] <- sum((((Yrep - fit0mat1[i, ]) ^ 2 / var(Y - fit0mat1[i, ]) - (Yrep - fit1mat1[i, ])) ^ 2/ sigmamat1[i]) )
  
}



##########################

ee1 <- colMeans(betamat1)
ee2 <- betamat1

pp <- dim(betamat1)[2]

for (i in 1:pp) {
  ee2[, i] <- (ee2[, i] - ee1[i]) / sqrt(var(betamat1[, i]))
  ee1[i] <- ee1[i] / sqrt(var(betamat1[, i]))
}
