library(GIGrvg)
library(breakfast)
library(pscl)


#dat <-read.csv(file = "/Users/yuihuiyao/Library/CloudStorage/Box-Box/Yuhui R21/Walker County De-Identified 2016-2021 Opioid ER Visits.csv")

dat <- read.csv(file = "C:/Users/Yuhui/Box/Yuhui R21/Walker County De-Identified 2016-2021 Opioid ER Visits.csv")

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

w <- 28

cntma <- rep(NA, dim(dat1)[1])

for (i in 1:dim(dat1)[1]) {
  
  if (i < w) {
    cntma[i] <- mean(as.numeric(dat1[1:i, 2]))
  } else {
    cntma[i] <- mean(as.numeric(dat1[(i - w + 1):i, 2]))
  }
  
}

#############
p <- 5
V <- getT(cntma, p)

check <- which(as.Date("2017-01-01") <= as.Date(dat1[, 1]) & as.Date(dat1[, 1]) <= as.Date("2020-12-31"))

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

########################

Y <- cntma;


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



getKLResiCSPPP <- function(Y, fit0, nsim, breaks = 10) {
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
    tmp <- getBinProb(tmpy, tmpx, breaks = breaks)
    out[i] <- KL(tmp, unit = "log")
  }
  out
}

findrootKLResiCSPPP <- function(FAP0, Y, fit0, nsim, nnsim, breaks = 20, interval = c(0, 5), tol = 1e-4) {
  
  rootfinding <- function(cc, FAP0, ref) {
    check <- ref <= cc
    ProbNSE <- mean(rowMeans(check))
    diff <- 1 - FAP0 - ProbNSE
    cat("cc:", cc, "and diff:", diff, '\n')
    return(diff)
  }
  
  ref <- matrix(NA, nrow = nnsim, ncol = length(Y))
  
  for (i in 1:nnsim) {
    ref[i, ] <- getKLResiCSPPP(Y, fit0, nsim, breaks) 
  }
  
  uniroot(rootfinding, interval = interval, FAP0 = FAP0, 
          ref = ref, tol = tol)$root
  
}

getKLResiCSPPPAlt <- function(Y, fit0, fit1, nsim, breaks = 10) {
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
    tmp <- getBinProb(tmpy, tmpx, breaks = breaks)
    out[i] <- KL(tmp, unit = "log")
  }
  out
}


fit0 <- t(cbind(1, V) %*% t(betamat1[, c(1, (dim(betamat1)[2] - p + 1):dim(betamat1)[2])]))

qq <- getKSlogpdfCSPPP(Y, fit0, nsim)
qq <- getKSResiCSPPP(Y, fit0, nsim)

qq <- getKLResiCSPPP(Y, fit0, nsim)


fit1 <- t(cbind(1, XShift, V) %*% t(betamat1))

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
