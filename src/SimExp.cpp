#include <RcppArmadillo.h>      // also pulls in Rcpp.h amd cmath

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


double getExp(double& lambda2, double& sigma2) {
  return Rcpp::rexp(1, sqrt(lambda2 / sigma2))[0];
}


double getTrunExp(double& lambda2, double& sigma2, double& UpperBound) {
  double randProb = Rcpp::runif(1, 0, 1)[0];
  return log(randProb * (1 - exp(- sqrt(lambda2 / sigma2) * UpperBound)) + 
             exp(-sqrt(lambda2 / sigma2) * UpperBound)) / (- sqrt(lambda2 / sigma2));
}

// [[Rcpp::export]]
Rcpp::NumericVector getBetaPlusOrMinus(int& p, double& lambda2, double& sigma2) {
  Rcpp::NumericVector BetaVector (p);
  double tmp;
  double upperBound;
  
  for (int i = 0; i < p; i++) {
    if (i == 0) {
      upperBound = 1;
      tmp = getExp(lambda2, sigma2);
    } else {
      upperBound = BetaVector(i - 1);
      tmp = getTrunExp(lambda2, sigma2, upperBound);
    }
    BetaVector(i) = tmp;
  }
  
  return BetaVector;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix getGammaPlusOrMinus(int& n, int& p1, int& p, double& lambda2, double& sigma2, 
                                        Rcpp::NumericVector& BetaVector) {
  //Rcpp::NumericVector GammaVector (p1 * (n - 2));
  int tmpp;
  double UpperBound;
  double tmpGammaDif;
  int k;
  
  if (p1 <= p) {
    tmpp = p1;
  } else {
    tmpp = p;
    // p1 must be smaller than p
  }
  
  Rcpp::NumericMatrix GammaMat(n, tmpp);
  
  for (int i = 1; i < (n - 1); i++) {
    if (i >= tmpp) {
      k = tmpp;
    } else {
      k = i;
    }
    for (int j = 0; j < k; j++) {
      tmpGammaDif = 0;
      if (j > 0) {
        for (int t = 0; t < i; t++) {
          tmpGammaDif = tmpGammaDif + GammaMat(t, j - 1) - GammaMat(t, j);
        }
        UpperBound = BetaVector(j - 1) - BetaVector(j) + tmpGammaDif + GammaMat(i, j - 1);
        GammaMat(i, j) = getTrunExp(lambda2, sigma2, UpperBound);
      } else if (j == 0) {
        UpperBound = 1 - BetaVector(j);
        GammaMat(i, j) = getExp(lambda2, sigma2);
      }
    }
  }
  
  
  return GammaMat;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix getDeltaPlusOrMinus(int& n, int& p2, int& p, 
                                        double& lambda2, double& sigma2) {
  //Rcpp::NumericVector GammaVector (p1 * (n - 2));
  int tmpp;
  int k;
  
  if (p2 <= p) {
    tmpp = p2;
  } else {
    tmpp = p;
    // p1 must be smaller than p
  }
  
  Rcpp::NumericMatrix DeltaMat(n, tmpp);
  
  for (int i = 0; i < n; i++) {
    if (i >= tmpp) {
      k = tmpp;
    } else {
      k = i;
    }

    for (int j = 0; j < k; j++) {
      DeltaMat(i, j) = getExp(lambda2, sigma2);
    }
  }
  
  
  return DeltaMat;
}

double getSign(double x) {
  double tmp;
  if (x > 0) {
    tmp = 1;
  } else if (x < 0) {
    tmp = -1;
  } else {
    tmp = 0;
  }
  return tmp;
}

// [[Rcpp::export]]
double getLaplace(double& lambda2, double& sigma2) {
  double randProb = Rcpp::runif(1, 0, 1)[0];
  double b = sqrt(lambda2 / sigma2);
  double tmpsign = getSign(randProb - 0.5);
  return -b * tmpsign * log(1 - 2 * abs(randProb - 0.5));
}

// [[Rcpp::export]]
Rcpp::NumericVector getGamma0(int& n, double& lambda2, double& sigma2) {
  //Rcpp::NumericVector GammaVector (p1 * (n - 2));
 
  Rcpp::NumericVector Gamma0Vector(n);
  
  for (int i = 1; i < (n - 1); i++) {
    Gamma0Vector(i) = getLaplace(lambda2, sigma2);
  }
  
  return Gamma0Vector;
}

// [[Rcpp::export]]
Rcpp::NumericVector getDelta0(int& n, double& lambda2, double& sigma2) {
  //Rcpp::NumericVector GammaVector (p1 * (n - 2));
  
  Rcpp::NumericVector Delta0Vector(n);
  
  for (int i = 0; i < n; i++) {
    Delta0Vector(i) = getLaplace(lambda2, sigma2);
  }
  
  return Delta0Vector;
}


// [[Rcpp::export]]
double getSigma2(double& a, double& b) {
  //Rcpp::NumericVector GammaVector (p1 * (n - 2));
  
  return 1 / Rcpp::rgamma(1, a, 1/b)(0);
  
}

// [[Rcpp::export]]
double getLambda2(double& a, double& b) {
  //Rcpp::NumericVector GammaVector (p1 * (n - 2));
  
  return Rcpp::rgamma(1, a, 1/b)(0);
  
}

// [[Rcpp::export]]
Rcpp::NumericVector getResidual(Rcpp::NumericVector& YTilde, 
                                int& p, int& p1, int& p2, 
                          Rcpp::NumericVector& BetaPlus, 
                          Rcpp::NumericVector& BetaMinus,
                          Rcpp::NumericVector& GammaPlus,
                          Rcpp::NumericVector& GammaMinus,
                          Rcpp::NumericVector& DeltaPlus,
                          Rcpp::NumericVector& DeltaMinus,
                          Rcpp::NumericVector& Gamma0,
                          Rcpp::NumericVector& Delta0,
                          double& lambda2, double& sigma2) {
  
  int n = YTilde.length();
  Rcpp::NumericVector ResiVector(n);
  Rcpp::NumericVector Ylag(p);
  Rcpp::NumericMatrix YlagMatrix(n, p);
  
  int k;
  Rcpp::NumericVector tmpGammaPlusSum(p1);
  Rcpp::NumericVector tmpGammaMinusSum(p1);
  Rcpp::NumericMatrix GammaPlusSum(n, p);
  Rcpp::NumericMatrix GammaMinusSum(n, p);
  Rcpp::NumericVector Gamma0Sum(n);
  double tmpGamma0Sum = 0;
  Rcpp::NumericMatrix DeltaPlusM(n, p);
  Rcpp::NumericMatrix DeltaMinusM(n, p);
  double tmp;
  
  for (int i = 0; i < n; i++) {

    tmpGamma0Sum = tmpGamma0Sum + Gamma0(i);
    Gamma0Sum(i) = tmpGamma0Sum;
    
    for (int j = 0; j < p; j++) {
      k = i - (j + 1) ;
      if (k >= 0) {
        Ylag(j) = YTilde(k);
      }
    }
    
    YlagMatrix(i, _) = Ylag;
    
    for (int j = 0; j < p1; j++) {
      tmpGammaPlusSum(j) = tmpGammaPlusSum(j) + GammaPlus(i, j);
      tmpGammaMinusSum(j) = tmpGammaMinusSum(j) + GammaMinus(i, j);
      GammaPlusSum(i, j) = tmpGammaPlusSum(j);
      GammaMinusSum(i, j) = tmpGammaMinusSum(j);
    }
    
    for (int j = 0; j < p2; j++) {
      DeltaPlusM(i, j) = DeltaPlus(i, j);
      DeltaMinusM(i, j) = DeltaMinus(i, j);
    }

    //tmp = 0;
    tmp = Gamma0Sum(i) + Delta0(i);
    for (int j = 0; j < p; j++) {
      tmp = tmp + YlagMatrix(i, j) * (BetaPlus(j) - BetaMinus(j));
      tmp = tmp + YlagMatrix(i, j) * (GammaPlusSum(i, j) - GammaMinusSum(i, j));
      tmp = tmp + YlagMatrix(i, j) * (DeltaPlusM(i, j) - DeltaMinusM(i, j));
    }
    
    ResiVector(i) = YTilde(i) - tmp;
    
  }
  
  return ResiVector;
  
}

// [[Rcpp::export]]
double getLikelihood(Rcpp::NumericVector& ResiVector, double& sigma2){
  double out = exp(sum(log(Rcpp::dnorm(ResiVector, 0, sqrt(sigma2)))));
  return out;
}





// [[Rcpp::export]]
Rcpp::List getProb(Rcpp::NumericVector& YTilde, int& p, int& p1, int& p2, 
                     double& a1, double& b1, double& a2, double& b2, int& nsim, 
                     double& tol){
  double lambda2;
  double sigma2;
  int n = YTilde.length();
  Rcpp::NumericVector BetaPlus(p);
  Rcpp::NumericVector BetaMinus(p);
  Rcpp::NumericMatrix GammaPlus(n, p1);
  Rcpp::NumericMatrix GammaMinus(n, p1);
  Rcpp::NumericMatrix DeltaPlus(n, p2);
  Rcpp::NumericMatrix DeltaMinus(n, p2);
  Rcpp::NumericVector Gamma0(n);
  Rcpp::NumericVector Delta0(n);
  Rcpp::NumericVector tmpResidualVector(n);
  Rcpp::NumericVector simLikelihood(nsim);
  Rcpp::NumericVector outlier(nsim);
  double tmp;
  double out;
  
  for (int v = 0; v < nsim; v++) {
    lambda2 = getLambda2(a1, b1);
    sigma2 = getSigma2(a2, b2);
    BetaPlus = getBetaPlusOrMinus(p, lambda2, sigma2);
    BetaMinus = getBetaPlusOrMinus(p, lambda2, sigma2);
    GammaPlus = getGammaPlusOrMinus(n, p1, p, lambda2, sigma2, BetaPlus);
    GammaMinus = getGammaPlusOrMinus(n, p1, p, lambda2, sigma2, BetaMinus);
    DeltaPlus = getDeltaPlusOrMinus(n, p2, p, lambda2, sigma2);
    DeltaMinus = getDeltaPlusOrMinus(n, p2, p, lambda2, sigma2);
    Gamma0 = getGamma0(n, lambda2, sigma2);
    Delta0 = getDelta0(n, lambda2, sigma2);
    tmpResidualVector = getResidual(YTilde, p, p1, p2, BetaPlus, BetaMinus,
                                    GammaPlus, GammaMinus, DeltaPlus, DeltaMinus,
                                    Gamma0, Delta0, lambda2, sigma2);
    simLikelihood(v) = getLikelihood(tmpResidualVector, sigma2);
    
    tmp = 0;
    for (int i = 0; i < n; i++) {
      if (abs(Gamma0(i)) > tol) {
        tmp = tmp + 1;
      }
      if (abs(Delta0(i)) > tol) {
        tmp = tmp + 1;
      }
      
      for (int j = 0; j < p1; j++) {
        if (abs(GammaPlus(i, j) - GammaMinus(i, j)) > tol) {
          tmp = tmp + 1;
        }
      }
      
      for (int j = 0; j < p2; j++) {
        if (abs(DeltaPlus(i, j) - DeltaMinus(i, j)) > tol) {
          tmp = tmp + 1;
        }
      }
    }
    
    if (tmp > 1) {
      outlier(v) = 1;
    }
    
  }
  
  out = 1 - mean(outlier * simLikelihood / mean(simLikelihood));
  
  return Rcpp::List::create(_["out"] = out, 
                            _["outlier"] = outlier,
                           _["simLikelihood"] = simLikelihood);
  
}


// [[Rcpp::export]]
Rcpp::NumericMatrix getW(int& n) {
  
  int p = 0;
  int m = n - 2 - p;
  Rcpp::NumericMatrix out(n, m);
  int k = -1;
  
  for (int j = 0; j < m; j++) {
  
    k++;
    
    for (int i = 0; i < n; i++) {

      if (i > (p + k)) {
        out(i, j) = 1;
      }

      
    }
    
  }

  return(out);
  
}

// [[Rcpp::export]]
Rcpp::NumericMatrix getV(int& n) {
  
  NumericMatrix out(n, n);
  
  for (int j = 0; j < n; j++) {
    
    for (int i = 0; i < n; i++) {
      
      if (i == j) {
        out(i, j) = 1;
      }
      
      
    }
    
  }
  
  return(out);
  
}


// [[Rcpp::export]]
Rcpp::NumericMatrix getYp(Rcpp::NumericVector& Y, int& n, int& p) {
  
  NumericMatrix out(n, p);
  int tmp;
  
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
 
      tmp = i - j - 1;
      if (tmp >= 0) {
        out(i, j) = Y(tmp);
      }
      
      
    }
    
  }
  
  return(out);
  
}

// [[Rcpp::export]]
Rcpp::NumericMatrix getXNull(Rcpp::NumericVector& Y, int& n, int& p) {
  
  NumericMatrix Yp = getYp(Y, n, p);
  
  int tmp = p;
  int k = 0;
  
  tmp = 2 * tmp;
  
  NumericMatrix out(n, tmp);
  
  k = 0;

  
  //laggy Y Plus
  for (int i = 0; i < p; i++) {
    out(_, k) = Yp(_, i);
    k++;
  }
  
  //laggy Y Minus
  for (int i = 0; i < p; i++) {
    out(_, k) = -Yp(_, i);
    k++;
  }
  
  return(out);
  
}

// [[Rcpp::export]]
Rcpp::NumericMatrix getX(Rcpp::NumericVector& Y, int& n, int& p, int& p1, int& p2) {
  
  NumericMatrix Yp = getYp(Y, n, p);

  int tmp = p;
  int k = 0;

    NumericMatrix W = getW(n);
    NumericMatrix V = getV(n);
    
    tmp = 2 * tmp + 2 * n - 2;

    int r1 = 0;
    
    if (p1 > 0) {
      for (int i = 0; i < p1; i++) {
        for (int j = (i + 1); j < (n - 2); j++) {
          r1++;
        }
      }
    }
    
    int r2 = 0;
    
    if (p2 > 0) {
      for (int i = 0; i < p2; i++) {
        for (int j = (i + 1); j < n; j++) {
          r2++;
        }
      }
    }

    tmp = tmp + 2 * r1 + r2;
    
    NumericMatrix out(n, tmp);
    
    k = 0;

    //Beta0 step shift
    for (int i = 0; i < (n - 2); i++) {
      out(_, k) = W(_, i);
      k++;
    }
    
    //Beta0 isolated shift
    for (int i = 0; i < n; i++) {
      out(_, k) = V(_, i);
      k++;
    }
  
  
    //laggy Y Plus
    for (int i = 0; i < p; i++) {
      out(_, k) = Yp(_, i);
      k++;
    }
    
    //laggy Y Minus
    for (int i = 0; i < p; i++) {
      out(_, k) = -Yp(_, i);
      k++;
    }
    
    if (p1 > 0) {
      //Beta1 step shift Plus
      for (int i = 0; i < p1; i++) {
        for (int j = (i + 1); j < (n - 2); j++) {
          out(_, k) = Yp(_, i) * W(_, j);
          k++;
        }
      }
      
      //Beta1 step shift Minus
      for (int i = 0; i < p1; i++) {
        for (int j = (i + 1); j < (n - 2); j++) {
          out(_, k) = -Yp(_, i) * W(_, j);
          k++;
        }
      }
    }
      
      
      
    if (p2 > 0) {
      //Beta1 isolated shift
      for (int i = 0; i < p2; i++) {
        for (int j = (i + 1); j < n; j++) {
          out(_, k) = Yp(_, i) * V(_, j);
          k++;
        }
      }
    }

  return(out);
  
}

// [[Rcpp::export]]
arma::vec convertNumericVectorToArmaVec (Rcpp::NumericVector& A) {
  int n = A.length();
  arma::vec A1(n);
  for (int i = 0; i < n; i++) {
    A1(i) = A(i);
  }
  return(A1);
}

// [[Rcpp::export]]
arma::mat convertNumericMatrixToArmaMat (Rcpp::NumericMatrix& A) {
  int n = A.nrow();
  int m = A.ncol();
  arma::mat A1(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      A1(i, j) = A(i, j);
    }
  }
  return(A1);
}

// [[Rcpp::export]]
double getSSE(arma::vec& Y, arma::mat& X, arma::vec& Beta) {
  
  arma::vec resi;
  arma::vec tmp(1);

  resi = Y - X * Beta;
  
  tmp = resi.t() * resi;
  
  return(tmp[0]);
  
}

// [[Rcpp::export]]
arma::vec getInitialBetaLSNull(Rcpp::NumericVector& Y, int& n, int& p) {
  // This function is using optimizer in R, 
  // and needs to be changed using C++ in the future.
  Rcpp::NumericMatrix Yp = getYp(Y, n, p);
  arma::mat Yp1 = convertNumericMatrixToArmaMat(Yp);
  arma::vec Y1 = convertNumericVectorToArmaVec(Y);
  arma::vec Beta1 = inv(Yp1.t() * Yp1) * Yp1.t() * Y1;
  arma::vec Beta1plus = Beta1;
  arma::vec Beta1minus = Beta1;
  //rebuild the beta vector
  int tmp = p;
  
  tmp = 2 * tmp;
  
  arma::vec BetaOut(tmp);
  BetaOut.zeros();
  
  for (int i = 0; i < p; i++) {
    if (Beta1plus(i) < 0) {
      Beta1plus(i) = 0;
    }
  }
  
  for (int i = 0; i < p; i++) {
    if (Beta1minus(i) > 0) {
      Beta1minus(i) = 0;
    }
  }
  
  int k1 = 0;
  int k2 = 0;
  for (int i = 0; i < tmp; i++) {
    
    if ((i >= 0) & (i < p)) {
      BetaOut(i) = Beta1plus(k1);
      k1++;
    } 
    
    if ((i >= p) & (i < tmp)) {
      BetaOut(i) = -Beta1minus(k2);
      k2++;
    }
    
  }
  
  return(BetaOut);
  
}

// [[Rcpp::export]]
arma::vec getInitialBetaLS(Rcpp::NumericVector& Y, int& n, int& p, int& p1, int& p2) {
  // This function is using optimizer in R, 
  // and needs to be changed using C++ in the future.
  Rcpp::NumericMatrix Yp = getYp(Y, n, p);
  arma::mat Yp1 = convertNumericMatrixToArmaMat(Yp);
  arma::vec Y1 = convertNumericVectorToArmaVec(Y);
  arma::vec Beta1 = inv(Yp1.t() * Yp1) * Yp1.t() * Y1;
  arma::vec Beta1plus = Beta1;
  arma::vec Beta1minus = Beta1;
  //rebuild the beta vector
  int tmp = p;
  
  tmp = 2 * tmp + 2 * n - 2;
  
  int r1 = 0;
  
  if (p1 > 0) {
    for (int i = 0; i < p1; i++) {
      for (int j = (i + 1); j < (n - 2); j++) {
        r1++;
      }
    }
  }
  
  int r2 = 0;
  
  if (p2 > 0) {
    for (int i = 0; i < p2; i++) {
      for (int j = (i + 1); j < n; j++) {
        r2++;
      }
    }
  }
  
  tmp = tmp + 2 * r1 + r2;
  
  arma::vec BetaOut(tmp);
  BetaOut.zeros();
  
  //Gamma0 and Delta0 goes first, then Beta1 plus and Beta1 minus, then Gamma1 plus
  //and Gamma1 minus, and then Delta1.
  for (int i = 0; i < p; i++) {
    if (Beta1plus(i) < 0) {
      Beta1plus(i) = 0;
    }
  }
  
  for (int i = 0; i < p; i++) {
    if (Beta1minus(i) > 0) {
      Beta1minus(i) = 0;
    }
  }
  
  int k1 = 0;
  int k2 = 0;
  for (int i = 0; i < tmp; i++) {
  
    if ((i >= (2 * n - 2)) & (i < (2 * n - 2 + p))) {
      BetaOut(i) = Beta1plus(k1);
      k1++;
    } 
    
    if ((i >= (2 * n - 2 + p)) & (i < (2 * n - 2 + 2 * p))) {
      BetaOut(i) = -Beta1minus(k2);
      k2++;
    }
    
  }
  
  return(BetaOut);
  
}

// [[Rcpp::export]]
arma::vec getOmega(arma::vec& Beta, double& lambda) {
  int n = Beta.n_elem ;
  double tmp = 0;
  arma::vec Omega(n);
  
  for (int i = 0; i < n; i++) {
    tmp = abs(Beta(i));
    Omega(i) = (1 / lambda) * tmp;
  }
  return Omega;
}

// [[Rcpp::export]]
arma::vec getBetaLASSO(arma::vec& Y, arma::mat& X, arma::vec& Omega) {
  int n = Omega.n_elem;
  arma::vec Beta(n);
  arma::mat OmegaMat(n, n);
  
  Beta.zeros();
  OmegaMat.zeros();

  for (int i = 0; i < n; i++) {
    OmegaMat(i, i) = sqrt(Omega(i));
  }
  
  Beta = OmegaMat * arma::inv(arma::eye(n, n) + OmegaMat * X.t() * X * OmegaMat) * 
    OmegaMat * X.t() * Y;
  
  return Beta;
}

// [[Rcpp::export]]
double invGaussian(double& mu, double& lambda) {
  double v = Rcpp::rnorm(1, 0, 1)(0);  // Sample from a normal distribution with a mean of 0 and 1 standard deviation
  double y = v * v;
  double x = mu + (mu * mu * y) / (2 * lambda) - (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu * mu * y * y);
  double test = Rcpp::runif(1, 0, 1)(0);  // Sample from a uniform distribution between 0 and 1
  if (test <= (mu) / (mu + x))
    return x;
  else
    return (mu * mu) / x;
}

// [[Rcpp::export]]
arma::vec rinvGaussian(int& n, double& mu, double& lambda) {
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    out(i) = invGaussian(mu, lambda);
  }
  
  return out;
}


// [[Rcpp::export]]
double invGamma(double& kappa, double& theta) {
  double out = 1 / Rcpp::rgamma(1, kappa, 1 / theta)(0);
  return out;
}

// [[Rcpp::export]]
arma::vec rinvGamma(int& n, double& kappa, double& theta) {
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    out(i) = invGamma(kappa, theta);
  }
  
  return out;
}

// [[Rcpp::export]]
double convertLambdaToLambdaPC(double& lambda, double& sigma2) {
  double lambdaPC = lambda / sqrt(sigma2);
  return lambdaPC;
}

// [[Rcpp::export]]
double convertLambdaPCToLambda(double& lambdaPC, double& sigma2) {
  double lambda = lambdaPC * sqrt(sigma2);
  return lambda;
}

// [[Rcpp::export]]
double getLambdaEM(arma::vec& Beta, double& lambda, double& sigma2, 
                   int& nsim, double& tol) {
  int c = Beta.n_elem;
  double lambdaPC = convertLambdaToLambdaPC(lambda, sigma2);
  arma::vec denom(c);
  double tmp;
  double tmp1 = 0;
  double tmp2 = 0;
  
  for (int i = 0; i < c; i++) {
    tmp = 0;
    if (abs(Beta(i)) <= tol) {
      tmp = 1 / (2 / (lambdaPC * lambdaPC));
      denom(i) = 1 / 2 * tmp * tmp;
    } else {
      tmp1 = lambdaPC * sqrt(sigma2) / abs(Beta(i));
      tmp2 = lambdaPC / sqrt(sigma2);
      for (int j = 0; j < nsim; j++) {
        tmp = tmp + 1 / 
          invGaussian(tmp1, tmp2);
      }
      denom(i) = tmp / nsim;
    }
  }

  double lambdaPCNew = sqrt(2 * c / sum(denom));
  double lambdaNew = convertLambdaPCToLambda(lambdaPCNew, sigma2);
  return lambdaNew;
}

// [[Rcpp::export]]
Rcpp::List getPosterior(arma::vec& Y, arma::mat& X, 
                       arma::vec& InitialOmega,
                       arma::vec& InitialBeta, 
                       double& InitialLambda, 
                       double& InitialSigma2, 
                       int& nsimInvGaussian, 
                       int& nsim, int& burnin, double& tol) {
  int n = Y.n_rows;
  int nsimTot = nsim + burnin;
  int p = InitialBeta.n_elem;
  arma::mat outOmega(nsim, p);
  arma::mat outBeta(nsim, p);
  arma::vec outSigma2(nsim);
  arma::vec outLambda(nsim);

  arma::vec Omega = InitialOmega;
  arma::vec Beta = InitialBeta;
  double sigma2 = InitialSigma2;
  double lambda = InitialLambda;
  
  int k = 0;
  for (int i = 0; i < (nsimTot); i++) {
    Omega = getOmega(Beta, lambda); 
    Beta = getBetaLASSO(Y, X, Omega);
    sigma2 = getSSE(Y, X, Beta) / n;
    lambda = getLambdaEM(Beta, lambda, sigma2, nsimInvGaussian, tol);
    
    if (i >= burnin) {
      for (int j = 0; j < p; j++) {
        outOmega(k, j) = Omega(j);
        outBeta(k, j) = Beta(j);
      } 
      outSigma2(k) = sigma2;
      outLambda(k) = lambda;
      k++;
    }
  }
  Rcpp::List out;
  return out = Rcpp::List::create(
    _["Omega"] = outOmega,
    _["Beta"] = outBeta,
    _["sigma2"] = outSigma2,
    _["lambda"] = outLambda
  );
  
}


//arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
//  int ncols = sigma.n_cols;
//  arma::mat Y = arma::randn(n, ncols);
//  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
//}



////Rcpp::List optim_rcpp(arma::vec& Beta,
////                     arma::vec& Y, arma::mat& X){
////  
////  // Extract R's optim function
////  Rcpp::Environment stats("package:stats"); 
////  Rcpp::Function optim = stats["optim"];
////  
////  // Call the optim function from R in C++ 
////  Rcpp::List opt_results = optim(Rcpp::_["par"]    = Beta,
////                                 // Make sure this function is not exported!
////                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&getSSE),
////                                 Rcpp::_["method"] = "BFGS",
////                                 // Pass in the other parameters as everything
////                                 // is scoped environmentally
////                                 Rcpp::_["Y"] = Y,
////                                 Rcpp::_["X"] = X);
////  
////  // Return estimated values
////  return opt_results;
////  
////}