// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#define ARMA_WARN_LEVEL 1
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace Rcpp;
using namespace arma;

arma::colvec rinvgaussiancpp(int n, double mu, double lambda){
 Rcpp::Environment pkg = Rcpp::Environment::namespace_env("VGAM");
 
 // Picking up Matrix() function from Matrix package
 Rcpp::Function rinvgaussian = pkg["rinv.gaussian"];
 
 Rcpp::NumericVector  tmp = rinvgaussian(n, mu, lambda);
 arma::colvec out = Rcpp::as<arma::colvec>(tmp); 
 return out;
}

// [[Rcpp::export]]
arma::mat getV(arma::colvec Y, int q) {
  int T = Y.n_elem;
  arma::mat out(T, q);
  out.zeros();
  
  int tmp;
  
  int i;
  int j;
  for (i = 1; i < T; i++) {
    
    if (i > q) {
      tmp = q;
    } else {
      tmp = i;
    }
    
    for (j = 0; j < tmp; j++) {
      out(i, j) = Y(i - j - 1);
    }
  }
  return(out);
}

arma::colvec rmvnorm(arma::colvec Mean, arma::mat Sigma) {
  
  int q = Mean.n_elem;
  arma::colvec Z = arma::randn(q);
  arma::colvec out = Mean + arma::chol(Sigma) * Z;
  return(out);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @param a is the lower bound in absolute value.
//' @param b is the upper bound in absolute value.
//' @param mean is a mean.
//' @param sd is a standard deviation.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec rtwosegnorm(int n, double a, double b, double mean, double sd) {
  arma::colvec out(n);
  arma::colvec U = arma::randu(n);
  
  double a1 = (-a - mean) / sd;
  double b1 = (-b - mean) / sd;
  double a2 = (a - mean) / sd;
  double b2 = (b - mean) / sd;
  
  double p1 = arma::normcdf(a1) - arma::normcdf(b1);
  double p2 = arma::normcdf(b2) - arma::normcdf(a2);
  double p = p1 + p2;
  
  int i;
  if (p > 0) {
    for (i = 0; i < n; i++) {
      if (U(i) <= p1 / p) {
        out(i) = R::qnorm5(U(i) * p + arma::normcdf(b1), 0.0, 1.0, 1, 0) * sd + mean;
      } else {
        out(i) = R::qnorm5(U(i) * p + arma::normcdf(a2) - p1, 0.0, 1.0, 1, 0) * sd + mean;
      }
    }
    
  }
  return(out);
  
}



// [[Rcpp::export]]
arma::mat getGMat(int T, int q) {
  arma::mat tmp(T, T);
  tmp.eye();
  
  arma::mat out(T - q, T);
  out = tmp.submat(q, 0, T - 1, T - 1);
  return(out);
}

// [[Rcpp::export]]
arma::mat getPhiMat(arma::colvec Phi, int T) {
  int q = Phi.n_elem;
  arma::mat tmp(T, T);
  tmp.zeros();
  
  int i;
  for (i = 0; i < q; i++) {
    tmp.diag(-i - 1).fill(Phi(i));
  }
  
  arma::mat out(T - q, T);
  out = tmp.submat(q, 0, T - 1, T - 1);
  return(out);
}

//' Design Matrix (MT)
//' 
//' gets a design matrix as that in MT
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatMT(100, 5)
// [[Rcpp::export]]
arma::mat getHMatMT(int T, int q) {
  arma::mat tmp(T, T);
  tmp.ones();
  arma::mat L = arma::trimatl(tmp);
  
  arma::mat out(T, T - q);
  out = L.submat(0, q, T - 1, T - 1);
  return(out);
  
}

//' Design Matrix for Sustained Shift (CM)
//' 
//' gets a design matrix for sustained shift as that in CM
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatSustained(100, 5)
// [[Rcpp::export]]
arma::mat getHMatSustained(int T, int q) {
  int w = 1;
  arma::mat tmp(T, T);
  tmp.ones();
  arma::mat L = arma::trimatl(tmp);
  arma::mat out = L.cols(q + w, T - 2 + 1);
  int nn = T - 2 - w + 1 - q;
  arma::mat out1(T, nn);
  int rr = 0;
  int kk = 0;
  int i;
  for (i = 0; i < nn; i++) {
    
    if (rr == 0) {
      out1.col(kk) = out.col(i);
      kk = kk + 1;
      rr = rr - w;
    }
    rr = rr + 1;
  }
  int tmpn = floor(nn / w);
  return(out1.cols(0, tmpn - 1));
}

//' Design Matrix for Isolated Shift (CM)
//' 
//' gets a design matrix for isolated shift as that in CM
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatIsolated(100, 5)
// [[Rcpp::export]]
arma::mat getHMatIsolated(int T, int q) {
  int w = 1;
  arma::mat tmp(T, T);
  
  int i;
  for (i = 0; i < w; i++) {
    tmp.diag(-i).ones();
  }
  
  arma::mat out = tmp.cols(q, T - 1);
  int nn = T - 1 - q + 1;
  arma::mat out1(T, nn);
  int rr = 0;
  int kk = 0;
  
  for (i = 0; i < nn; i++) {
    
    if (rr == 0) {
      out1.col(kk) = out.col(i);
      kk = kk + 1;
      rr = rr - w;
    }
    rr = rr + 1;
  }
  int tmpn = floor(nn / w);
  return(out1.cols(0, tmpn - 1));
}

//' gets a design matrix for gradual shift
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatGradual(100, 5)
// [[Rcpp::export]]
arma::mat getHMatGradual(int T, int q) {
  arma::mat tmp(T, T);
  int w = 1;
  int i;
  for (i = 0; i < T; i++) {
    tmp.diag(-i).fill(i);
  }
  
  arma::mat out = tmp.cols(q, T - 2);
  int nn = T - 2 - q + 1;
  arma::mat out1(T, nn);
  int rr = 0;
  int kk = 0;
  
  for (i = 0; i < nn; i++) {
    
    if (rr == 0) {
      out1.col(kk) = out.col(i);
      kk = kk + 1;
      rr = rr - w;
    }
    rr = rr + 1;
  }
  int tmpn = floor(nn / w);
  return(out1.cols(0, tmpn - 1));
  
}

//' gets a design matrix for Fourier series
//'
//' @param T is length of a process.
//' @param s is the number of period
//' @param n is the number of Fourier series
//' @export
//' @examples
//' getXSeasonalityFS(100, 15, 10)
// [[Rcpp::export]]
arma::mat getXSeasonalityFS(int T, double s, int n) {
  double pi = 3.141592653589793238462643383280;
  
  arma::mat X(T, 1);
  arma::mat out(T, 2 * n);
  
  int i;
  for (i = 0; i < T; i++) {
    X(i, 0) = i + 1.0;
  }
  
  int k = 0;
  int j;
  for (j = 0; j < n; j++) {
    out.col(k) = arma::cos(2.0 * pi * (j + 1.0) * X / s);
    k = k + 1;
    
    out.col(k) = arma::sin(2.0 * pi * (j + 1.0) * X / s);
    k = k + 1;
  }
  return(out);
}

arma::mat getInv(arma::mat A) {
  int nn = A.n_cols;
  arma::mat out(nn, nn);
  bool success = true;
  
  success = arma::inv_sympd(out, A);
  if (success == false) {
    success = arma::inv(out, A);
    if (success == false) {
      success = arma::pinv(out, A);
    }
  } 
  
  return(out);
}

/////////////////////////////

arma::mat removeRow(arma::mat A, int a) {
  
  arma::mat out(A.n_rows - 1, A.n_cols);
  int n = A.n_rows;
  int j = 0;
  for (int i = 0; i < n; i++) {
    if (i != a) {
      out.row(j) = A.row(i);
      j = j + 1;
    }
  }
  return(out);
} 


arma::mat removeCol(arma::mat A, int a) {
  
  arma::mat out(A.n_rows, A.n_cols - 1);
  int n = A.n_cols;
  int j = 0;
  for (int i = 0; i < n; i++) {
    if (i != a) {
      out.col(j) = A.col(i);
      j = j + 1;
    }
  }
  return(out);
} 


arma::mat removeCols(arma::mat A, arma::uvec a) {
  
  int mm = a.n_elem;
  
  arma::mat out(A.n_rows, A.n_cols - mm);
  int n = A.n_cols;
  int j = 0;
  int k = 0;
  int i = 0;
  double tmp;
  for (i = 0; i < n; i++) {
    //Rcpp::Rcout << k << std::endl;
    //Rcpp::Rcout << a(k) << std::endl;
    tmp = i + 0.0;
    if (tmp != a(k)) {
      out.col(j) = A.col(i);
      j = j + 1;
    } else {
      if (k + 1 < mm - 1) {
        k = k + 1;
      } else {
        k = mm - 1;
      }
    }
  }
  return(out);
} 


arma::mat removeRows(arma::mat A, arma::uvec a) {
  
  int mm = a.n_elem;
  
  arma::mat out(A.n_rows - mm, A.n_cols);
  int n = A.n_rows;
  int j = 0;
  int k = 0;
  int i = 0;
  double tmp;
  for (i = 0; i < n; i++) {
    tmp = i + 0.0;
    if (tmp != a(k)) {
      out.row(j) = A.row(i);
      j = j + 1;
    } else {
      if (k + 1 < mm - 1) {
        k = k + 1;
      } else {
        k = mm - 1;
      }
    }
  }
  return(out);
} 


Rcpp::List arimacpp(arma::colvec Y, int q){
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats");
  
  // Picking up Matrix() function from Matrix package
  Rcpp::Function arima = pkg["arima"];
  
  Rcpp::NumericVector ord = Rcpp::NumericVector::create(q, 0, 0);
  
  return arima(Y, Rcpp::Named("order") = ord, Rcpp::Named("method") = "CSS");
}

arma::mat checkSym(arma::mat S) {
  arma::mat out;
  if (S.is_symmetric() == false) {
    S = (S.t() + S) / 2;
  }
  out = S;
  return(out);
}

arma::colvec updatePhi(arma::mat V, arma::mat Vas, 
                       arma::mat A, arma::colvec oldPhi, 
                       double sigma2, arma::mat inveta2mat, 
                       double bound0, double boundqplus1,
                       int MonoFlg, Rcpp::String method) {
  
  //int n = V.n_rows;
  int q = Vas.n_cols;
  
  // Initialize 
  arma::mat tVasVas(q, q);
  arma::mat invtVasVas(q, q);
  arma::mat Phihat(q, 1);
  arma::mat Phi = oldPhi;
  arma::mat S(q, q);
  arma::mat tmpS(q, q);
  arma::colvec M(q);
  
  arma::mat Sgg(1, q);
  arma::mat Snotgg(q - 1, q);
  arma::mat Snotnot(q - 1, q - 1);
  arma::mat invSnotgg(q - 1, q - 1);
  
  arma::mat tmpMat(1, 1); 
  arma::mat tmpSS(q, q);
  int gg;
  
  double mj;
  double sj;
  
  double bound1;
  double bound2;
  
  arma::colvec tmp; 
  
  //Rcpp::Rcout << Vas << std::endl;
  
  // Get tV V
  tVasVas = Vas.t() * Vas;
  //Rcpp::Rcout << tVasVas << std::endl;
  
  // Get Phi hat
  invtVasVas = getInv(tVasVas);
  Phihat = invtVasVas* Vas.t() * V;
  
  // update Phi
  if (MonoFlg == 0) {
    if (method == "MT") {
      tmpSS = tVasVas / sigma2 + A;
      tmpSS = checkSym(tmpSS);
      S = getInv(tmpSS);
      S = checkSym(S);
      M = S * ((tVasVas / sigma2) * Phihat);
    } else if (method == "regression") {
      tmpSS = tVasVas + A;
      tmpSS = checkSym(tmpSS);
      tmpS = getInv(tmpSS);
      S = tmpS * sigma2;
      S = checkSym(S);
      M = tmpS * (tVasVas * Phihat);
      
    } else if ((method == "LASSO") || (method == "ALASSO")) {
      tmpSS = tVasVas + inveta2mat;
      tmpSS = checkSym(tmpSS);
      tmpS = getInv(tmpSS);
      S = tmpS * sigma2;
      S = checkSym(S);
      M = tmpS * (tVasVas * Phihat);
    }
    
    Phi = rmvnorm(M, S);
  } else {
    if ((method == "MonoLASSO") || (method == "MonoALASSO")) {
      tmpSS = tVasVas + inveta2mat;
      tmpSS = checkSym(tmpSS);
      tmpS = getInv(tmpSS);
      S = tmpS * sigma2;
      S = checkSym(S);
      M = tmpS * (tVasVas * Phihat);
      for (gg = 0; gg < q; gg++) {
        Sgg = S.row(gg);
        Snotgg = removeRow(S, gg);
        Snotnot = removeCol(Snotgg, gg);
        invSnotgg = getInv(Snotnot);
        
        tmpMat = M(gg) + removeCol(Sgg, gg) * invSnotgg * (removeRow(Phi, gg) - removeRow(M, gg));
        mj = tmpMat(0);
        
        tmpMat = Sgg(gg) - removeCol(Sgg, gg) * invSnotgg * Snotgg.col(gg);
        sj = tmpMat(0);
        
        if (gg == 0) {
          bound1 = abs(bound0);
          bound2 = abs(Phi(1));
        } else if (gg == q - 1) {
          bound1 = abs(Phi(q - 2));
          bound2 = abs(boundqplus1);
        } else {
          bound1 = abs(Phi(gg - 1));
          bound2 = abs(Phi(gg + 1));
        }
        tmp = rtwosegnorm(1, bound2, bound1, mj, sqrt(sj));
        Phi(gg) = tmp(0);
      }
    }
  }
  return(Phi);
}


double updateSigma2(arma::mat resi, arma::colvec Phi, arma::mat inveta2mat, int T, int q, 
                    arma::mat A, double a, double b, Rcpp::String method) {
  double sigma2 = 0.0;
  
  arma::mat tResiResi = resi.t() * resi;
  double tmptResiResi = tResiResi(0);
  
  arma::mat tPhiVarPhi; 
  double tmptPhiVarPhi;
  // update sigma2
  if (method == "MT") {
    sigma2 = 1.0 / R::rgamma((T - q) / 2.0 + a, 1.0 / (tmptResiResi / 2.0 + b));
  } else if (method == "regression") {
    tPhiVarPhi = Phi.t() * getInv(A) * Phi;
    tmptPhiVarPhi = tPhiVarPhi(0);
    sigma2 = 1.0 / R::rgamma(T / 2.0 + a, 1.0 / (tmptResiResi / 2.0 + tmptPhiVarPhi / 2.0 + b));
  } else if ((method == "LASSO") || (method == "ALASSO") || (method == "MonoLASSO") || (method == "MonoALASSO")) {
    tPhiVarPhi = Phi.t() * inveta2mat * Phi;
    tmptPhiVarPhi = tPhiVarPhi(0);
    sigma2 = 1.0 / R::rgamma(T / 2.0 + a, 1.0 / (tmptResiResi / 2.0 + tmptPhiVarPhi / 2.0 + b));
  }
  
  return(sigma2);
}

arma::mat updateinveta2(arma::colvec Phi, double sigma2, arma::mat lambda2, int q, double tol) {
  arma::mat Phim = arma::conv_to<arma::mat>::from(Phi); 
  arma::mat inveta2(q, 1);
  arma::mat muPrime = arma::sqrt(sigma2 * (lambda2 % arma::pow(Phim, -2)));
  arma::colvec tmp; 
  int gg;
  
  double tmpmuPrime;
  double tmplambda2;
  
  for (gg = 0; gg < q; gg++) {
    tmpmuPrime = muPrime(gg, 0);
    tmplambda2 = lambda2(gg, 0);
    //tmp = rrinvgauss(1, tmpmuPrime, tmplambda2);
    if ((-tol < Phim(gg, 0)) && (Phim(gg, 0) < tol)) {
      inveta2(gg, 0) = 1 / R::rgamma(1.0 / 2.0, 2.0 / tmplambda2);
    } else {
      tmp = rinvgaussiancpp(1, tmpmuPrime, tmplambda2);
      inveta2(gg, 0) = tmp(0);
    }
    
    if (inveta2(gg, 0) < tol) {
      inveta2(gg, 0) = tol;
    }
    
  }
  return(inveta2);
}

arma::mat updatelambda2(arma::mat eta2, int q, double alpha, double beta, Rcpp::String method){
  arma::mat lambda2(q, 1); 
  double shape;
  double scale;
  arma::vec tmpsum; 
  int gg;
  // update lambda2
  if ((method == "LASSO") || (method == "MonoLASSO")) {
    shape = alpha + q;
    tmpsum = arma::sum(eta2);
    scale = beta + tmpsum(0)/2.0;
    lambda2.fill(R::rgamma(shape, scale));
  } else if ((method == "ALASSO") || (method == "MonoALASSO")) {
    shape = alpha + 1;
    for (gg = 0; gg < q; gg++) {
      scale = beta + eta2(gg)/2.0;
      lambda2(gg) = R::rgamma(shape, scale);
    }
  }
  return(lambda2);
}

Rcpp::List updateTauGamma(arma::colvec Y, arma::colvec Phi, arma::mat Tau, arma::mat Gamma, 
                          double muq, double sigma2, double pho, double xi2,
                          int T, int q, arma::mat D, arma::mat H_, int Hflg, int m){
  
  arma::mat Ht(T, 1); 
  arma::mat DHt(T, m);
  arma::mat tHtDHt(m, 1); 
  arma::mat Hnot(T, m - 1);
  arma::mat tmpDHt;
  arma::mat tmptHtDHt; 
  
  arma::mat Taunot(m - 1, 1);
  arma::mat Taut(1, 1);
  arma::mat Gammanot(m - 1, 1);
  arma::mat Gammat(1, 1);
  
  arma::mat zetanot(T, 1);
  arma::mat zetat(T, 1);
  
  arma::mat pvec(m, 1); 
  
  double tmpzetanot;
  double tmpzetat;
  arma::mat tmp;
  double p;
  
  arma::mat Gammathat(m, 1);
  double st;
  double mt; 
  
  int jj;
  
  if (Hflg == 1) {
    
    // get DHt and tHtDHt
    
    //for (jj = 0; jj < m; jj++){
    //  Ht = H_.col(jj);
    //  DHt.col(jj) = D * Ht;
    //  tmp = Ht.t() * DHt.col(jj);
    //  tHtDHt(jj) = tmp(0);
    //}
    
    // update Tau and Gamma
    for (jj = 0; jj < m; jj++) {
      Hnot = removeCol(H_, jj);
      Ht = H_.col(jj);
      
      //tmpDHt = DHt.col(jj);
      //tmptHtDHt = tHtDHt(jj);
      
      //############
      Taunot = removeRow(Tau, jj);
      Taut = Tau.row(jj);
      
      Gammanot = removeRow(Gamma, jj);
      Gammat = Gamma.row(jj);
      
      //update Tau
      zetanot = Y - muq - Hnot * (Taunot % Gammanot);
      zetat = zetanot - Ht * Gammat;
      
      tmp = arma::exp(-1.0 / 2.0 / sigma2 * zetanot.t() * D * zetanot);
      tmpzetanot = tmp(0);
      
      tmp = arma::exp(-1.0 / 2.0 / sigma2 * zetat.t() * D * zetat);
      tmpzetat = tmp(0);
      
      p = pho * tmpzetat / (pho * tmpzetat + (1 - pho) * tmpzetanot);
      
      Tau(jj) = R::rbinom(1, p);
      pvec(jj) = p;
      //############
      if (Tau(jj) == 1) {
        
        DHt.col(jj) = D * Ht;
        tmp = Ht.t() * DHt.col(jj);
        tHtDHt(jj) = tmp(0);
        
        tmpDHt = DHt.col(jj);
        tmptHtDHt = tHtDHt(jj);
        
        //#update Gamma
        tmp = tmpDHt.t() * zetanot / tmptHtDHt;
        Gammathat = tmp(0);
        tmp = 1.0 / (tmptHtDHt / sigma2 + 1 / xi2);
        st = tmp(0);
        tmp = st * (tmptHtDHt * Gammathat) / sigma2;
        mt = tmp(0);
      } else {
        mt = 0;
        st = xi2;
      }
      
      Gamma(jj) = R::rnorm(mt, sqrt(st));
    }
    
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["p"] = pvec
  );
  return(out);
}

Rcpp::List updateMuqMu(arma::colvec Y, arma::mat Tau, arma::mat Gamma, double sigma2,
                       arma::mat One, arma::mat D, arma::mat H_, int Hflg, int T, double tol){
  
  arma::colvec zeta;
  arma::mat HTauGamma;
  arma::mat tOneD = One.t() * D;
  arma::mat tmp = tOneD * One;
  double tOneDOne = tmp(0);
  double muqhat;
  double sq;
  double muq;
  arma::mat Mu(T, 1);
  
  if (tOneDOne < tol) {
    tOneDOne = tol;
  }
  
  if (Hflg == 0) {
    zeta = Y;
  } else {
    HTauGamma = H_ * (Tau % Gamma);
    zeta = Y - HTauGamma;
  }
  
  //#cat("tOneDOne:", tOneDOne, "\n")
  tmp = tOneD * zeta / tOneDOne;
  muqhat = tmp(0);
  sq = sigma2 / tOneDOne;
  muq = R::rnorm(muqhat, sqrt(sq));
  
  if (Hflg == 0) {
    Mu.fill(muq);
  } else {
    Mu = muq + HTauGamma;
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["muq"] = muq,
    Rcpp::_["Mu"] = Mu
  ); 
  
  return(out);
}

// [[Rcpp::export]]
Rcpp::List GibbsRFLSMcpp(arma::colvec& Y,int& q, 
                      arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                      double& theta1, double& theta2, double& xi2,
                      Rcpp::String& method, double& bound0, double& boundqplus1,
                      int& nsim, int& by, int& burnin,
                      double& tol, Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  /////////////////////////////////
  arma::mat H_;
  
  // Calculate H
  int Hflg = 1;
  int m = 1;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
  } else {
    Hflg = 0;
  }
  
  int T = Y.n_elem;
  
  // Calculate G
  arma::mat G = getGMat(T, q);
  
  // Initialize ones
  arma::mat One(T, 1);
  One.ones();
  
  // Initialize the output
  arma::mat Phiout(q, nsim);
  Phiout.zeros();
  
  arma::mat sigma2out(1, nsim);
  sigma2out.zeros();
  
  arma::mat Tauout(m, nsim);
  Tauout.zeros();
  
  arma::mat Gammaout(m, nsim);
  Gammaout.zeros();
  
  arma::mat pout(m, nsim);
  pout.zeros();
  
  arma::mat muqout(1, nsim);
  muqout.zeros();
  
  arma::mat Muout(T, nsim);
  Muout.zeros();
  
  arma::mat phoout(1, nsim);
  phoout.zeros();
  
  arma::mat eta2out(q, nsim);
  eta2out.zeros();
  
  arma::mat lambda2out(q, nsim);
  lambda2out.zeros();
  
  // Is it mono?
  int MonoFlg = 0;
  if ((method == "MonoLASSO") || (method == "MonoALASSO")) {
    MonoFlg = 1;
  }
  
  // Initialize the learning
  Rcpp::List model0 = arimacpp(Y, q);
  
  Rcpp::NumericVector coef = model0["coef"];
  
  double muq = coef[q];
  
  arma::mat Phihat(q, 1);
  Rcpp::NumericMatrix varcoef = model0["var.coef"];
  double tmpphi;
  int ii;
  for (ii = 0; ii < q; ii++) {
    tmpphi = coef[ii];
    Phihat(ii) = tmpphi;
  }
  
  arma::mat Phi = Phihat;
  
  
  double bound1;
  double bound2;
  
  double tmpPhi;
  double tmpPhiVar;
  
  arma::colvec tmp;
  
  int gg;
  if (MonoFlg == 1) {
    for (gg = 0; gg < q; gg++) {
      if (gg == 0) {
        bound1 = abs(bound0);
        bound2 = abs(Phi(1));
      } else if (gg == q - 1) {
        bound1 = abs(Phi(q - 2));
        bound2 = abs(boundqplus1);
      } else {
        bound1 = abs(Phi(gg - 1));
        bound2 = abs(Phi(gg + 1));
      }
      tmpPhi = Phi(gg);
      tmpPhiVar = varcoef(gg, gg);
      if (!((bound2 <= abs(tmpPhi)) && 
          (abs(tmpPhi) <= bound1))) {
        tmp = rtwosegnorm(1, boundqplus1, bound1, 
                          tmpPhi, sqrt(tmpPhiVar));
        Phi(gg) = tmp(0);
      }
    }
  }
  
  arma::mat Mu(T, 1);
  Mu.fill(muq);
  double sigma2 = model0["sigma2"];
  
  arma::mat Tau(m, 1);
  Tau.zeros();
  
  arma::mat Gamma(m, 1);
  Gamma.zeros();
  
  arma::mat pvec(m, 1);
  pvec.zeros();
  
  double pho = R::rbeta(theta1, theta2);
  
  arma::mat eta2(q, 1);
  eta2.zeros();
  
  arma::mat lambda2(q, 1);
  lambda2.zeros();
  
  arma::mat inveta2(q, 1);
  
  eta2 = arma::pow(Phi, 2);
  inveta2 = arma::pow(eta2, -1);
  
  arma::mat inveta2mat(q, q);
  inveta2mat.diag() = inveta2;
  
  if ((method == "LASSO") || (method == "MonoLASSO")) {
    lambda2.fill(pow(q * sqrt(sigma2) / arma::accu(arma::abs(Phi)), 2));
  } else if ((method == "ALASSO") || (method == "MonoALASSO")) {
    for (gg = 0; gg < q; gg++) {
      lambda2(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
    }
  }
  
  arma::mat DHt(T, T);
  DHt.zeros();
  arma::mat tHtDHt(T, 1);
  tHtDHt.zeros();
  
  int rr = 0;
  
  int TotalSim = nsim * by + burnin;
  
  //outputseq <- seq(burnin + 1, TotalSim, step)
  
  arma::mat V_(T, 1);
  V_.zeros();
  
  arma::mat V(T - q, 1);
  V.zeros();
  
  arma::mat Vas_(T, q);
  Vas_.zeros();
  arma::mat Vas(T - q, q);
  Vas.zeros();
  
  
  arma::mat VasPhi(T - q, 1);
  arma::mat resi(T - q, 1);
  
  
  arma::mat PhiMat(T - q, T);
  arma::mat C;
  arma::mat D;
  
  Rcpp::List TauGamma; 
  Rcpp::List MuqMu; 
  
  arma::mat tmpSumTau; 
  
  for (ii = 0; ii < TotalSim; ii++) {
    
    if (ii % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((ii + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
    //update V
    V_ = Y - Mu;
    V = V_.rows(q, T - 1);
    Vas_ = getV(V_, q);
    Vas = Vas_.rows(q, T - 1);
    
    //Rcpp::Rcout << Mu << std::endl;
    //Rcpp::Rcout << V << std::endl;
    
    // update Phi
    Phi = updatePhi(V, Vas, A, 
                    Phi, sigma2, inveta2mat, 
                    bound0, boundqplus1,
                    MonoFlg, method);
    
    // Get residuals
    VasPhi = Vas * Phi;
    resi = V - VasPhi;
    
    // update sigma2
    sigma2 = updateSigma2(resi, Phi, inveta2mat, T, q, 
                          A, a, b, method);
    
    // update eta2
    inveta2 = updateinveta2(Phi, sigma2, lambda2, q, tol);
    eta2 = arma::pow(inveta2, -1);
    inveta2mat.diag() = inveta2;
    
    // update lambda2
    lambda2 = updatelambda2(eta2, q, alpha, beta, method);
    
    ///////////////////////////////////////////////////
    //update the random level shift model
    ///////////////////////////////////////////////////
    
    //Calculate Phi Matrix 
    PhiMat = getPhiMat(Phi, T);
    
    //Calculate C Matrix
    C = G - PhiMat;
    
    //Calculate D Matrix
    D = C.t() * C;
    
    //#update Tau and Gamma
    
    TauGamma = updateTauGamma(Y, Phi, Tau, Gamma, 
                              muq, sigma2, pho, xi2,
                              T, q, D, H_, Hflg, m);
    
    Tau = Rcpp::as<arma::mat>(TauGamma["Tau"]);
    Gamma = Rcpp::as<arma::mat>(TauGamma["Gamma"]);
    pvec = Rcpp::as<arma::mat>(TauGamma["p"]);
    
    //#update muq and Mu
    
    MuqMu = updateMuqMu(Y, Tau, Gamma, sigma2,
                        One, D, H_, Hflg, T, tol);
    
    muq = MuqMu["muq"];
    //Rcpp::Rcout << muq << std::endl;
    
    Mu = Rcpp::as<arma::mat>(MuqMu["Mu"]);
    
    //#update pho
    if (Hflg == 1) {
      tmpSumTau = arma::sum(Tau);
      pho = R::rbeta(theta1 + tmpSumTau(0), theta2 + m - tmpSumTau(0));
    }
    
    if (ii >= burnin) {
      if (ii % by == 0) {
        Phiout.col(rr) = Phi;
        sigma2out.col(rr) = sigma2;
        Tauout.col(rr) = Tau;
        Gammaout.col(rr) = Gamma;
        pout.col(rr) = pvec;
        muqout.col(rr) = muq;
        Muout.col(rr) = Mu;
        phoout.col(rr) = pho;
        eta2out.col(rr) = eta2;
        lambda2out.col(rr) = lambda2;
        rr = rr + 1;
      }
    }
    
  }
  
  /////////////////////////////////
  
  Rcpp::Rcout <<"Training: 100%" << std::endl;
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  
  /////////////////////////////////
  
  Rcpp::List out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["p"] = pout,
    _["muq"] = muqout,
    _["Mu"] = Muout,
    //_["pho"] = phoout,
    //_["eta2"] = eta2out,
    _["lambda2"] = lambda2out
  );
  
  return(out);
  
}
 