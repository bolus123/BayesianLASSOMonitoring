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

arma::colvec rmvnorm(arma::colvec Mean,arma::mat Sigma) {
  
  int q = Mean.n_elem;
 arma::colvec Z =arma::randn(q);
 arma::colvec out = Mean +arma::chol(Sigma) * Z;
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
  arma::colvec U =arma::randu(n);
   
   double a1 = (-a - mean) / sd;
   double b1 = (-b - mean) / sd;
   double a2 = (a - mean) / sd;
   double b2 = (b - mean) / sd;
   
   double p1 =arma::normcdf(a1) -arma::normcdf(b1);
   double p2 =arma::normcdf(b2) -arma::normcdf(a2);
   double p = p1 + p2;
   
   int i;
   if (p > 0) {
     for (i = 0; i < n; i++) {
       if (U(i) <= p1 / p) {
         out(i) = R::qnorm5(U(i) * p +arma::normcdf(b1), 0.0, 1.0, 1, 0) * sd + mean;
       } else {
         out(i) = R::qnorm5(U(i) * p +arma::normcdf(a2) - p1, 0.0, 1.0, 1, 0) * sd + mean;
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
  arma::mat L =arma::trimatl(tmp);
   
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
  arma::mat L =arma::trimatl(tmp);
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
     out.col(k) =arma::cos(2.0 * pi * (j + 1.0) * X / s);
     k = k + 1;
     
     out.col(k) =arma::sin(2.0 * pi * (j + 1.0) * X / s);
     k = k + 1;
   }
   return(out);
 }

arma::mat getInv(arma::mat A) {
  int nn = A.n_cols;
 arma::mat out(nn, nn);
  bool success = true;
  
  success =arma::inv_sympd(out, A);
  if (success == false) {
    success =arma::inv(out, A);
    if (success == false) {
      success =arma::pinv(out, A);
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


arma::mat removeCols(arma::mat A,arma::uvec a) {
  
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


arma::mat removeRows(arma::mat A,arma::uvec a) {
  
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

Rcpp::List arimaxcpp(arma::colvec Y, int q, arma::mat X) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats");
  
  // Picking up Matrix() function from Matrix package
  Rcpp::Function arima = pkg["arima"];
  
  Rcpp::NumericVector ord = Rcpp::NumericVector::create(q, 0, 0);
  
  return arima(Y, Rcpp::Named("order") = ord, Rcpp::Named("xreg") = X, Rcpp::Named("method") = "CSS");
}

arma::mat checkSym(arma::mat S) {
 arma::mat out;
  if (S.is_symmetric() == false) {
    S = (S.t() + S) / 2;
  }
  out = S;
  return(out);
}

arma::colvec updatePhi(arma::mat V,arma::mat Vas, 
                      arma::mat A,arma::colvec oldPhi, 
                       double sigma2,arma::mat inveta2mat, 
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


double updateSigma2(arma::mat resi,arma::colvec Phi,arma::mat inveta2mat, int T, int q, 
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

arma::mat updateinveta2(arma::colvec Phi, double sigma2,arma::mat lambda2, int q, double tol) {
 arma::mat Phim =arma::conv_to<arma::mat>::from(Phi); 
 arma::mat inveta2(q, 1);
 arma::mat muPrime =arma::sqrt(sigma2 * (lambda2 %arma::pow(Phim, -2)));
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
    tmpsum =arma::sum(eta2);
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

Rcpp::List updateTauGamma(arma::colvec Y,arma::colvec Phi,arma::mat Tau,arma::mat Gamma, 
                          double muq, double sigma2, double pho, double xi2,
                          int T, int q,arma::mat D,arma::mat H_, int Hflg, int m){
  
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
  
 arma::mat muGamma = Gamma;
 arma::mat sigma2Gamma = Gamma;
  
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
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetanot.t() * D * zetanot);
      tmpzetanot = tmp(0);
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetat.t() * D * zetat);
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
      muGamma(jj) = mt;
      sigma2Gamma(jj) = st;
    }
    
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["p"] = pvec,
    Rcpp::_["muGamma"] = muGamma,
    Rcpp::_["sigma2Gamma"] = sigma2Gamma
  );
  return(out);
}

Rcpp::List updateMuqMu(arma::colvec Y,arma::mat Tau,arma::mat Gamma, double sigma2,
                      arma::mat One,arma::mat D,arma::mat H_, int Hflg, int T, double tol){
  
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
  
 arma::mat muGammaout(m, nsim);
  muGammaout.zeros();
 arma::mat sigma2Gammaout(m, nsim);
  sigma2Gammaout.zeros();
  
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
  
 arma::mat muGamma(m, 1);
  muGamma.zeros();
 arma::mat sigma2Gamma(m, 1);
  sigma2Gamma.zeros();
  
 arma::mat pvec(m, 1);
  pvec.zeros();
  
  double pho = R::rbeta(theta1, theta2);
  
 arma::mat eta2(q, 1);
  eta2.zeros();
  
 arma::mat lambda2(q, 1);
  lambda2.zeros();
  
 arma::mat inveta2(q, 1);
  
  eta2 =arma::pow(Phi, 2);
  inveta2 =arma::pow(eta2, -1);
  
 arma::mat inveta2mat(q, q);
  inveta2mat.diag() = inveta2;
  
  if ((method == "LASSO") || (method == "MonoLASSO")) {
    lambda2.fill(pow(q * sqrt(sigma2) /arma::accu(arma::abs(Phi)), 2));
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
    eta2 =arma::pow(inveta2, -1);
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
    muGamma = Rcpp::as<arma::mat>(TauGamma["muGamma"]);
    sigma2Gamma = Rcpp::as<arma::mat>(TauGamma["sigma2Gamma"]);
    
    //#update muq and Mu
    
    MuqMu = updateMuqMu(Y, Tau, Gamma, sigma2,
                        One, D, H_, Hflg, T, tol);
    
    muq = MuqMu["muq"];
    //Rcpp::Rcout << muq << std::endl;
    
    Mu = Rcpp::as<arma::mat>(MuqMu["Mu"]);
    
    //#update pho
    if (Hflg == 1) {
      tmpSumTau =arma::sum(Tau);
      pho = R::rbeta(theta1 + tmpSumTau(0), theta2 + m - tmpSumTau(0));
    }
    
    if (ii >= burnin) {
      if (ii % by == 0) {
        Phiout.col(rr) = Phi;
        sigma2out.col(rr) = sigma2;
        Tauout.col(rr) = Tau;
        Gammaout.col(rr) = Gamma;
        muGammaout.col(rr) = muGamma;
        sigma2Gammaout.col(rr) = sigma2Gamma;
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
    _["muGamma"] = muGammaout,
    _["sigma2Gamma"] = sigma2Gammaout,
    _["p"] = pout,
    _["muq"] = muqout,
    _["Mu"] = Muout,
    //_["pho"] = phoout,
    //_["eta2"] = eta2out,
    _["lambda2"] = lambda2out
  );
  
  return(out);
  
}


// [[Rcpp::export]]
Rcpp::List GibbsRFLSMUpdatecpp(arma::colvec Y,int q, 
                        arma::mat A, double a, double b, double alpha, double beta, 
                         double theta1, double theta2, double xi2,
                         Rcpp::String method, double bound0, double boundqplus1,
                         int nsim, int by, int burnin,
                         double tol, 
                         Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                         Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  //auto start = std::chrono::system_clock::now();
  //std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  //Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
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
 arma::mat G_;
  
  if (G.isNotNull()) {
    G_ = Rcpp::as<arma::mat>(G);
  } else {
    G_ = getGMat(T, q);
  }
  
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
  
 arma::mat muGammaout(m, nsim);
  muGammaout.zeros();
 arma::mat sigma2Gammaout(m, nsim);
  sigma2Gammaout.zeros();
  
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
  
  Rcpp::List model0;
  Rcpp::NumericVector coef;
  double muq;
 arma::mat Phihat(q, 1);
  Rcpp::NumericMatrix varcoef;
  double tmpphi;
  int ii;
 arma::mat Phi;
  double bound1;
  double bound2;
  
  double tmpPhi;
  double tmpPhiVar;
  
 arma::colvec tmp;
  int gg;
  
 arma::mat Mu(T, 1);
  double sigma2;
  
 arma::mat Tau(m, 1);
  Tau.zeros();
  
 arma::mat Gamma(m, 1);
  Gamma.zeros();
  
 arma::mat muGamma(m, 1);
  muGamma.zeros();
 arma::mat sigma2Gamma(m, 1);
  sigma2Gamma.zeros();
  
 arma::mat pvec(m, 1);
  pvec.zeros();
  
  double pho;
  
 arma::mat eta2(q, 1);
  eta2.zeros();
  
 arma::mat lambda2(q, 1);
  lambda2.zeros();
  
 arma::mat inveta2(q, 1);
 arma::mat inveta2mat(q, q);
  
  Rcpp::List oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
  
  //Rcpp::Rcout << 1 << std::endl;
  
  if (oldpars.isNotNull()) {
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
    Tau = Rcpp::as<arma::mat>(oldpars_["Tau"]);
    Gamma = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
    pho = oldpars_["pho"];
    
    //Rcpp::Rcout << 2 << std::endl;
    
    if ((method == "LASSO") || (method == "ALASSO") || (method == "MonoLASSO") || (method == "MonoALASSO")) {
      eta2 = Rcpp::as<arma::mat>(oldpars_["eta2"]);
      inveta2 =arma::pow(eta2, -1);
      inveta2mat.diag() = inveta2;
      lambda2 = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
    }
    
    muq = oldpars_["muq"];
  } else {
    model0 = arimacpp(Y, q);
    
    coef = model0["coef"];
    
    muq = coef[q];
    
    varcoef = Rcpp::as<Rcpp::NumericMatrix>(model0["var.coef"]);
    for (ii = 0; ii < q; ii++) {
      tmpphi = coef[ii];
      Phihat(ii) = tmpphi;
    }
    
    Phi = Phihat;
    
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
    
    Mu.fill(muq);
    sigma2 = model0["sigma2"];
    
    pho = R::rbeta(theta1, theta2);
    
    eta2 =arma::pow(Phi, 2);
    inveta2 =arma::pow(eta2, -1);
    
    inveta2mat.diag() = inveta2;
    
    if ((method == "LASSO") || (method == "MonoLASSO")) {
      lambda2.fill(pow(q * sqrt(sigma2) /arma::accu(arma::abs(Phi)), 2));
    } else if ((method == "ALASSO") || (method == "MonoALASSO")) {
      for (gg = 0; gg < q; gg++) {
        lambda2(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
      }
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
    
    //if (ii % 100 == 0) {
   //  Rcpp::Rcout <<"Training: " << ((ii + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    //}
    
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
    eta2 =arma::pow(inveta2, -1);
    inveta2mat.diag() = inveta2;
    
   // update lambda2
    lambda2 = updatelambda2(eta2, q, alpha, beta, method);
    
    ///////////////////////////////////////////////////
    //update the random level shift model
    ///////////////////////////////////////////////////
    
    //Calculate Phi Matrix 
    PhiMat = getPhiMat(Phi, T);
    
    //Calculate C Matrix
    C = G_ - PhiMat;
    
    //Calculate D Matrix
    D = C.t() * C;
    
    //#update Tau and Gamma
    
    TauGamma = updateTauGamma(Y, Phi, Tau, Gamma, 
                              muq, sigma2, pho, xi2,
                              T, q, D, H_, Hflg, m);
    
    Tau = Rcpp::as<arma::mat>(TauGamma["Tau"]);
    Gamma = Rcpp::as<arma::mat>(TauGamma["Gamma"]);
    pvec = Rcpp::as<arma::mat>(TauGamma["p"]);
    muGamma = Rcpp::as<arma::mat>(TauGamma["muGamma"]);
    sigma2Gamma = Rcpp::as<arma::mat>(TauGamma["sigma2Gamma"]);
    
    //#update muq and Mu
    
    MuqMu = updateMuqMu(Y, Tau, Gamma, sigma2,
                        One, D, H_, Hflg, T, tol);
    
    muq = MuqMu["muq"];
    //Rcpp::Rcout << muq << std::endl;
    
    Mu = Rcpp::as<arma::mat>(MuqMu["Mu"]);
    
    //#update pho
    if (Hflg == 1) {
      tmpSumTau =arma::sum(Tau);
      pho = R::rbeta(theta1 + tmpSumTau(0), theta2 + m - tmpSumTau(0));
    }
    
    if (ii >= burnin) {
      if (ii % by == 0) {
        Phiout.col(rr) = Phi;
        sigma2out.col(rr) = sigma2;
        Tauout.col(rr) = Tau;
        Gammaout.col(rr) = Gamma;
        muGammaout.col(rr) = muGamma;
        sigma2Gammaout.col(rr) = sigma2Gamma;
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
  
  //Rcpp::Rcout <<"Training: 100%" << std::endl;
  //
  //auto end = std::chrono::system_clock::now();
  //std::chrono::duration<double> elapsed_seconds = end-start;
  //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  //
  //Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
 //            << "Elapsed time: " << elapsed_seconds.count() << "s"
 //            << std::endl;
  
  /////////////////////////////////
  
  Rcpp::List out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    //_["muGamma"] = muGammaout,
    //_["sigma2Gamma"] = sigma2Gammaout,
    //_["p"] = pout,
    _["muq"] = muqout,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out
  );
  
  return(out);
  
}


//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::mat lhf(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2) {
  
  double pi = 3.14159265359;
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Y.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V_;
 arma::mat V; 
 arma::mat Vas_;
 arma::mat Vas;
 arma::mat VasPhi;
 arma::mat resi; 
  
  V_ = Y - Mu;
  V = V_.rows(q, T - 1);
  Vas_ = getV(V_, q);
  Vas = Vas_.rows(q, T - 1);
  
  //Rcpp::Rcout << V << std::endl;
  //Rcpp::Rcout << Vas << std::endl;
  
  VasPhi = Vas * Phi;
  
  //Rcpp::Rcout << VasPhi << std::endl;
  
  resi = V - VasPhi;
  
  //Rcpp::Rcout << resi << std::endl;
  
 arma::mat lh = sqrt(1.0 / 2.0 / pi / sigma2) * exp(- 1.0 / 2.0 *arma::pow(resi, 2.0) / sigma2);
  return(lh);
  
}

// [[Rcpp::export]]
arma::colvec boxcoxtr(arma::colvec Y, double theta) {
  
 arma::colvec Ybc = (arma::pow(Y, theta) - 1.0) / theta;
  return(Ybc);
  
}

// [[Rcpp::export]]
arma::colvec invboxcoxtr(arma::colvec Ybc, double theta) {
  
 arma::colvec Y =arma::pow(Ybc * theta + 1.0, 1.0 / theta);
  return(Y);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec invyeojohnsontr(arma::colvec Yyj, double theta, double eps) {
  int T = Yyj.n_elem;
 arma::colvec ans(T);
  
  for (int i = 0; i < T; i++) {
    if ((Yyj(i) >= 0.0) && (abs(theta) > eps)) {
      ans(i) = pow(Yyj(i) * theta + 1.0, 1.0 / theta) - 1.0;
    }
    else if ((Yyj(i) >= 0.0) && (abs(theta) <= eps)) {
      ans(i) = exp(Yyj(i));
    }
    else if ((Yyj(i) < 0.0) && (abs(theta - 2.0) > eps)) {
      ans(i) = 1.0 - pow(-(2.0 - theta) * Yyj(i) + 1.0, 1.0 / (2.0 - theta));
    }
    else if ((Yyj(i) < 0.0) && (abs(theta - 2.0) <= eps)) {
      ans(i) = -exp(-Yyj(i));
    }
  }
  return(ans);
} 

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec yeojohnsontr(arma::colvec Y, double theta, double eps) {
  int T = Y.n_elem;
 arma::colvec ans(T);
  
  for (int i = 0; i < T; i++) {
    if ((Y(i) >= 0) && (abs(theta) > eps)) {
      ans(i) = (pow(Y(i) + 1, theta) - 1) / theta;
    }
    else if ((Y(i) >= 0) && (abs(theta) <= eps)) {
      ans(i) = log(Y(i));
    }
    else if ((Y(i) < 0) && (abs(theta - 2) > eps)) {
      ans(i) = -((pow(-Y(i) + 1, 2 - theta) - 1) / (2 - theta));
    }
    else if ((Y(i) < 0) && (abs(theta - 2) <= eps)) {
      ans(i) = -log(-Y(i));
    }
  }
  
  return(ans);
} 


// [[Rcpp::export]]
arma::mat lhBCf(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
 arma::colvec Ybc = boxcoxtr(Y, theta);
 arma::mat lh = lhf(Ybc, Phi, Mu, sigma2);
 arma::mat lhBC = lh %arma::pow(Y.rows(q, T - 1), theta - 1.0);
  
  return(lhBC);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::mat lhYJf(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
 arma::colvec Yyj = yeojohnsontr(Y, theta, 0.000001);
 arma::mat lh = lhf(Yyj, Phi, Mu, sigma2);
  
 arma::mat lhYJ = lh;
  
  for (int i = 0; i < (T - q); i++) {
    lhYJ(i) = lhYJ(i) * pow(abs(Y(i)) + 1, (theta - 1) * sign(Y(i)));
  }
  
  return(lhYJ);
  
}



// [[Rcpp::export]]
arma::mat thetaBoxCoxMH(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                  double oldtheta, int burnin, int nsim, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double oldtheta_ = oldtheta;
  double thetaas;
  double A;
  
 arma::mat thetaout(nsim, 1); 
  
 arma::mat oldlhBC = lhBCf(Y, Phi, Mu, sigma2, oldtheta_);
  //double sumoldllhBC =arma::accu(oldllhBC);
 arma::uvec ind0 =arma::find(oldlhBC <= tol); 
  oldlhBC(ind0).fill(tol);
  
 arma::uvec ind1; 
  
 arma::mat newlhBC; 
  //double sumnewllhBC;
  
 arma::mat tmp; 
  double pd;
  
  int i;
  int j = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    u = R::runif(0.0, 1.0);
    thetaas = R::rnorm(oldtheta_, 0.1);
    newlhBC = lhBCf(Y, Phi, Mu, sigma2, thetaas);
    ind1 =arma::find(newlhBC <= tol);
    newlhBC(ind1).fill(tol);
    
    tmp =arma::cumprod(newlhBC %arma::pow(oldlhBC, -1));
    tmp = tmp * (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * thetaas * thetaas)) /
      (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * oldtheta_ * oldtheta_));
      
    //Rcpp::Rcout << tmp << std::endl;
    pd = tmp(T - q - 1);
    //Rcpp::Rcout << pd << std::endl;
    A = std::min(1.0, pd);
    //Rcpp::Rcout << tmp(T - q - 1) << std::endl;
    //Rcpp::Rcout << A << std::endl;
    
    if (u < A) {
      oldtheta_ = thetaas;
      oldlhBC = newlhBC;
    } 
    
    //Rcpp::Rcout << oldtheta_ << std::endl;
    
    if (i >= burnin) {
      thetaout(j, 0) = oldtheta_;
      j = j + 1;
    }
  }
  
  return(thetaout);
  
}

// [[Rcpp::export]]
arma::mat thetaYeoJohnsonMH(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                        double oldtheta, int burnin, int nsim, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double oldtheta_ = oldtheta;
  double thetaas;
  double A;
  
 arma::mat thetaout(nsim, 1); 
  
 arma::mat oldlhYJ = lhYJf(Y, Phi, Mu, sigma2, oldtheta_);
  //double sumoldllhBC =arma::accu(oldllhBC);
 arma::uvec ind0 =arma::find(oldlhYJ <= tol); 
  oldlhYJ(ind0).fill(tol);
  
 arma::uvec ind1; 
  
 arma::mat newlhYJ; 
  //double sumnewllhBC;
  
 arma::mat tmp; 
  double pd;
  
  int i;
  int j = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    u = R::runif(0.0, 1.0);
    thetaas = R::rnorm(oldtheta_, 0.1);
    newlhYJ = lhYJf(Y, Phi, Mu, sigma2, thetaas);
    ind1 =arma::find(newlhYJ <= tol);
    newlhYJ(ind1).fill(tol);
    
    tmp =arma::cumprod(newlhYJ %arma::pow(oldlhYJ, -1));
    tmp = tmp * (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * thetaas * thetaas)) /
      (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * oldtheta_ * oldtheta_));
    
    //Rcpp::Rcout << tmp << std::endl;
    pd = tmp(T - q - 1);
    //Rcpp::Rcout << pd << std::endl;
    A = std::min(1.0, pd);
    //Rcpp::Rcout << tmp(T - q - 1) << std::endl;
    //Rcpp::Rcout << A << std::endl;
    
    if (u < A) {
      oldtheta_ = thetaas;
      oldlhYJ = newlhYJ;
    } 
    
    //Rcpp::Rcout << oldtheta_ << std::endl;
    
    if (i >= burnin) {
      thetaout(j, 0) = oldtheta_;
      j = j + 1;
    }
  }
  
  return(thetaout);
  
}

// [[Rcpp::export]]
Rcpp::List GibbsRFLSMBoxCoxcpp(arma::colvec& Y,int& q, 
                        arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                         double& theta1, double& theta2, double& xi2,
                         Rcpp::String& method, double& bound0, double& boundqplus1,
                         int updateBC, double& theta,
                         int& nsim, int& by, int& burnin,
                         double& tol, 
                         Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                         Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
 arma::mat H_;
  
 // Calculate H
  int m = 1;
  if (H.isNotNull()) {
    m = H_.n_cols;
  } 
  
  int T = Y.n_elem;
  
  int TotalSim = nsim * by + burnin;
  
  Rcpp::List GibbsRFLSMModel; 
 arma::colvec Ybc;
  
  Rcpp::List oldpars_;
  
 arma::mat Phi(q, 1);
  Phi.zeros();
  
 arma::mat Mu(T - q, 1);
  Mu.zeros();
  
  double sigma2 = 1.0;
  double theta_ = theta;
 arma::mat tmp; 
  
  Ybc = boxcoxtr(Y, theta_);
  
  if (oldpars.isNotNull()) {
    oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  } else {
    
    oldpars_ = GibbsRFLSMUpdatecpp(Ybc, q, 
                             A, a, b, alpha, beta, 
                             theta1, theta2, xi2, 
                             method, bound0, boundqplus1, 
                             1, 1, TotalSim / 10, tol,
                             G, oldpars, H);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
  Rcpp::NumericMatrix GG;
 arma::mat G_;
  
  if (G.isNotNull()) {
    GG = Rcpp::as<Rcpp::NumericMatrix>(G);
  } else {
    G_ = getGMat(T, q);
    GG = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(G_));
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  
  
  //////////////////////////
  
 arma::mat Phiout(q, nsim); 
 arma::mat sigma2out(1, nsim);
 arma::mat Tauout(m, nsim);
 arma::mat Gammaout(m, nsim);
 arma::mat muqout(1, nsim);
 arma::mat Muout(T, nsim);
 arma::mat phoout(1, nsim);
 arma::mat eta2out(q, nsim);
 arma::mat lambda2out(q, nsim);
 arma::mat thetaout(1, nsim);
  
  //////////////////////////
  
  int i;
  int rr = 0;
  
  for (i = 0; i < TotalSim; i++) {
    
    if (i % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((i + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
    if (updateBC == 1) {
      tmp = thetaBoxCoxMH(Y, Phi, Mu, sigma2, 
                    theta_, 1, 1, tol);
      theta_ = tmp(0);
      //Rcpp::Rcout << theta_ << std::endl;
      Ybc = boxcoxtr(Y, theta_);
    }
    
 //  Rcpp::Rcout << 4 << std::endl;
 //  
    oldpars_ = GibbsRFLSMUpdatecpp(Ybc, q, 
                             A, a, b, alpha, beta, 
                             theta1, theta2, xi2, 
                             method, bound0, boundqplus1, 
                             1, 1, 0, tol,
                             GG, oldpars_, H);
 //  
 //  Rcpp::Rcout << 5 << std::endl;
    
    if (i >= burnin) {
      if (i % by == 0) {
        Phiout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Phi"]);
        sigma2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["sigma2"]);
        Tauout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Tau"]);
        Gammaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
        muqout.col(rr) = Rcpp::as<arma::mat>(oldpars_["muq"]);
        Muout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Mu"]);
        phoout.col(rr) = Rcpp::as<arma::mat>(oldpars_["pho"]);
        eta2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["eta2"]);
        lambda2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
        thetaout.col(rr) = theta_;
        rr = rr + 1;
      }
    }
    
  }
  
  Rcpp::List out; 
  out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["muq"] = muqout,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out,
    _["theta"] = thetaout
  );
  
  return(out);
  
}


//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec rtrnorm(int n, double mean, double sd, double lower, double upper) {
 arma::colvec out(n);
 arma::colvec U =arma::randu(n);
  
  double slower = (lower - mean) / sd;
  double supper = (upper - mean) / sd;
  
  double Z = R::pnorm5(supper, 0.0, 1.0, 1, 0) - R::pnorm5(slower, 0.0, 1.0, 1, 0);
  
  int i;
  
  for (i = 0; i < n; i++) {
    out(i) = R::qnorm5(U(i) * Z + R::pnorm5(slower, 0.0, 1.0, 1, 0), 0.0, 1.0, 1, 0) * sd + mean;
  }
  
  return(out);
  
}


arma::colvec getucY(arma::colvec Yyj,arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Y.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::colvec ucY = Y;
 arma::colvec ucYyj = Yyj; 
  
 arma::mat fit(T - q, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (int i = 0; i < T; i++) {
    V(i, 0) = ucYyj(i) - Mu(i, 0);
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i, 0) + VasPhi(0);
      
      if (Y(i) == 0) { 
        //Rcpp::Rcout << 3.01 << std::endl;
        tmp = rtrnorm(1, fit(i - q, 0), sqrt(sigma2), (-1.0) *arma::datum::inf, 0.0);
        //Rcpp::Rcout << 3.1 << std::endl;
        ucYyj(i) = tmp(0); 
        //Rcpp::Rcout << 3.2 << std::endl;
      }
      //Rcpp::Rcout << ucYyj(i)  << std::endl;
    }
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(ucYyj);
  
}



// [[Rcpp::export]]
Rcpp::List GibbsRFLSMYeoJohnsonZcpp(arma::colvec& Y,int& q, 
                              arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                               double& theta1, double& theta2, double& xi2,
                               Rcpp::String& method, double& bound0, double& boundqplus1,
                               int updateYJ, double& theta,
                               int updateZ, double eps,
                               int& nsim, int& by, int& burnin,
                               double& tol, 
                               Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                               Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  ///////////////////////////////////
  
 arma::mat H_;
  
 // Calculate H
  int m = 1;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
  } 
  
  int T = Y.n_elem;
  
  int TotalSim = nsim * by + burnin;
  
  Rcpp::List GibbsRFLSMModel; 
 arma::colvec Yyj = Y;
  
  Rcpp::List oldpars_;
  
 arma::mat Phi(q, 1);
  Phi.zeros();
  
 arma::mat Mu(T - q, 1);
  Mu.zeros();
  
  double sigma2 = 1.0;
  double theta_ = theta;
 arma::mat tmp; 
  
  if (updateYJ == 1) {
    Yyj = yeojohnsontr(Y, theta_, eps);
  }
  
  
  if (oldpars.isNotNull()) {
    oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  } else {
    
    oldpars_ = GibbsRFLSMUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, bound0, boundqplus1, 
                                   1, 1, TotalSim / 10, tol,
                                   G, oldpars, H);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  }
  
  if (updateZ == 1) {
    Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
  }
  
  Rcpp::NumericMatrix GG;
 arma::mat G_;
  
  if (G.isNotNull()) {
    GG = Rcpp::as<Rcpp::NumericMatrix>(G);
  } else {
    G_ = getGMat(T, q);
    GG = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(G_));
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  
  
  //////////////////////////
  
 arma::mat Phiout(q, nsim); 
 arma::mat sigma2out(1, nsim);
 arma::mat Tauout(m, nsim);
 arma::mat Gammaout(m, nsim);
 arma::mat muqout(1, nsim);
 arma::mat Muout(T, nsim);
 arma::mat phoout(1, nsim);
 arma::mat eta2out(q, nsim);
 arma::mat lambda2out(q, nsim);
 arma::mat thetaout(1, nsim);
 arma::mat Yyjout(T, nsim);
  
  //////////////////////////
  
  
  //////////////////////////
  
  int i;
  int rr = 0;
  
  for (i = 0; i < TotalSim; i++) {
    
    if (i % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((i + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
   //  Rcpp::Rcout << 4 << std::endl;
   //  
    oldpars_ = GibbsRFLSMUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   GG, oldpars_, H);
   //  
   //  Rcpp::Rcout << 5 << std::endl;
    
    if (updateYJ == 1) {
      tmp = thetaYeoJohnsonMH(Y, Phi, Mu, sigma2, 
                              theta_, 1, 1, tol);
      theta_ = tmp(0);
      //Rcpp::Rcout << theta_ << std::endl;
      Yyj = yeojohnsontr(Y, theta_, eps);
    }
    
    if (updateZ == 1) {
      Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
    }
    
    if (i >= burnin) {
      if (i % by == 0) {
        Phiout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Phi"]);
        sigma2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["sigma2"]);
        Tauout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Tau"]);
        Gammaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
        muqout.col(rr) = Rcpp::as<arma::mat>(oldpars_["muq"]);
        Muout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Mu"]);
        phoout.col(rr) = Rcpp::as<arma::mat>(oldpars_["pho"]);
        eta2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["eta2"]);
        lambda2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
        thetaout.col(rr) = theta_;
        Yyjout.col(rr) = Yyj;
        rr = rr + 1;
      }
    }
    
  }
  
  Rcpp::Rcout <<"Training: 100%" << std::endl;
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  
  
  Rcpp::List out; 
  out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["muq"] = muqout,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out,
    _["theta"] = thetaout,
    _["Yyj"] = Yyjout
  );
  
  return(out);
  
}


//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec getfityj(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double eps) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyj.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::mat fit(T - q, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (int i = 0; i < T; i++) {
    V(i, 0) = Yyj(i) - Mu(i, 0);
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i, 0) + VasPhi(0);
    }
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(fit);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec getfit(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double theta, double eps) {
  
 arma::colvec fittr = getfityj(Yyj, Phi, Mu, eps);
 arma::colvec fit = invyeojohnsontr(fittr, theta, eps);
  return(fit);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYyjph1(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double sigma2) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyj.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::mat fit(T - q, 1);
 arma::mat simYyjph1(T - q, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (int i = 0; i < T; i++) {
    V(i, 0) = Yyj(i) - Mu(i, 0);
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i, 0) + VasPhi(0);
      simYyjph1(i - q, 0) = R::rnorm(fit(i - q, 0), sqrt(sigma2));
    }
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(simYyjph1);
  
}


//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYph1(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps) {
  
 arma::colvec Yyjph1 = simYyjph1(Yyj, Phi, Mu, sigma2); 
 arma::colvec Yph1 = invyeojohnsontr(Yyjph1, theta, eps); 
 arma::uvec ind0 =arma::find(Yph1 <= 0.0); 
  Yph1(ind0).fill(0.0);
  
  
  return(Yph1);
  
}


//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYyjph2(int h,arma::colvec Yyjph1,arma::mat Phi,arma::mat Mu, double sigma2) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyjph1.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T + h, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::mat fit(T - q + h, 1);
 arma::mat simYyjph2(T - q + h, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  //for (int i = 0; i < (T + h); i++) {
  for (int i = 0; i < (T + h); i++) {
    
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i, 0) + VasPhi(0);
      simYyjph2(i - q, 0) = R::rnorm(fit(i - q, 0), sqrt(sigma2));
    }
    
    if (i < T) {
      V(i, 0) = Yyjph1(i) - Mu(i, 0);
    } else if (i >= T) {
      V(i, 0) = simYyjph2(i - q, 0) - Mu(i, 0);
    }
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(simYyjph2);

}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYph2(int h,arma::colvec Yyjph1,arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps) {
  
 arma::colvec Yyjph2 = simYyjph2(h, Yyjph1, Phi, Mu, sigma2); 
 arma::colvec Yph2 = invyeojohnsontr(Yyjph2, theta, eps); 
 arma::uvec ind0 =arma::find(Yph2 <= 0.0); 
  Yph2(ind0).fill(0.0);
  
  
  return(Yph2);
  
}




arma::colvec updateCoef(arma::mat Y,arma::mat X, 
                        arma::mat A, arma::colvec oldCoef, 
                        double sigma2, arma::mat inveta2mat, 
                        double bound0, double boundqplus1,
                        int mono, Rcpp::String method) {
  
  //int n = X.n_rows;
  int q = X.n_cols;
  
  // Initialize 
  arma::mat tXX(q, q);
  arma::mat invtXX(q, q);
  arma::mat Coefhat(q, 1);
  arma::mat Coef = oldCoef;
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
  tXX = X.t() * X;
  //Rcpp::Rcout << tVasVas << std::endl;
  
  // Get coef hat
  invtXX = getInv(tXX);
  Coefhat = invtXX * X.t() * Y;
  
  // update coef
  if (method == "MT") {
    tmpSS = tXX / sigma2 + A;
    tmpSS = checkSym(tmpSS);
    S = getInv(tmpSS);
    S = checkSym(S);
    M = S * ((tXX / sigma2) * Coefhat);
  } else if (method == "regression") {
    tmpSS = tXX + A;
    tmpSS = checkSym(tmpSS);
    tmpS = getInv(tmpSS);
    S = tmpS * sigma2;
    S = checkSym(S);
    M = tmpS * (tXX * Coefhat);
    
  } else if ((method == "LASSO") || (method == "ALASSO")) {
    tmpSS = tXX + inveta2mat;
    tmpSS = checkSym(tmpSS);
    tmpS = getInv(tmpSS);
    S = tmpS * sigma2;
    S = checkSym(S);
    M = tmpS * (tXX * Coefhat);
  }
  
  if (mono == 0) {
    Coef = rmvnorm(M, S);
  } else {
    for (gg = 0; gg < q; gg++) {
      Sgg = S.row(gg);
      Snotgg = removeRow(S, gg);
      Snotnot = removeCol(Snotgg, gg);
      invSnotgg = getInv(Snotnot);
      
      tmpMat = M(gg) + removeCol(Sgg, gg) * invSnotgg * (removeRow(Coef, gg) - removeRow(M, gg));
      mj = tmpMat(0);
      
      tmpMat = Sgg(gg) - removeCol(Sgg, gg) * invSnotgg * Snotgg.col(gg);
      sj = tmpMat(0);
      
      if (gg == 0) {
        bound1 = abs(bound0);
        bound2 = abs(Coef(1));
      } else if (gg == q - 1) {
        bound1 = abs(Coef(q - 2));
        bound2 = abs(boundqplus1);
      } else {
        bound1 = abs(Coef(gg - 1));
        bound2 = abs(Coef(gg + 1));
      }
      tmp = rtwosegnorm(1, bound2, bound1, mj, sqrt(sj));
      Coef(gg) = tmp(0);
    }
  }
  
  return(Coef);
}


double updateSigma2X(arma::mat resi, arma::mat phi, arma::mat beta, 
                     arma::mat invphieta2mat, arma::mat invbetaeta2mat, 
                     int T, int q, int p, 
                     arma::mat phiA, arma::mat betaA, double a, double b, 
                     Rcpp::String method, int Xflg)  {
  double sigma2 = 0.0;
  
  arma::mat tResiResi = resi.t() * resi;
  double tmptResiResi = tResiResi(0);
  
  arma::mat tCoefVarCoef; 
  double tmptCoefVarCoef;
  
  double tmpa = T / 2.0 + a;
  double tmpb = tmptResiResi / 2.0 + b;
  
  
  arma::mat tmpcoef; 
  
  // update phi contribution
  if (method == "regression") {
    tCoefVarCoef = phi.t() * getInv(phiA) * phi;
    tmptCoefVarCoef = tCoefVarCoef(0);
    tmpa = tmpa + q / 2.0;
    tmpb = tmpb + tmptCoefVarCoef / 2.0;
  } else if ((method == "LASSO") || (method == "ALASSO")) {
    tCoefVarCoef = phi.t() * invphieta2mat * phi;
    tmptCoefVarCoef = tCoefVarCoef(0);
    tmpa = tmpa + q / 2.0;
    tmpb = tmpb + tmptCoefVarCoef / 2.0;
  }
  
  // update beta contribution
  if (Xflg == 1) {
    if (method == "regression") {
      tCoefVarCoef = beta.t() * getInv(betaA) * beta;
      tmptCoefVarCoef = tCoefVarCoef(0);
      tmpa = tmpa + p / 2.0;
      tmpb = tmpb + tmptCoefVarCoef / 2.0;
    } else if ((method == "LASSO") || (method == "ALASSO")) {
      tCoefVarCoef = beta.t() * invbetaeta2mat * beta;
      tmptCoefVarCoef = tCoefVarCoef(0);
      tmpa = tmpa + p / 2.0;
      tmpb = tmpb + tmptCoefVarCoef / 2.0;
    }
  }
  
  sigma2 = 1.0 / R::rgamma(tmpa, 1.0 / (tmpb));
  
  return(sigma2);
}

arma::mat updatelambda2X(arma::mat phieta2, arma::mat betaeta2, int q, int p, double alpha, double beta, 
                         Rcpp::String method, int Xflg){
  arma::mat lambda2(q + p, 1); 
  lambda2.zeros();
  
  double shape;
  double scale;
  
  arma::vec tmpsum; 
  int gg;
  // update lambda2
  if (method == "LASSO") {
    shape = alpha + q + p;
    tmpsum = arma::sum(phieta2);
    if (Xflg == 1) {
      tmpsum = tmpsum + arma::sum(betaeta2);
    }
    scale = beta + tmpsum(0) / 2.0;
    lambda2.fill(R::rgamma(shape, scale));
  } else if (method == "ALASSO") {
    shape = alpha + 1;
    for (gg = 0; gg < q; gg++) {
      scale = beta + phieta2(gg) / 2.0;
      lambda2(gg) = R::rgamma(shape, scale);
    }
    if (Xflg == 1) {
      for (gg = 0; gg < p; gg++) {
        scale = beta + betaeta2(gg) / 2.0;
        lambda2(q + gg) = R::rgamma(shape, scale);
      }
    }
  }
  return(lambda2);
}


arma::mat updateinveta2X(arma::colvec coef, double sigma2, arma::mat lambda2, int q, double tol) {
  arma::mat coefm =arma::conv_to<arma::mat>::from(coef); 
  arma::mat inveta2(q, 1);
  arma::mat muPrime =arma::sqrt(sigma2 * (lambda2 %arma::pow(coefm, -2)));
  arma::colvec tmp; 
  int gg;
  
  double tmpmuPrime;
  double tmplambda2;
  
  for (gg = 0; gg < q; gg++) {
    tmpmuPrime = muPrime(gg, 0);
    tmplambda2 = lambda2(gg, 0);
    //tmp = rrinvgauss(1, tmpmuPrime, tmplambda2);
    if ((-tol < coefm(gg, 0)) && (coefm(gg, 0) < tol)) {
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

arma::mat updateResi(arma::mat B, arma::colvec Phi, int q) {
  int T = B.n_rows;
  int m = B.n_cols;
  arma::mat tmpB = B;
  arma::mat tmpBcol(T, 1);
  arma::mat tmp;
  
  for (int gg = 0; gg < m; gg++) {
    tmpBcol = B.col(gg);
    tmp = getV(tmpBcol, q);
    tmpBcol = tmpBcol - tmp * Phi;
    tmpB.col(gg) = tmpBcol;
  }
  
  return(tmpB);
} 


Rcpp::List updateTauGammaX(arma::colvec Y, arma::mat X_, arma::colvec Phi, arma::colvec Beta, 
                           arma::mat Tau, arma::mat Gamma, 
                           double mu0, double sigma2, double pho, double xi2,
                           int T, arma::mat H_, int m, int Xflg, int q){
  
  //double pi = 3.14159265359;
  
  arma::mat Ht(T, 1); 
  arma::mat tmpHt = Ht; 
  arma::mat tmptHtHt; 
  arma::mat tmptHtHtoversigma2; 
  arma::mat Hnot(T, m - 1);
  
  arma::mat Taunot(m - 1, 1);
  arma::mat Taut(1, 1);
  arma::mat Gammanot(m - 1, 1);
  arma::mat Gammat(1, 1);
  
  arma::mat zetanot(T, 1);
  arma::mat zetat(T, 1);
  
  arma::mat probvec(m, 1); 
  
  arma::mat muGamma = Gamma;
  arma::mat sigma2Gamma = Gamma;
  
  //double tmpzetanot;
  double tmpzetat;
  arma::mat tmp;
  double prob;
  
  arma::mat Gammathat(m, 1);
  double st;
  double mt; 
  
  int jj;
  
  arma::mat tmpY = Y - mu0;
  
  if (Xflg == 1) {
    tmpY = tmpY - X_ * Beta; 
  }
  
  //Rcpp::Rcout << "Y:" << Y << std::endl;
  //Rcpp::Rcout << "tmpY:" << tmpY << std::endl;
  //Rcpp::Rcout << "Tau and Gamma:" << arma::join_rows(Tau, Gamma) << std::endl;
  //Rcpp::Rcout << "Gamma:" << Gamma << std::endl;
  
  // update Tau and Gamma
  for (jj = 0; jj < m; jj++) {
    
    Hnot = removeCol(H_, jj);
    Ht = H_.col(jj);
    
    //Rcpp::Rcout << "Hnot:" << Hnot << std::endl;
    //Rcpp::Rcout << "Ht:" << Ht << std::endl;
    
    ////////////////
    
    Taunot = removeRow(Tau, jj);
    Taut = Tau.row(jj);
    
    Gammanot = removeRow(Gamma, jj);
    Gammat = Gamma.row(jj);
    
    //Rcpp::Rcout << "Taut:" << Taut << std::endl;
    //Rcpp::Rcout << "Gammat:" << Gammat << std::endl;
    
    //Rcpp::Rcout << "Taunot:" << Taunot << std::endl;
    //Rcpp::Rcout << "Gammanot:" << Gammanot << std::endl;
    
    ////////////////
    
    //Rcpp::Rcout << "Taunot:" << Taunot << std::endl;
    
    //update Tau
    zetanot = tmpY - Hnot * (Taunot % Gammanot);
    zetat = zetanot - Ht * Gammat;
    
    //Rcpp::Rcout << "a0:" << Taunot << std::endl;
    //Rcpp::Rcout << "a:" << Taunot % Gammanot << std::endl;
    //Rcpp::Rcout << "b:" << Ht * Gammat << std::endl;
    
    //Rcpp::Rcout << "zetanot:" << zetanot.t() * zetanot << std::endl;
    //Rcpp::Rcout << "zetat:" << zetat.t() * zetat << std::endl;
    
    zetanot = updateResi(zetanot, Phi, q);
    zetat = updateResi(zetat, Phi, q);
    
    //Rcpp::Rcout << "resizetanot:" << zetanot << std::endl;
    //Rcpp::Rcout << "resizetat:" << zetat << std::endl;
    
    //tmp = arma::exp(-1.0 / 2.0 / sigma2 * zetanot.t() * zetanot);
    //tmpzetanot = tmp(0);
    
    //tmp = arma::exp(-1.0 / 2.0 / sigma2 * zetat.t() * zetat);
    //tmpzetat = tmp(0);
    
    //prob = pho * tmpzetat / (pho * tmpzetat + (1 - pho) * tmpzetanot);
    
    tmp = arma::exp(-1.0 / 2.0 / sigma2 * zetanot.t() * zetanot - 
      (-1.0 / 2.0 / sigma2 * zetat.t() * zetat));
    
    tmpzetat = tmp(0);
    
    prob = 1.0 / (1.0 + (1.0 - pho) / pho * tmpzetat);
    
    //Rcpp::Rcout << "tmpzetat:" << tmpzetat << std::endl;
    //Rcpp::Rcout << "prob:" << prob << std::endl;
    
    Tau(jj) = R::rbinom(1, prob);
    probvec(jj) = prob;
    
    //Rcpp::Rcout << prob << std::endl;
    //Rcpp::Rcout << Tau(jj) << std::endl;
    
    //############
    if (Tau(jj) == 1) {
      
      tmpHt = updateResi(Ht, Phi, q);
      
      //#update Gamma
      tmptHtHt = tmpHt.t() * tmpHt;
      tmptHtHtoversigma2 = tmptHtHt / sigma2;
      Gammathat = (1.0 / tmptHtHt) * tmpHt.t() * zetanot;
      
      tmp = 1.0 / (tmptHtHtoversigma2 + 1.0 / xi2);
      st = tmp(0);
      tmp = tmp * tmptHtHtoversigma2 * Gammathat;
      mt = tmp(0);
    } else {
      mt = 0;
      st = xi2;
    }
    
    Gamma(jj) = R::rnorm(mt, sqrt(st));
    muGamma(jj) = mt;
    sigma2Gamma(jj) = st;
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["prob"] = probvec,
    Rcpp::_["muGamma"] = muGamma,
    Rcpp::_["sigma2Gamma"] = sigma2Gamma
  );
  return(out);
  
}

// [[Rcpp::export]]
Rcpp::List initGibbsRFLSMXcpp(arma::colvec Y, Rcpp::List bset, double tol,
                                Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue, 
                                Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue, 
                                Rcpp::Nullable<Rcpp::NumericMatrix> lambda2 = R_NilValue) {
  
  
  
  int phiq = bset["phiq"];
  double gammaxi2 = bset["gammaxi2"];
  Rcpp::String method = bset["method"];
  int updatelambda2 = bset["updatelambda2"];
  
  /////////////////////
  
  int T = Y.n_rows;
  
  Rcpp::List tmpmodel; 
 
  Rcpp::NumericVector tmpcoef;
  arma::mat Phi;
  arma::mat Beta;
  double mu0;
  arma::mat Mu; 
  double sigma2;
   
  arma::mat One(T, 1); 
  
  /////////////////////
  int m;
  arma::mat Tau;
  arma::mat Gamma;
  arma::mat H_;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
    Gamma.randn(m, 1);
    Gamma = Gamma * sqrt(gammaxi2);
    Tau.zeros(m, 1);
  }
  
  /////////////////////
  
  arma::mat X_;
  int Xflg = 0;
  int p = 0;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    tmpmodel = arimaxcpp(Y, phiq, X_);
    Xflg = 1;
    p = X_.n_cols;
  } else {
    tmpmodel = arimacpp(Y, phiq);
  }
  
  arma::mat eta2;
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    eta2.set_size(phiq + p, 1);
  }
  
  arma::mat coef(phiq + p + 1, 1);
  
  
  /////////////////////
  
  sigma2 = tmpmodel["sigma2"];
  tmpcoef = tmpmodel["coef"];
  
  for (int gg = 0; gg < (phiq + p + 1); gg++) {
    coef(gg) = tmpcoef(gg);  
  }
  
  /////////////////////  

  mu0 = coef(phiq);
  Mu = One * mu0;
  
  /////////////////////  
  
  Phi = coef.rows(0, phiq - 1);
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    eta2.rows(0, phiq - 1) = arma::pow(coef.rows(0, phiq - 1), 2);
  }
  
  /////////////////////    
  
  if (Xflg == 1) {
  
    Beta = coef.rows(phiq + 1, phiq + p);
    Mu = Mu + X_ * Beta;
    
    if ((method == "LASSO") || (method == "ALASSO")) {
      eta2.rows(phiq, phiq + p - 1) = arma::pow(coef.rows(phiq + 1, phiq + p), 2);
    }
    
  }
  
  /////////////////////
  
  arma::mat lambda2_;
  double tmpval = 0.0;
  int gg;
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    lambda2_.set_size(phiq + p, 1);
    if (lambda2.isNotNull()) {
    lambda2_ = Rcpp::as<arma::mat>(lambda2);
  } else {
    if (updatelambda2 == 1) {
      if ((method == "LASSO")) {
        tmpval = arma::accu(arma::abs(Phi));
        if (Xflg == 1) {
          tmpval = tmpval + arma::accu(arma::abs(Beta));
        }
        lambda2_.fill(pow((phiq + p) * sqrt(sigma2) / tmpval, 2));
      } else if ((method == "ALASSO")) {
        for (gg = 0; gg < (phiq); gg++) {
          lambda2_(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
        }
        if (Xflg == 1) {
          for (gg = phiq; gg < (phiq + p); gg++) {
            lambda2_(gg) = pow((sqrt(sigma2) / abs(Beta(gg - phiq))), 2);
          }
        }
      }
    }
  }
  }
  
  //lambda2_.fill(tol);
  
  
  
  
  //Rcpp::Rcout << lambda2_ << std::endl;
  
  
  /////////////////////
  
  Rcpp::List out;
  
  out = Rcpp::List::create(
     Rcpp::_["Phi"] = Phi,
     Rcpp::_["Beta"] = Beta,
     Rcpp::_["Tau"] = Tau,
     Rcpp::_["Gamma"] = Gamma,
     Rcpp::_["mu0"] = mu0,
     Rcpp::_["Mu"] = Mu,
     Rcpp::_["eta2"] = eta2,
     Rcpp::_["sigma2"] = sigma2,
     Rcpp::_["lambda2"] = lambda2_
  );
  
  return(out);
  
}



// [[Rcpp::export]]
Rcpp::List simpleinitGibbsRFLSMXcpp(arma::colvec Y, Rcpp::List bset, double tol,
                                Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue, 
                                Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue, 
                                Rcpp::Nullable<Rcpp::NumericMatrix> lambda2 = R_NilValue) {
  
  
  
  int phiq = bset["phiq"];
  //double gammaxi2 = bset["gammaxi2"];
  Rcpp::String method = bset["method"];
  int updatelambda2 = bset["updatelambda2"];
  //int phimono = bset["phimono"];
  
  /////////////////////
  
  int T = Y.n_rows;
  
  Rcpp::List tmpmodel; 
 
  Rcpp::NumericVector tmpcoef;
  arma::mat Phi;
  arma::mat Beta;
  double mu0;
  arma::mat Mu; 
  double sigma2;
   
  arma::mat One(T, 1); 
  One.ones();
  
  /////////////////////
  int m;
  arma::mat Tau;
  arma::mat Gamma;
  arma::mat H_;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
    //Gamma.randn(m, 1);
    //Gamma = Gamma * sqrt(gammaxi2);
    Gamma.set_size(m, 1);
    Gamma.fill(tol);
    Tau.zeros(m, 1);
  }
  
  /////////////////////
  
  //Rcpp::Rcout << 3.1 << std::endl;
  
  arma::mat X_;
  int Xflg = 0;
  int p = 0;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    Xflg = 1;
    p = X_.n_cols;
    Beta.set_size(p, 1);
    Beta.fill(tol);
  } 
  
  int nn = 0;
  
  //if (phimono == 2) {
  //  nn = phiq * 2 + p;
  //} else {
    nn = phiq + p;
  //}
  
  /////////////////////
  
  arma::mat eta2;
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    eta2.set_size(nn, 1);
  }
  
  /////////////////////  

  //mu0 = coef(phiq);
  mu0 = arma::accu(Y) / T;
  Mu = One * mu0;
  
  //Rcpp::Rcout << 3.22 << std::endl;
  
  /////////////////////  
  
  //if (phimono == 2) {
    
  //  Phi.set_size(phiq * 2, 1);
    
  //} else {
    
    Phi.set_size(phiq, 1);
  
  //}
  Phi.fill(tol);
  
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    //if (phimono == 2) {
    //  eta2.rows(0, phiq * 2 - 1) = arma::pow(Phi, 2);
    //} else {
      eta2.rows(0, phiq - 1) = arma::pow(Phi, 2);
    //}
  }
  
  
  //Rcpp::Rcout << 3.23 << std::endl;
  
  /////////////////////    
  
  if (Xflg == 1) {
  
    Mu = Mu + X_ * Beta;
    
    if ((method == "LASSO") || (method == "ALASSO")) {
      //if (phimono == 2) {
      //  eta2.rows(phiq * 2, phiq * 2 + p - 1) = arma::pow(Beta, 2);
      //} else {
        eta2.rows(phiq, phiq + p - 1) = arma::pow(Beta, 2);
      //}
    }
    
  }
  
  sigma2 = arma::accu(arma::pow(Y - Mu, 2)) / T;
  
  //Rcpp::Rcout << 3.3 << std::endl;
  
  /////////////////////
  
  arma::mat lambda2_;
  double tmpval = 0.0;
  int gg;
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    lambda2_.set_size(phiq + p, 1);
    if (lambda2.isNotNull()) {
    lambda2_ = Rcpp::as<arma::mat>(lambda2);
  } else {
    if (updatelambda2 == 1) {
      if ((method == "LASSO")) {
        
          tmpval = arma::accu(arma::abs(Phi));
          if (Xflg == 1) {
            tmpval = tmpval + arma::accu(arma::abs(Beta));
          }
          lambda2_.fill(pow((nn) * sqrt(sigma2) / tmpval, 2));
        
      } else if ((method == "ALASSO")) {
        
        //if (phimono == 2) {
        //  
        //  for (gg = 0; gg < (phiq * 2); gg++) {
        //    lambda2_(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
        //  }
        //  if (Xflg == 1) {
        //    for (gg = (phiq * 2); gg < (nn); gg++) {
        //      lambda2_(gg) = pow((sqrt(sigma2) / abs(Beta(gg - (phiq * 2)))), 2);
        //    }
        //  }
        //  
        //} else {
          for (gg = 0; gg < (phiq); gg++) {
            lambda2_(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
          }
          if (Xflg == 1) {
            for (gg = phiq; gg < (phiq + p); gg++) {
              lambda2_(gg) = pow((sqrt(sigma2) / abs(Beta(gg - phiq))), 2);
            }
          }
        //}
        
      }
    }
  }
  }
  
  //lambda2_.fill(tol);
  
  
  //Rcpp::Rcout << 3.4 << std::endl;
  
  //Rcpp::Rcout << lambda2_ << std::endl;
  
  
  /////////////////////
  
  Rcpp::List out;
  
  out = Rcpp::List::create(
     Rcpp::_["Phi"] = Phi,
     Rcpp::_["Beta"] = Beta,
     Rcpp::_["Tau"] = Tau,
     Rcpp::_["Gamma"] = Gamma,
     Rcpp::_["mu0"] = mu0,
     Rcpp::_["Mu"] = Mu,
     Rcpp::_["eta2"] = eta2,
     Rcpp::_["sigma2"] = sigma2,
     Rcpp::_["lambda2"] = lambda2_
  );
  
  return(out);
  
}


// [[Rcpp::export]]
Rcpp::List GibbsRFLSMXUpdatecpp(arma::colvec Y, Rcpp::List pars, Rcpp::List bset,
                                double tol, 
                                Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  int T = Y.n_rows;
  
  arma::mat One(T, 1);
  One.ones();
  
  arma::mat X_;
  arma::mat Beta;
  int Xflg = 0;
  
  // Calculate X
  int p = 0;
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    p = X_.n_cols;
    Beta = Rcpp::as<arma::mat>(pars["Beta"]);
    Xflg = 1;
  } 
  
  arma::mat H_;
  arma::mat Tau;
  arma::mat Gamma; 
  int Hflg = 0;
  
  // Calculate H
  int m = 0;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
    Tau = Rcpp::as<arma::mat>(pars["Tau"]);
    Gamma = Rcpp::as<arma::mat>(pars["Gamma"]);
    Hflg = 1;
  } 
  
  arma::mat Phi = Rcpp::as<arma::mat>(pars["Phi"]); 
  double sigma2 = pars["sigma2"];
  arma::mat lambda2 = Rcpp::as<arma::mat>(pars["lambda2"]);
  
  arma::mat Mu = Rcpp::as<arma::mat>(pars["Mu"]);
  arma::mat eta2 = Rcpp::as<arma::mat>(pars["eta2"]);
  arma::mat inveta2 = eta2; 
  
  double mu0 = pars["mu0"];
  
  
  // read bset
  Rcpp::String method = bset["method"];
  int phimono = bset["phimono"];
  int phiq = bset["phiq"];
  arma::mat phiA = Rcpp::as<arma::mat>(bset["phiA"]);
  double phibound0 = bset["phibound0"];
  double phiboundqplus1 = bset["phiboundqplus1"];
  arma::mat betaA = Rcpp::as<arma::mat>(bset["betaA"]);
  double gammaxi2 = bset["gammaxi2"];
  double tautheta1 = bset["tautheta1"];
  double tautheta2 = bset["tautheta2"];
  double sigma2a = bset["sigma2a"];
  double sigma2b = bset["sigma2b"];
  double lambda2alpha = bset["lambda2alpha"];
  double lambda2beta = bset["lambda2beta"];
  int updatelambda2 = bset["updatelambda2"];
  
  // initialize the computation
  arma::mat tmpResiY;
  arma::mat tmpResiX; 
  
  arma::mat phieta2;
  arma::mat invphieta2;
  arma::mat invphieta2mat;
  
  arma::mat betaeta2; 
  arma::mat invbetaeta2;
  arma::mat invbetaeta2mat;
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    phieta2 = eta2.rows(0, phiq - 1);
    invphieta2 = arma::pow(phieta2.rows(0, phiq - 1), -1);
    invphieta2mat.set_size(phiq, phiq); 
    invphieta2mat.zeros();
    invphieta2mat.diag() = invphieta2;
    
    invbetaeta2mat.set_size(p, p);
    invbetaeta2mat.zeros();
    if (Xflg == 1) {
      betaeta2 = eta2.rows(phiq, phiq + p - 1);
      invbetaeta2 = arma::pow(betaeta2, -1);
      invbetaeta2mat.diag() = invbetaeta2;
    }
  }
  
  Rcpp::List TauGamma;
  arma::mat tmpSumTau;
  double pho;
  arma::colvec coef(phiq + p, 1); 
  
  double mu0hat;
  
  /////////////////////////////////////   
  // update Phi
  tmpResiY = Y - Mu;
  tmpResiX = getV(tmpResiY, phiq);
  
  Phi = updateCoef(tmpResiY, tmpResiX, 
                   phiA, Phi, 
                   sigma2, invphieta2mat, 
                   phibound0, phiboundqplus1,
                   phimono, method);
  
  coef.rows(0, phiq - 1) = Phi;
  
  // update beta
  if (Xflg == 1) {
    tmpResiY = Y - mu0;
    if (Hflg == 1) {
      tmpResiY = tmpResiY - H_ * (Tau % Gamma);
    }
    tmpResiY = updateResi(tmpResiY, Phi, phiq);
    tmpResiX = updateResi(X_, Phi, phiq);
    
    Beta = updateCoef(tmpResiY, tmpResiX, 
                      betaA, Beta, 
                      sigma2, invbetaeta2mat, 
                      phibound0, phiboundqplus1,
                      0, method);
    
    coef.rows(phiq, phiq + p - 1) = Beta;
  }
  
  
  // update Tau and Gamma
  if (Hflg == 1) {
    //#update pho
    tmpSumTau = arma::sum(Tau);
    pho = R::rbeta(tautheta1 + tmpSumTau(0), tautheta2 + m - tmpSumTau(0));
    
    
    //Rcpp::Rcout << Y << std::endl;
    //Rcpp::Rcout << X_ << std::endl;
    //Rcpp::Rcout << Phi << std::endl;
    //Rcpp::Rcout << Beta << std::endl;
    //Rcpp::Rcout << Tau << std::endl;
    //Rcpp::Rcout << Gamma << std::endl;
    //Rcpp::Rcout << mu0 << std::endl;
    //Rcpp::Rcout << sigma2 << std::endl;
    //Rcpp::Rcout << pho << std::endl;
    //Rcpp::Rcout << gammaxi2 << std::endl;
    //Rcpp::Rcout << m << std::endl;
    //Rcpp::Rcout << phiq << std::endl;
    
    TauGamma = updateTauGammaX(Y, X_, Phi, Beta, 
                               Tau, Gamma, 
                               mu0, sigma2, pho, gammaxi2,
                               T, H_, m, Xflg, phiq);
    
    Tau = Rcpp::as<arma::mat>(TauGamma["Tau"]);
    Gamma = Rcpp::as<arma::mat>(TauGamma["Gamma"]);
    
    //Rcpp::Rcout << Tau << std::endl;
    //Rcpp::Rcout << Gamma << std::endl;
    
  }
  
  
  // update mu0
  tmpResiY = Y;
  if (Xflg == 1) {
    tmpResiY = tmpResiY - X_ * Beta;
  }
  if (Hflg == 1) {
    tmpResiY = tmpResiY - H_ * (Tau % Gamma);
  }
  tmpResiY = updateResi(tmpResiY, Phi, phiq);
  mu0hat = arma::accu(tmpResiY) / (T + 0.0);
  mu0 = R::rnorm(mu0hat, sqrt(sigma2 / (T + 0.0)));
  
  
  // update Mu
  
  Mu = mu0 * One;
  
  if (Hflg == 1) {
    Mu = Mu + H_ * (Tau % Gamma);
  }
  
  if (Xflg == 1) {
    Mu = Mu + X_ * Beta;
  }
  
  
  // update sigma2
  // update eta2
  
  if ((method == "LASSO") || (method == "ALASSO")) {
    inveta2 = updateinveta2X(coef, sigma2, lambda2, phiq + p, tol);
    eta2 = arma::pow(inveta2, -1);
    
    invphieta2 = inveta2.rows(0, phiq - 1);
    
    phieta2 = eta2.rows(0, phiq - 1);
    
    invphieta2mat.diag() = invphieta2;
    
    
    if(Xflg == 1) {
      invbetaeta2 = inveta2.rows(phiq, phiq + p - 1);
      betaeta2 = eta2.rows(phiq, phiq + p - 1);
      invbetaeta2mat.diag() = invbetaeta2;
    }
  }
  
  
  
  tmpResiY = Y - Mu;
  tmpResiY = updateResi(tmpResiY, Phi, phiq);
  sigma2 = updateSigma2X(tmpResiY, Phi, Beta, 
                         invphieta2mat, invbetaeta2mat, 
                         T, phiq, p, 
                         phiA, betaA, sigma2a, sigma2b, 
                         method, Xflg);
  
  
  // update lambda2
  
  if (updatelambda2 == 1) {
    if ((method == "LASSO") || (method == "ALASSO")) {
      lambda2 = updatelambda2X(phieta2, betaeta2, phiq, p, lambda2alpha, lambda2beta, method, Xflg);
    }
    
  }
  
  /////////////////////////////////////   
  // output the result
  Rcpp::List out;
  
  out = Rcpp::List::create(
    Rcpp::_["Phi"] = Phi,
    Rcpp::_["Beta"] = Beta,
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["mu0"] = mu0, 
    Rcpp::_["Mu"] = Mu, 
    Rcpp::_["eta2"] = eta2,
    Rcpp::_["sigma2"] = sigma2,
    Rcpp::_["lambda2"] = lambda2
  );
  
  return(out);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::mat lhfX(arma::colvec Y, arma::mat Phi, arma::mat Mu, double sigma2) {
 
   double pi = 3.14159265359;
   
   int q = Phi.n_rows;
   
   arma::mat V; 
   arma::mat Vas;
   arma::mat VasPhi;
   arma::mat resi; 
   
   V = Y - Mu;
   Vas = getV(V, q);
   
   VasPhi = Vas * Phi;
   
   resi = V - VasPhi;
   
   arma::mat tmp = sqrt(1 / 2.0 / pi / sigma2) * exp(- 1.0 / 2.0 * arma::pow(resi, 2) / sigma2);
   
   return(tmp);
   
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::mat llhfX(arma::colvec Y, arma::mat Phi, arma::mat Mu, double sigma2) {
 
   double pi = 3.14159265359;
   
   int q = Phi.n_rows;
   //int T = Y.n_elem;
   
   arma::mat V; 
   arma::mat Vas;
   arma::mat VasPhi;
   arma::mat resi; 
   
   V = Y - Mu;
   Vas = getV(V, q);
   
   VasPhi = Vas * Phi;
   
   resi = V - VasPhi;
   
   arma::mat tmp = sqrt(1 / 2.0 / pi / sigma2) * exp(- 1.0 / 2.0 * arma::pow(resi, 2) / sigma2);
   tmp = arma::log(tmp);
   
   return(tmp);
   
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::mat lhYJfX(arma::colvec Y, arma::mat Phi, arma::mat Mu, double sigma2, double theta, double eps) {
 
 int T = Y.n_elem;
 
 arma::colvec Yyj = yeojohnsontr(Y, theta, eps);
 arma::mat lhYJ = lhfX(Yyj, Phi, Mu, sigma2);
 double tmp = 1.0;
 
 for (int i = 0; i < T; i++) {
   tmp = tmp * pow(abs(Y(i)) + 1, (theta - 1) * sign(Y(i)));
 }
 
 //Rcpp::Rcout << lhYJ << std::endl;
 //Rcpp::Rcout << tmp << std::endl;
 
 
 
 return(Yyj);
 
}
 
//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
double llhYJfX(arma::colvec Y, arma::mat Phi, arma::mat Mu, double sigma2, double theta, double eps) {
 
 int T = Y.n_elem;
 
 arma::colvec Yyj = yeojohnsontr(Y, theta, eps);
 arma::mat llhYJ = llhfX(Yyj, Phi, Mu, sigma2);
 double tmp = arma::accu(llhYJ);
 
 for (int i = 0; i < T; i++) {
   tmp = tmp + log(pow(abs(Y(i)) + 1, (theta - 1) * sign(Y(i))));
 }
 
 
 return(tmp);
 
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
double llhYJfXt(arma::colvec Y, int t, arma::mat Phi, arma::mat Mu, 
                double sigma2, double theta, double eps) {
 
 int q = Phi.n_rows;
 int T = Y.n_elem;
 
 arma::colvec Yyj = yeojohnsontr(Y, theta, eps);
 arma::mat llhYJ = llhfX(Yyj, Phi, Mu, sigma2);
 
 //int m = t + q;
 //if (m > (T - 1)) {
 //  m = T - 1;
 //}
 
 //arma::mat llhYJt = llhYJ.rows(t, m); 
 //double tmp = arma::accu(llhYJt);
 
 double tmp = 0;
 
 //for (int i = t; i <= m; i++) {
 for (int i = 0; i < T; i++) {
   tmp = llhYJ(t) + log(pow(abs(Y(i)) + 1, (theta - 1) * sign(Y(i))));
 }
 
 
 return(tmp);
 
}
 

// [[Rcpp::export]]
arma::mat thetaYeoJohnsonMHX(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                        double oldtheta, int burnin, int nsim, double tol) {
  
  //double pi = 3.14159265359;
  
  //int T = Y.n_elem;
  //int q = Phi.n_rows;
  
  double u;
  double oldtheta_ = oldtheta;
  double thetaas;
  double A;
  
 arma::mat thetaout(nsim, 1); 
  
 double oldllhYJ = llhYJfX(Y, Phi, Mu, sigma2, oldtheta_, tol);

 double newllhYJ; 
 arma::uvec ind1;
  
 double tmp; 
  double pd;
  
  int i;
  int j = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    u = R::runif(0.0, 1.0);
    thetaas = R::rnorm(oldtheta_, 0.1);
    newllhYJ = llhYJfX(Y, Phi, Mu, sigma2, thetaas, tol);

    //tmp = newllhYJ - oldllhYJ;
    
    tmp = newllhYJ - oldllhYJ + (R::dnorm4(thetaas, 1, 0.1, 1) - R::dnorm4(oldtheta_, 1, 0.1, 1));
    
    //Rcpp::Rcout << newllhYJ << std::endl;
    //Rcpp::Rcout << oldllhYJ << std::endl;
    
    tmp = exp(tmp);
    
    //Rcpp::Rcout << tmp << std::endl;
    
    //Rcpp::Rcout << tmp << std::endl;
    pd = tmp;
    //Rcpp::Rcout << pd << std::endl;
    A = std::min(1.0, pd);
    //Rcpp::Rcout << tmp(T - q - 1) << std::endl;
    //Rcpp::Rcout << A << std::endl;
    
    if (u < A) {
      oldtheta_ = thetaas;
      oldllhYJ = newllhYJ;
    } 
    
    //Rcpp::Rcout << oldtheta_ << std::endl;
    
    if (i >= burnin) {
      thetaout(j, 0) = oldtheta_;
      j = j + 1;
    }
  }
  
  return(thetaout);
  
}

double dtrnorm(double x, double mean, double sd, double lower, double upper) {
  double beta = (upper - mean) / sd;
  double alpha = (lower - mean) / sd;
  double z = R::pnorm5(beta, 0.0, 1.0, 1, 0) - R::pnorm5(alpha, 0.0, 1.0, 1, 0);
  double xi = (x - mean) / sd;
  double out = 0.0;
  if ((alpha <= xi) && (xi <= beta)) {
    out = R::dnorm4(xi, 0, 1, 0) / sd / z;
  } 
  
  return(out);
}


double llYZt(arma::colvec YZ, arma::mat Phi,arma::mat Mu, 
             double sigma2, double theta, double tol, int t) {
  int q = Phi.n_rows;
  arma::colvec tmpYZ(q + 1);
  tmpYZ.zeros();
  arma::colvec tmpYZyj(q + 1);
  tmpYZyj.zeros();
  arma::colvec resi(q + 1);
  resi.zeros();

  int m;
  
  arma::colvec tmp; 
  
  for (int i = 0; i <= q; i++) {
    m = t - i;
    if (m >= 0) {
      tmpYZ(i) = YZ(m);
      tmp = yeojohnsontr(tmpYZ.row(i), theta, tol);
      tmpYZyj(i) = tmp(0);
      resi(i) = tmpYZyj(i) - Mu(m);
    }
  }
  
  tmp = resi(0) - resi.rows(1, q) * Phi;
  tmp = tmp / sqrt(sigma2);
  
  double tmpllt = R::dnorm4(tmp(0), 0, 1, 1) + log(pow(abs(YZ(t)) + 1, (theta - 1) * sign(YZ(t))));
  
  return(tmpllt);
  
}

double llf(arma::colvec resi, arma::colvec YZ, double sigma2, double theta){
  
  int T = YZ.n_rows;
  arma::colvec tmp(T); 
  
  for (int i = 0; i < T; i++) {
    tmp(i) = R::dnorm4(resi(i), 0.0, sqrt(sigma2), 1) + 
      (theta - 1.0) * sign(YZ(i)) * log(abs(YZ(i)) + 1.0);
  }
  
  double out = arma::accu(tmp);
  return(out);
  
} 

double updateYZt(arma::colvec YZ, 
                arma::mat Phi,arma::mat Mu, double sigma2, double theta,
                int t, double lb, double ub, double tol, int burnin) {
  
  int T = YZ.n_rows;
  int q = Phi.n_rows;
  
  arma::colvec oldYZ = YZ; 
  arma::colvec oldYZyj = yeojohnsontr(oldYZ, theta, tol);
  arma::colvec oldresi = oldYZyj - Mu; 
  oldresi = updateResi(oldresi, Phi, q); 
  
  //Rcpp::Rcout << 2.01  << std::endl;
  //Rcpp::Rcout << "oldresi:" << oldresi  << std::endl;
  //Rcpp::Rcout << "oldYZ:" << oldYZ  << std::endl;
  
  double oldll = llf(oldresi, oldYZ, sigma2, theta);
  
  //Rcpp::Rcout << "oldll:" << oldll  << std::endl;
  
  arma::colvec newYZ = oldYZ; 
  arma::colvec newYZyj = oldYZyj;
  
  arma::colvec YZtas;
  arma::colvec YZyjtas;
  arma::colvec resias;
  
  double newll;
  
  double pd;
  double A;
  double u;
  
  //Rcpp::Rcout << 2.1  << std::endl;
  
  for (int i = 0; i < (burnin + 1); i++) {
    
    u = R::runif(0, 1);
    
    YZtas = rtrnorm(1, oldYZ(t), 0.1, lb, ub);
    //Rcpp::Rcout << "YZtas:" << YZtas  << std::endl;
    
    YZyjtas = yeojohnsontr(YZtas, theta, tol);
    //Rcpp::Rcout << "YZyjtas:" << YZyjtas  << std::endl;
    
    newYZyj.row(t) = YZyjtas;
    newYZ.row(t) = YZtas;
    
    resias = newYZyj - Mu;
    resias = updateResi(resias, Phi, q); 
    
    newll = llf(resias, newYZ, sigma2, theta);
   
    //Rcpp::Rcout << 2.2  << std::endl;
   
    pd = newll - log(dtrnorm(newYZ(t), oldYZ(t), 0.1, lb, ub)) - 
      (oldll - log(dtrnorm(oldYZ(t), newYZ(t), 0.1, lb, ub)));
    
    pd = exp(pd);
    
    A = std::min(1.0, pd);
    if (u < A) {
      oldll = newll;
      oldYZyj = newYZyj;
      oldYZ = newYZ;
    } 
    
    //Rcpp::Rcout << "lb:" << lb  << std::endl;
    //Rcpp::Rcout << "ub:" << ub  << std::endl;
    //Rcpp::Rcout << "newYZyj.row(t):" << newYZyj.row(t)  << std::endl;
    //Rcpp::Rcout << "newYZ.row(t):" << newYZ.row(t)  << std::endl;
    //Rcpp::Rcout << "oldYZyj(t):" << oldYZyj(t)  << std::endl;
    //Rcpp::Rcout << "oldYZ(t):" << oldYZ(t)  << std::endl;
    //Rcpp::Rcout << "pd:" << pd  << std::endl;
    
    //Rcpp::Rcout << "oldll:" << oldll  << std::endl;
    
  }
  
  //Rcpp::Rcout << "oldYZyj(t):" << oldYZyj(t)  << std::endl;
  //Rcpp::Rcout << "oldYZ(t):" << oldYZ(t)  << std::endl;
  //Rcpp::Rcout << "theta:" << theta  << std::endl;
  
  double out = oldYZ(t);
  return(out);
  
}



// [[Rcpp::export]]
arma::mat getYZMHX(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta,
                        arma::colvec oldZ, int leftcensoring, int rounding, 
                        int burnin, int nsim, double tol) {
  
  int T = Y.n_elem;
  
  arma::mat oldYZ = Y + oldZ;
  arma::mat newYZ = oldYZ;
  
  double tmpYZt;
  arma::mat YZout(T, nsim);
  
  double lbr;
  double ubr;
  
  double lbl;
  double ubl;
  
  double ub;
  double lb;
  
  arma::colvec tmp; 
  
  int flgr;
  int flgl;
  
  int i;
  int t = 0;

  //Rcpp::Rcout << 1  << std::endl;
  
    for (i = 0; i < nsim; i++) {
      for (t = 0; t < T; t++) {
        
         if ((rounding == 1) || (leftcensoring == 1)) {
        
            flgr = 0;
            flgl = 0 ;
            
            if (rounding == 1) {
              tmp = -0.5 + Y.row(t);
              //tmp = yeojohnsontr(tmp, theta, tol);
              lbr = tmp(0);
              
              tmp = 0.5 + Y.row(t);
              //tmp = yeojohnsontr(tmp, theta, tol);
              ubr = tmp(0);
              flgr = 1;
            }
            
            if (leftcensoring == 1) {
              if (Y(t) <= 0.0) {
                lbl = (-1) * arma::datum::inf;
                ubl = 0.0;
                flgl = 1;
              }
            }
            
            
            
            if ((flgr == 0) && (flgl == 0)) {
              newYZ(t) = oldYZ(t);
              //newYZ(t) = Y(t) + oldZ(t);
            } else {
              if ((flgr == 1) && (flgl == 1)) {
                lb = lbl;
                ub = ubr;
              } else if ((flgr == 1) && (flgl == 0)) {
                lb = lbr;
                ub = ubr;
              } else if ((flgr == 0) && (flgl == 1)) {
                lb = lbl;
                ub = ubl;
              }
              
              //Rcpp::Rcout << "lb:" << lb << std::endl;
              //Rcpp::Rcout << "ub:"  << ub  << std::endl;
              //Rcpp::Rcout << "flgr:"<< flgr   << std::endl;
              //Rcpp::Rcout << "flgl:"<< flgl   << std::endl;
              
              //Rcpp::Rcout << 2  << std::endl;
              
              tmpYZt = updateYZt(newYZ, Phi, Mu, sigma2, theta,
                t, lb, ub, tol, burnin);
              
              //Rcpp::Rcout << "t:" << t  << std::endl;
              //Rcpp::Rcout << "tmpYZt:" << tmpYZt  << std::endl;
              //Rcpp::Rcout << 3  << std::endl;
              
              newYZ(t) = tmpYZt;
              
            }
      
      
      }
    }
    
    //YZout.col(i) = Y + newZ;
    YZout.col(i) = newYZ;
  }
      
  return(YZout);
  
}


double updateZt(arma::colvec Z, arma::colvec Y, 
                arma::mat Phi,arma::mat Mu, double sigma2, double theta,
                int t, double lb, double ub, double tol, int burnin) {
  
  int T = Y.n_rows;
  int q = Phi.n_rows;
  
  arma::colvec oldZ = Z; 
  
  arma::colvec oldYZ = Y + Z; 
  arma::colvec oldYZyj = yeojohnsontr(oldYZ, theta, tol);
  arma::colvec oldresi = oldYZyj - Mu; 
  oldresi = updateResi(oldresi, Phi, q); 
  
  double oldll = llf(oldresi, oldYZ, sigma2, theta);
  
  arma::colvec newYZ = oldYZ; 
  arma::colvec newYZyj = oldYZyj;
  
  arma::colvec Ztas;
  arma::colvec Zyjtas;

  arma::colvec YZtas;
  arma::colvec YZyjtas;
  arma::colvec resias;
  
  arma::colvec tmp; 
  
  double newll;
  
  double pd;
  double A;
  double u;
  
  //Rcpp::Rcout << 2.1  << std::endl;
  
  for (int i = 0; i < (burnin + 1); i++) {
    
    u = R::runif(0, 1);
    
    Ztas = rtrnorm(1, oldZ(t), 0.1, lb, ub);
    YZtas = Y(t) + Ztas;
    
    YZyjtas = yeojohnsontr(YZtas, theta, tol);
    //Rcpp::Rcout << "YZyjtas:" << YZyjtas  << std::endl;
    
    newYZyj.row(t) = YZyjtas;
    newYZ.row(t) = YZtas;
    
    resias = newYZyj - Mu;
    resias = updateResi(resias, Phi, q); 
    
    newll = llf(resias, newYZ, sigma2, theta);
   
    //Rcpp::Rcout << 2.2  << std::endl;
   
    pd = newll + log(dtrnorm(Ztas(0), 0, 0.1, lb, ub)) - log(dtrnorm(Ztas(0), oldZ(t), 0.1, lb, ub)) - 
      (oldll + log(dtrnorm(oldZ(t), 0, 0.1, lb, ub)) - log(dtrnorm(oldZ(t), Ztas(0), 0.1, lb, ub)));
    
    pd = exp(pd);
    
    A = std::min(1.0, pd);
    if (u < A) {
      oldZ.row(t) = Ztas;
      oldll = newll;
      oldYZyj = newYZyj;
      oldYZ = newYZ;
    } 
    
    
  }
  
  double out = oldZ(t);
  return(out);
  
}



// [[Rcpp::export]]
arma::mat getZZ(arma::colvec Y, arma::mat Phi,arma::mat Mu, double sigma2, double theta,
                        arma::colvec oldZ, int leftcensoring, int rounding, 
                        int burnin, int nsim, double tol) {
  
  int T = Y.n_elem;
  
  arma::mat oldYZ = Y + oldZ;
  arma::mat newYZ = oldYZ;
  
  arma::mat newZ = oldZ;
  
  //double tmpYZt;
  double tmpZt;
  arma::mat Zout(T, nsim);
  
  double lbr;
  double ubr;
  
  double lbl;
  double ubl;
  
  double ub;
  double lb;
  
  arma::colvec tmp; 
  
  int flgr;
  int flgl;
  
  int i;
  int t = 0;

  //Rcpp::Rcout << 1  << std::endl;
  
    for (i = 0; i < nsim; i++) {
      for (t = 0; t < T; t++) {
        
         if ((rounding == 1) || (leftcensoring == 1)) {
        
            flgr = 0;
            flgl = 0 ;
            
            if (rounding == 1) {
              lbr = -0.5;
              
              ubr = 0.5;
              flgr = 1;
            }
            
            if (leftcensoring == 1) {
              if (Y(t) <= 0.0) {
                lbl = (-1.0) * arma::datum::inf;
                ubl = 0.0;
                flgl = 1;
              }
            }
            
            
            
            if ((flgr == 0) && (flgl == 0)) {
              newZ(t) = oldZ(t);
              //newYZ(t) = Y(t) + oldZ(t);
            } else {
              if ((flgr == 1) && (flgl == 1)) {
                lb = lbl;
                ub = ubr;
              } else if ((flgr == 1) && (flgl == 0)) {
                lb = lbr;
                ub = ubr;
              } else if ((flgr == 0) && (flgl == 1)) {
                lb = lbl;
                ub = ubl;
              }
              
              newZ(t) = updateZt(newZ, Y, Phi, Mu, sigma2, theta, t, lb, ub, tol, burnin);
              
              
            }
      
      
      }
    }
    
    //YZout.col(i) = Y + newZ;
    Zout.col(i) = newZ;
  }
      
  return(Zout);
  
}



// [[Rcpp::export]]
arma::colvec getYZ(arma::colvec Yyj,arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                   double theta, double eps, int leftcensoring, int rounding, int updateYJ) {
  
  double lowerbound = 0.0;
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Y.n_elem;
  
  arma::colvec YZ = Y;
  arma::colvec YZyj = Yyj; 
  
  arma::mat V = YZyj - Mu; 
  arma::mat Vas(1, q);
  arma::mat VasPhi;
  
  arma::mat fit(T, 1);
  arma::colvec tmp(1);
  
  double lbr;
  double ubr;
  
  double lbl;
  double ubl;
  
  double lb;
  double ub;
  
  arma::colvec lbvec(T, 1);
  arma::colvec ubvec(T, 1);
  
  int flgr = 0;
  int flgl = 0;
  
  //Rcpp::Rcout << "leftcensoring:" << leftcensoring << std::endl;
  //Rcpp::Rcout << "rounding:" << rounding  << std::endl;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  //Rcpp::Rcout << "Phi:" << Phi  << std::endl;
  
  int ii;
  if ((leftcensoring == 1) || (rounding == 1)) {
    for (int i = 0; i < T; i++) {
      
      /////////////////////////////////
      for (int j = 0; j < q; j++) {
        ii = i - 1 - j;
        if (ii >= 0) {
          Vas(0, j) = V(ii, 0);
        } else {
          Vas(0, j) = 0;
        }
      }
      
      VasPhi = Vas * Phi;
      fit(i, 0) = Mu(i, 0) + VasPhi(0);
      
      //Rcpp::Rcout << "Mu:" << Mu(i, 0) << std::endl;
      //Rcpp::Rcout << "VasPhi:" << VasPhi(0) << std::endl;
      
      /////////////////////////////////
      
      flgr = 0;
      
      if (rounding == 1) {
        if (updateYJ == 1) {
          tmp = yeojohnsontr(YZ.row(i) - 0.5, theta, eps);
          lbr = tmp(0);
          tmp = yeojohnsontr(YZ.row(i) + 0.5, theta, eps);
          ubr = tmp(0);
        } else {
          lbr = YZ(i) - 0.5;
          ubr = YZ(i) + 0.5;
        }
        flgr = 1;
      }
      
      flgl = 0;
      
      if (leftcensoring == 1) {
        lbl = (-1.0) *arma::datum::inf;
        if (YZ(i) <= lowerbound) {
          if (updateYJ == 1) {
            tmp(0) = lowerbound;
            tmp = yeojohnsontr(tmp, theta, eps);
            ubl = tmp(0);
          } else{
            ubl = lowerbound;
          }
          
          flgl = 1;
        }
      }
      
      //////// find the union
      
      if (flgr == 1) {
        lb = lbr;
        ub = ubr;
      }
      
      if (flgl == 1) {
        lb = lbl;
        ub = ubl;
        if (flgr == 1) {
          if (ubr > ubl) {
            ub = ubr;
          }
        } 
      }
        
        //lbvec(i) = lb;
        //ubvec(i) = ub;
        
        
        
        if ((flgr == 1) || (flgl == 1)) {
          tmp = rtrnorm(1, fit(i, 0), sqrt(sigma2), lb, ub);
          YZyj(i) = tmp(0);
          tmp = invyeojohnsontr(YZyj.row(i), theta, eps);
          tmp.replace(arma::datum::nan, (-1.0) * arma::datum::inf);
          YZ.row(i) = tmp;
        }
        
        //Rcpp::Rcout << "i:" << i  << std::endl;
        //Rcpp::Rcout << "flgr + flgl:" << flgr + flgl  << std::endl;
        //Rcpp::Rcout << "fit(i, 0):" << fit(i, 0)  << std::endl;
        //Rcpp::Rcout << "sigma2:" << sigma2  << std::endl;
        //Rcpp::Rcout << "lb:" << lb  << std::endl;
        //Rcpp::Rcout << "ub:" << ub  << std::endl;
        //Rcpp::Rcout << "tmp:" << tmp  << std::endl;
        //Rcpp::Rcout << "YZ:" << YZ(i)  << std::endl;
        //Rcpp::Rcout << "YZyj:" << YZyj(i)  << std::endl;
      
      //////////////////////////// 
      
      V(i, 0) = YZyj(i) - Mu(i, 0);
      
    }
  }
  
  return(YZ);
  
}




// [[Rcpp::export]]
Rcpp::List GibbsRFLSMXcpp(arma::colvec Y, 
                          Rcpp::List bset, double tol, 
                          int nsim, int thin, int burnin, 
                          int verbose,
                          Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericMatrix> lambda2 = R_NilValue,
                          Rcpp::Nullable<double> theta = R_NilValue) {
  
  //int YJ = bset["YJ"];
  int leftcensoring = bset["leftcensoring"];
  //int lowerbound = bset["lowerbound"];
  //double lowerbound = 0.0;
  int rounding = bset["rounding"];
  int updateYJ = bset["updateYJ"];
  int phiq = bset["phiq"];
  Rcpp::String method = bset["method"];
  int phimono = bset["phimono"];
  int updatelambda2 = bset["updatelambda2"];
  
  //Rcpp::Rcout << 1 << std::endl;
  
  /////////////////////////////////////
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  Rcpp::String monoword;
  
  if (phimono == 1) {
    monoword = " with constraints on Phi";
  } else {
    monoword = " without constraints on Phi";
  }
  
  if (verbose == 1) {
    Rcpp::Rcout << "Start training using " << method.get_cstring() << monoword.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  }
  
  /////////////////////////////////////
  
  arma::mat Phi;
  arma::mat Mu;
  double sigma2;
  
  arma::mat YZ;
  arma::mat Z; 
  
  arma::mat tmp;
  int T = Y.n_rows;
  
  //Rcpp::Rcout << 2 << std::endl;
  
  /////////////////////////////////////
  
  arma::mat X_;
  int p = 0;
  int Xflg = 0;
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    p = X_.n_cols;
    Xflg = 1;
  } 
  
  arma::mat H_;
  int m = 0;
  int Hflg = 0;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
    Hflg = 1;
  }
  
  
  double theta_;
  if (theta.isNotNull()) {
    theta_ = Rcpp::as<double>(theta);
  } else {
    theta_ = 1.0;
  }
  
  //Rcpp::Rcout << 3 << std::endl;
  
  /////////////////////////////////////
  // initialize the iteration
  arma::mat Yyj;
  
  if (updateYJ == 1) {
    Yyj = yeojohnsontr(Y, theta_, tol);
  } else {
    Yyj = Y;
  }
  
  //Rcpp::List iter = initGibbsRFLSMXcpp(Yyj, bset, tol, X, H, lambda2);
  Rcpp::List iter = simpleinitGibbsRFLSMXcpp(Yyj, bset, tol, X, H, lambda2);
  
  //Rcpp::Rcout << iter << std::endl;
  
  /////////////////////////////////////
  
  arma::mat Phimat(phiq, nsim);
  arma::mat Betamat;
  
  if (Xflg == 1) {
    Betamat.set_size(p, nsim);
  }
  
  arma::mat Taumat;
  arma::mat Gammamat;
  
  if (Hflg == 1) {
    Taumat.set_size(m, nsim);
    Gammamat.set_size(m, nsim);
  }
  
  arma::mat mu0mat(1, nsim);
  arma::mat Mumat(T, nsim);
  
  arma::mat eta2mat;
  if ((method == "LASSO") || (method == "ALASSO")) {
    eta2mat.set_size(phiq + p, nsim);
  }
  
  arma::mat sigma2mat(1, nsim);
  
  arma::mat lambda2mat;
  if ((updatelambda2 == 1) && ((method == "LASSO") || (method == "ALASSO"))) {
    lambda2mat.set_size(phiq + p, nsim);
  }
  
  arma::mat thetamat;
  if (updateYJ == 1) {
    thetamat.set_size(1, nsim);
  }
  
  arma::mat Zmat;
  
  if ((leftcensoring == 1) || (rounding == 1)) {
    Zmat.set_size(T, nsim);
  }
  
  //Rcpp::Rcout << 5 << std::endl;
  
  Z.set_size(T, 1); 
  Z.zeros();
  /////////////////////////////////////
  
  int tot_num = burnin + nsim * thin;
  int k = 0;
  for (int i = 0; i < tot_num; i++) {
    if (verbose == 1) {
      if (i % 100 == 0) {
        Rcpp::Rcout <<"Training: " << ((i + 0.0) / (tot_num + 0.0) * 100.0) << '%' << std::endl;
      }
    }
    
    Phi = Rcpp::as<arma::mat>(iter["Phi"]);
    Mu = Rcpp::as<arma::mat>(iter["Mu"]);
    sigma2 = iter["sigma2"];
    
    //Rcpp::Rcout << "i:" << i << std::endl;
    
    if ((leftcensoring == 1) || (rounding == 1)) {
      //YZ = getYZ(Yyj, Y, Phi, Mu, sigma2, 
      //           theta_, tol, leftcensoring, rounding, updateYJ);
      
      //Rcpp::Rcout <<"leftcensoring: " << leftcensoring << std::endl;
      //Rcpp::Rcout <<"rounding: " << rounding << std::endl;
      
      //YZ = getYZMHX(Y, Phi, Mu, sigma2, theta_, Z,  
      //              leftcensoring, rounding, 
      //              0, 1, tol);
      
      Z = getZZ(Y, Phi, Mu,  sigma2,  theta_,
           Z,  leftcensoring,  rounding, 
                        0, 1, tol);
      
      //YZ = getYZZ(Y, Phi, Mu, sigma2, theta_, Z,  
      //              leftcensoring, rounding, 
      //              0, 1, tol);
      
      YZ = Y + Z;
    } else {
      YZ = Y;
    }
    
    //Rcpp::Rcout << "leftcensoring:" << leftcensoring << std::endl;
    //Rcpp::Rcout << "rounding:" << rounding << std::endl;
    
    //Rcpp::Rcout << "Y:" << Y << std::endl;
    //Rcpp::Rcout << "YZ:" << YZ << std::endl;
    
    //Rcpp::Rcout << 7 << std::endl;
    
    if (updateYJ == 1) {
      tmp = thetaYeoJohnsonMHX(YZ, Phi, Mu, sigma2, theta_, 0, 1, tol);
      theta_ = tmp(0);
      Yyj = yeojohnsontr(YZ, theta_, tol);
    } else {
      Yyj = YZ;
    }
    
    //Rcpp::Rcout << "Yyj:" << Yyj << std::endl;
    
    //Rcpp::Rcout << 8 << std::endl;
    
    //if (YJ == 1) {
      
    //} else {
    //  Yyj = YZ;
    //}
    
    //Rcpp::Rcout << 9 << std::endl;
    
    iter = GibbsRFLSMXUpdatecpp(Yyj, iter, bset, tol, X, H);
      
    //Rcpp::Rcout << 10 << std::endl;
      
    if (i >= burnin) {
      if (i % thin == 0) {
         Phimat.col(k) = Rcpp::as<arma::mat>(iter["Phi"]);
      
        if (Xflg == 1) {
          Betamat.col(k) = Rcpp::as<arma::mat>(iter["Beta"]);
        }
        
        if (Hflg == 1) {
          Taumat.col(k) = Rcpp::as<arma::mat>(iter["Tau"]);
          Gammamat.col(k) = Rcpp::as<arma::mat>(iter["Gamma"]);
        }
        
        mu0mat(k) = iter["mu0"];
        
        if ((method == "LASSO") || (method == "ALASSO")){
          eta2mat.col(k) = Rcpp::as<arma::mat>(iter["eta2"]);
        }
        
        sigma2mat(k) = iter["sigma2"];
        
        if ((updatelambda2 == 1) && ((method == "LASSO") || (method == "ALASSO"))){
          lambda2mat.col(k) = Rcpp::as<arma::mat>(iter["lambda2"]);
        }
        
        if (updateYJ == 1) {
          thetamat(k) = theta_;
        }
        
        if ((leftcensoring == 1) || (rounding == 1)) {
          Zmat.col(k) = Z;
        }
        
        k++;
      }

    }
    
  }
  
  ///////////////////////////
  if (verbose == 1) {
    Rcpp::Rcout <<"Training: 100%" << std::endl;
  }
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  if (verbose == 1) {
    Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  }
  ///////////////////////////
  
  Rcpp::List out;
  
  out = Rcpp::List::create(
    Rcpp::_["Phi"] = Phimat,
    Rcpp::_["Beta"] = Betamat,
    Rcpp::_["Tau"] = Taumat,
    Rcpp::_["Gamma"] = Gammamat,
    Rcpp::_["mu0"] = mu0mat,
    Rcpp::_["eta2"] = eta2mat,
    Rcpp::_["sigma2"] = sigma2mat,
    Rcpp::_["lambda2"] = lambda2mat, 
    Rcpp::_["theta"] = thetamat,
    Rcpp::_["Z"] = Zmat
  );
  
  return(out);
  
}



//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYyjXph1(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double sigma2) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyj.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(1, q);
 arma::mat VasPhi;
  
 arma::mat fit(T, 1);
 arma::mat simYyjph1(T, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  int ii = 0;
  for (int i = 0; i < T; i++) {
    V(i, 0) = Yyj(i) - Mu(i, 0);
    
    //Rcpp::Rcout << i << std::endl;
    
    for (int j = 0; j < q; j++) {
        ii = i - 1 - j;
        if (ii >= 0) {
          Vas(0, j) = V(ii, 0);
        } else {
          Vas(0, j) = 0;
        }
    }
    
    VasPhi = Vas * Phi;
    fit(i, 0) = Mu(i, 0) + VasPhi(0);
    simYyjph1(i, 0) = R::rnorm(fit(i, 0), sqrt(sigma2));
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(simYyjph1);
  
}


//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYXph1(arma::colvec Y, arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps, 
                      int leftcensoring, int rounding, 
                      Rcpp::Nullable<Rcpp::NumericMatrix> Z = R_NilValue) {
  
  double lowerbound = 0.0;
  
  arma::mat Z_; 
  if (Z.isNotNull()) {
    Z_ = Rcpp::as<arma::mat>(Z);
  }
 
  arma::mat YZ; 
  if ((leftcensoring == 1) || (rounding == 1)) {
    YZ = Y + Z_;
  } else {
    YZ = Y;
  }
  
  arma::mat Yyj = yeojohnsontr(YZ, theta, eps);
  
   arma::colvec Yyjph1 = simYyjXph1(Yyj, Phi, Mu, sigma2); 
   arma::colvec Yph1 = invyeojohnsontr(Yyjph1, theta, eps); 
   arma::uvec ind0;
   if (leftcensoring == 1) {
     ind0 =arma::find(Yph1 <= lowerbound); 
     Yph1(ind0).fill(lowerbound);
   }
   
   if (rounding == 1) {
     Yph1 = arma::round(Yph1);
   }
    return(Yph1);
  
}


//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYyjXph2(int h,arma::colvec Yyjph1,arma::mat Phi,arma::mat Mu, double sigma2) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyjph1.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T + h, 1); 
 arma::mat Vas(1, q);
 arma::mat VasPhi;
  
 arma::mat fit(T + h, 1);
 arma::mat simYyjph2(T + h, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  int ii = 0;
  //for (int i = 0; i < (T + h); i++) {
  for (int i = 0; i < (T + h); i++) {
    
    
    //Rcpp::Rcout << i << std::endl;
    
    for (int j = 0; j < q; j++) {
        ii = i - 1 - j;
        if (ii >= 0) {
          Vas(0, j) = V(ii, 0);
        } else {
          Vas(0, j) = 0;
        }
      }
    
      VasPhi = Vas * Phi;
      fit(i, 0) = Mu(i, 0) + VasPhi(0);
      simYyjph2(i, 0) = R::rnorm(fit(i, 0), sqrt(sigma2));
    
    
    if (i < T) {
      V(i, 0) = Yyjph1(i) - Mu(i, 0);
    } else if (i >= T) {
      V(i, 0) = simYyjph2(i, 0) - Mu(i, 0);
    }
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(simYyjph2);

}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec simYXph2(int h,arma::colvec Y1,arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps, 
                      int leftcensoring, int rounding, 
                      Rcpp::Nullable<Rcpp::NumericMatrix> Z1 = R_NilValue) {
  
  double lowerbound = 0.0;
  
  arma::mat Z1_; 
  if (Z1.isNotNull()) {
    Z1_ = Rcpp::as<arma::mat>(Z1);
  }
 
  arma::mat YZ1; 
  if ((leftcensoring == 1) || (rounding == 1)) {
    YZ1 = Y1 + Z1_;
  } else {
    YZ1 = Y1;
  }
  
  arma::mat Yyjph1 = yeojohnsontr(YZ1, theta, eps);
  
 arma::colvec Yyjph2 = simYyjXph2(h, Yyjph1, Phi, Mu, sigma2); 
 arma::colvec Yph2 = invyeojohnsontr(Yyjph2, theta, eps); 
 arma::uvec ind0;
 if (leftcensoring == 1) {
   ind0 =arma::find(Yph2 <= lowerbound); 
   Yph2(ind0).fill(lowerbound);
 }
 
 if (rounding == 1) {
   Yph2 = arma::round(Yph2);
 }
  
  return(Yph2);
  
}

