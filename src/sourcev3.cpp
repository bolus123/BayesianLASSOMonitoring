#include <RcppArmadillo.h>   
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec getuDelta(arma::mat U, arma::vec Delta, int n) {
  arma::vec uDelta(n);
  arma::mat tmp(1, 1);
  
  for (int i = 0; i < n; i++) {
    tmp = U.row(i) * Delta;
    if (i > 0) {
      uDelta(i) = uDelta(i - 1) + tmp(0, 0);
    } else if (i == 0) {
      uDelta(i) = tmp(0, 0);
    }
  }
  
  return(uDelta);
}

// [[Rcpp::export]]
double getbeta00(arma::vec V, arma::vec Beta10, arma::vec Beta20, arma::vec UDelta, arma::mat U,
                 arma::mat X, arma::mat T, double tau200, double sigma2) {
  
  int n = V.n_elem;
  arma::vec R00(n);
  
  double mean = 0;
  double var = 0;
  
  double out = 0;
  
  R00 = V - UDelta - X * Beta10 - T * Beta20;
  
  mean = 1 / (n + 1 / tau200 ) * arma::accu(R00);
  var = sigma2 / (n + 1 / tau200);
  
  out = R::rnorm(mean, sqrt(var));
  
  return(out);
  
}

// [[Rcpp::export]]
arma::mat getInvTau2(arma::vec Tau2) {
  
  int q = Tau2.n_elem;
  arma::mat out(q, q);
  out.zeros();
  
  for (int i = 0; i < q; i++) {
    
    out(i, i) = 1 / Tau2(i);
    
  }
  
  return(out);
  
}

// [[Rcpp::export]]
arma::vec rmvnormCpp(arma::vec Mu, arma::mat Sigma) {
  
  int m = Mu.n_elem;
  arma::mat SigmaSq(m, m);
  arma::vec Z(m);
  arma::vec out(m);
  
  arma::sqrtmat_sympd(SigmaSq, Sigma);
  
  //if (SigmaSq == false)  {
  //  arma::sqrtmat(SigmaSq, Sigma);
  //}
  
  for (int i = 0; i < m; i++) {
    Z(i) = R::rnorm(0, 1);
  }
  
  out = SigmaSq * Z + Mu;
  
  return(out);
  
}

// [[Rcpp::export]]
arma::vec getBeta00(double beta00, int n) {
  
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    out(i) = beta00;
  }

  return(out);
  
}

// [[Rcpp::export]]
arma::vec getBeta10(arma::vec V, arma::vec Beta00, arma::vec Beta20, arma::vec UDelta, 
                 arma::mat U, arma::mat X, arma::mat T, arma::vec tau210, double sigma2) {
  
  int n = V.n_elem;
  int q = X.n_cols;
  
  arma::vec R10(n);

  arma::mat InvTau210(q, q) ;
  
  arma::mat XXplusInvTau210(q, q);
  arma::mat InvXXplusInvTau210(q, q);
  
  arma::vec Mean(q);
  arma::mat Var(q, q);
  
  arma::vec out(q);
  
  R10 = V - Beta00 - UDelta - T * Beta20;
  
  InvTau210 = getInvTau2(tau210); 
  
  XXplusInvTau210 = arma::trans(X) * X + InvTau210;
  
  inv_sympd(InvXXplusInvTau210, XXplusInvTau210);
  
  Mean = InvXXplusInvTau210 * arma::trans(X) * R10;
  Var = sigma2 * InvXXplusInvTau210;
  
  out = rmvnormCpp(Mean, Var);
  
  return(out);
  
}

// [[Rcpp::export]]
double rtrnormCpp(double mean, double var, double lower, double upper) {
  double tmp = R::runif(0, 1);
  double sd = sqrt(var);
  double a = (lower - mean) / sd;
  double b = (upper - mean) / sd;
  double Z = R::pnorm5(b, 0, 1, true, false) - R::pnorm5(a, 0, 1, true, false);
  double out;
  out = R::qnorm5(tmp * Z + R::pnorm5(a, 0, 1, true, false), 0, 1, true, false) * sd + mean;
  return(out);
}


// [[Rcpp::export]]
arma::vec getBeta20(arma::vec V, arma::vec Beta00, arma::vec Beta10, arma::vec UDelta, 
                    arma::mat U, arma::mat X, arma::mat T, arma::vec tau220, double sigma2) {
  
  int n = V.n_elem;
  int p = T.n_cols;
  
  arma::vec R20(n);
  
  arma::mat InvTau220(p, p) ;
  
  arma::mat TTplusInvTau220(p, p);
  arma::mat InvTTplusInvTau220(p, p);
  
  arma::vec Mean(p);
  arma::mat Var(p, p);
  
  double mean;
  double var;
  arma::vec tmpmean;
  arma::vec tmpvar;
  
  arma::colvec SubVecMean;
  arma::rowvec SubVecVar;
  arma::mat SubMatVar;
  arma::mat InvSubMatVar;
  
  arma::vec out(p);
  arma::colvec SubOut;
  
  R20 = V - Beta00 - UDelta - X * Beta10;
  
  InvTau220 = getInvTau2(tau220); 
  
  TTplusInvTau220 = arma::trans(T) * T + InvTau220;
  
  inv_sympd(InvTTplusInvTau220, TTplusInvTau220);
  
  Mean = InvTTplusInvTau220 * arma::trans(T) * R20;
  Var = sigma2 * InvTTplusInvTau220;
  
  for (int i = 0; i < p; i++) {
    
    if (i > 0) {
      SubVecMean = Mean.subvec(0, i - 1);
      SubVecVar = Var.submat(i, 0, i, i - 1);
      SubMatVar = Var.submat(0, 0, i - 1, i - 1);
      inv_sympd(InvSubMatVar, SubMatVar);
      SubOut = out.subvec(0, i - 1);
      
      tmpmean = Mean(i) + SubVecVar * InvSubMatVar * (SubOut - SubVecMean);
      tmpvar = Var(i, i) - SubVecVar * InvSubMatVar * arma::trans(SubVecVar);
      
      mean = tmpmean(0);
      var = tmpvar(0);
      out(i) = rtrnormCpp(mean, var, -abs(out(i - 1)), abs(out(i - 1)));
    } else {
      mean = Mean(i);
      var = Var(i, i);
      out(i) = R::rnorm(mean, sqrt(var));
    }
    
  }
  
  return(out);
  
}

// [[Rcpp::export]]
arma::mat getSMat(int n) {
  arma::mat tmp(n, n);
  arma::mat out(n, n);
  tmp.ones();
  out = arma::trimatu(tmp);
  return(out);
}

// [[Rcpp::export]]
arma::mat getSSMat(int n) {
  int tmp;
  arma::mat out(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i > j) {
        tmp = i;
      } else {
        tmp = j;
      }
      out(i, j) = n - tmp;
    }
  }
  return(out);
}

// [[Rcpp::export]]
arma::vec getDelta(arma::vec V, arma::vec Beta00, arma::vec Beta10, arma::vec Beta20, 
                   arma::mat U, arma::mat X, arma::mat T, arma::vec tauDelta, double sigma2) {
  
  int n = V.n_elem;
  int K = U.n_cols;
  
  arma::vec RDelta(n);
  
  arma::mat InvTauDelta(K, K) ;
  
  arma::mat SMat(n, n);
  arma::mat SSMat(n, n);
  
  arma::mat USSUplusInvTauDelta(K, K);
  arma::mat InvUSSUplusInvTauDelta(K, K);
  
  arma::vec Mean(K);
  arma::mat Var(K, K);
  
  arma::vec out(K);
  
  RDelta = V - Beta00 - X * Beta10 - T * Beta20;
  
  InvTauDelta = getInvTau2(tauDelta); 
  
  SMat = getSMat(n);
  SSMat = getSSMat(n);
  
  USSUplusInvTauDelta = arma::trans(U) * SSMat * U + InvTauDelta;
  
  inv_sympd(InvUSSUplusInvTauDelta, USSUplusInvTauDelta);
  
  Mean = InvUSSUplusInvTauDelta * arma::trans(U) * SMat * RDelta;
  Var = sigma2 * InvUSSUplusInvTauDelta;
  
  out = rmvnormCpp(Mean, Var);
  
  return(out);
}

// [[Rcpp::export]]
arma::mat getU(arma::vec V, double beta00, arma::vec Beta10, arma::vec Beta20, arma::vec Delta,
               arma::vec UDelta, arma::mat X, arma::mat T, double sigma2, arma::vec gamma) {
  
  int n = V.n_elem;
  int K = gamma.n_elem;
  
  arma::vec R(n);
  arma::vec tmpR;
  
  arma::mat tmp(n, K);
  arma::mat prob(n, K);
  double randProb;
  double cursor;
  arma::vec tmpprobcumsum(K);
  arma::mat out(n, K);
  out.zeros();
  
  for (int i = 0; i < n; i++) {
    if (i > 0) {
      tmpR = V(i) - beta00 - UDelta(i - 1) - X.row(i) * Beta10 - T.row(i) * Beta20;
    } else {
      tmpR = V(i) - beta00 - X.row(i) * Beta10 - T.row(i) * Beta20;
    }
    R(i) = tmpR(0);
    tmp.row(i) = arma::exp(- 1 / 2 / sigma2 * (Delta % Delta - 2 * Delta * R(i)) + arma::log(gamma));
    prob.row(i) = tmp.row(i) / arma::accu(tmp.row(i));
    tmpprobcumsum = arma::cumsum(prob.row(i));
    randProb = R::runif(0, 1);
    cursor = 0;
    for (int k = 0; k < K; k++) {
      if (randProb > tmpprobcumsum(k)) {
        cursor = cursor + 1;
      }
    }
    out(i, cursor) = 1;
  }
  return(out);
}

// [[Rcpp::export]]
double rinvgammaCpp(double shape, double rate) {
  
  double scale = 1 / rate;
  double tmp = R::rgamma(shape, scale);
  double out = 1 / tmp;
  return(out);
  
}

// [[Rcpp::export]]
double getSigma2(arma::vec V, double beta00, arma::vec Beta00, arma::vec Beta10, arma::vec Beta20, arma::vec Delta, 
                 arma::vec UDelta, arma::mat X, arma::mat T, double tau200, arma::vec tau210, 
                 arma::vec tau220, arma::vec tauDelta, double a1, double a2) {
  
  int n = V.n_elem;
  int q = Beta10.n_elem;
  int p = Beta20.n_elem;
  int K = Delta.n_elem;
  
  arma::vec resi(n);
  
  double shape = (n + 1 + q + p + K) / 2 + a1;
  arma::vec tmprate;
  double rate;
  double out;
  
  arma::mat InvTau210 = getInvTau2(tau210);
  arma::mat InvTau220 = getInvTau2(tau220);
  arma::mat InvTauDelta = getInvTau2(tauDelta);
  
  resi = V - Beta00 - UDelta - X * Beta10 - T * Beta20;
  tmprate = arma::accu(resi % resi) / 2 + (beta00 * beta00) / tau200 + 
    arma::trans(Beta10) * InvTau210 * Beta10 + arma::trans(Beta20) * InvTau220 * Beta20 +
    arma::trans(Delta) * InvTauDelta * Delta + 2 * a2;
  rate = tmprate(0);
  
  out = rinvgammaCpp(shape, rate);
  return(out);
  
}

// [[Rcpp::export]]
double rgammaCpp(double shape, double rate) {
  double scale = 1 / rate;
  double out = R::rgamma(shape, scale);
  return(out);
}

// [[Rcpp::export]]
double getLambda2(double tau200, arma::vec tau210, 
                  arma::vec tau220, arma::vec tauDelta, double b1, double b2) {
  int q = tau210.n_elem;
  int p = tau220.n_elem;
  int K = tauDelta.n_elem;
  
  double shape = K + p + q + 1 + b1;
  double rate = 1 / 2 * (tau200 * tau200 + arma::accu(tau210 % tau210) + arma::accu(tau220 % tau220) + 
    arma::accu(tauDelta % tauDelta) + 2 / b2);
  
  double out = rgammaCpp(shape, rate);
  return(out);
}

// [[Rcpp::export]]
arma::colvec chol_solve(arma::mat& M, arma::colvec& V) {
  arma::mat R = arma::chol(M);
  arma::colvec b_star = arma::solve(trimatl(R.t()), V);
  return(arma::solve(trimatu(R), b_star));
}

// //' Randomly generate a generalized inverse gaussian random variable.
// //'
// //' Randomly generates one draw from a generalized inverse gaussian distribution.
// //' @param chi A positive double.
// //' @param psi A positive double.
// //' @param lambda A non-negative double.
// //' @return A random draw from the generalized inverse gaussian distribution with parameters chi, psi, and lambda (double).
// [[Rcpp::export]]
double rgig_cpp(double chi, double psi, double lambda) {
  double final_draw = 0;
  double alpha = sqrt(psi / chi);  //psi = b, chi = a, lambda = p
  double beta = sqrt(chi*psi);
  if ((lambda > 1) || (beta > 1)) {
    double m = (sqrt(pow(lambda - 1.0, 2) + pow(beta, 2)) + (lambda - 1.0)) / beta;
    double a = -2.0*(lambda + 1.0) / beta - m;
    double b = 2.0*(lambda - 1.0)*m / beta - 1.0;
    double c = m;
    double p = b - pow(a, 2) / 3.0;
    double q = 2.0*pow(a, 3) / 27.0 - a*b / 3.0 + c;
    double phi = acos(-(q / 2.0)*sqrt(-27.0 / pow(p, 3)));
    double x_minus = sqrt(-(4.0 / 3.0)*p)*cos(phi / 3.0 + (4.0 / 3.0)*M_PI) - a / 3.0;
    double x_plus = sqrt(-(4.0 / 3.0)*p)*cos(phi / 3.0) - a / 3.0;
    double v_plus = sqrt(pow(m, lambda - 1.0)*exp(-(beta / 2.0)*(m + 1.0 / m)));
    double u_minus = (x_minus - m)*sqrt(pow(x_minus, lambda - 1.0)*exp(-(beta / 2.0)*(x_minus + 1.0 / x_minus)));
    double u_plus = (x_plus - m)*sqrt(pow(x_plus, lambda - 1.0)*exp(-(beta / 2.0)*(x_plus + 1.0 / x_plus)));
    bool keep_looping = true;
    double u_draw; double v_draw; double x_draw;
    while (keep_looping) {
      u_draw = R::runif(u_minus, u_plus);
      v_draw = R::runif(0, v_plus);
      x_draw = u_draw / v_draw + m;
      if ((pow(v_draw, 2) <= pow(x_draw, lambda - 1.0)*exp(-(beta / 2.0)*(x_draw + 1.0 / x_draw))) && (x_draw > 0)) {
        final_draw = x_draw;
        keep_looping = false;
      }
    }
  }
  else if (lambda >= 0 && lambda <= 1 && beta >= std::min(1.0 / 2.0, (2.0 / 3.0)*sqrt(1.0 - lambda)) && beta <= 1) {
    double m = beta / ((1.0 - lambda) + sqrt(pow(1.0 - lambda, 2) + pow(beta, 2)));
    double x_plus = ((1.0 + lambda) + sqrt(pow(1 + lambda, 2) + pow(beta, 2))) / beta;
    double v_plus = sqrt(pow(m, lambda - 1.0)*exp(-(beta / 2.0)*(m + 1.0 / m)));
    double u_plus = x_plus*sqrt(pow(x_plus, lambda - 1.0)*exp(-(beta / 2.0)*(x_plus + 1.0 / x_plus)));
    bool keep_looping = true;
    double u_draw; double v_draw; double x_draw;
    while (keep_looping) {
      u_draw = R::runif(0, u_plus);
      v_draw = R::runif(0, v_plus);
      x_draw = u_draw / v_draw;
      if (pow(v_draw, 2) <= pow(x_draw, lambda - 1.0)*exp(-(beta / 2.0)*(x_draw + 1.0 / x_draw))) {
        final_draw = x_draw;
        keep_looping = false;
      }
    }
  }
  else if (lambda >= 0 && lambda < 1 && beta > 0 && beta <= (2.0 / 3.0)*sqrt(1.0 - lambda)) {
    double m = beta / ((1.0 - lambda) + sqrt(pow(1.0 - lambda, 2) + pow(beta, 2)));
    double x0 = beta / (1.0 - lambda);
    double x_star = std::max(x0, 2.0 / beta);
    double k1 = pow(m, lambda - 1.0)*exp(-(beta / 2.0)*(m + 1.0 / m));
    double A1 = k1*x0;
    double A2; double k2;
    if (x0 < 2.0 / beta) {
      k2 = exp(-beta);
      if (lambda == 0) {
        A2 = k2*log(2.0 / pow(beta, 2));
      }
      else {
        A2 = k2*(pow(2.0 / beta, lambda) - pow(x0, lambda)) / lambda;
      }
    }
    else {
      k2 = 0;
      A2 = 0;
    }
    double k3 = pow(x_star, lambda - 1.0);
    double A3 = 2.0*k3*exp(-x_star*beta / 2.0) / beta;
    double A = A1 + A2 + A3;
    bool keep_looping = true;
    double u_draw; double v_draw; double x_draw; double h;
    while (keep_looping) {
      u_draw = R::runif(0, 1);
      v_draw = R::runif(0, A);
      if (v_draw <= A1) {
        x_draw = x0*v_draw / A1;
        h = k1;
      }
      else if (v_draw <= A1 + A2) {
        v_draw = v_draw - A1;
        if (lambda == 0) {
          x_draw = beta*exp(v_draw*exp(beta));
        }
        else {
          x_draw = pow(pow(x0, lambda) + v_draw*lambda / k2, 1.0 / lambda);
        }
        h = k2*pow(x_draw, lambda - 1.0);
      }
      else {
        v_draw = v_draw - (A1 + A2);
        x_draw = -2.0 / beta*log(exp(-x_star*beta / 2.0) - v_draw*beta / (2.0*k3));
        h = k3*exp(-x_draw*beta / 2.0);
      }
      if (u_draw*h <= pow(x_draw, lambda - 1.0)*exp(-(beta / 2.0)*(x_draw + 1.0 / x_draw))) {
        final_draw = x_draw;
        keep_looping = false;
      }
    }
  }
  return final_draw / alpha;
}

// [[Rcpp::export]]
double getTau2(double beta, double sigma2, double lambda2) {
  double out = rgig_cpp(1 / 2, lambda2, beta * beta / sigma2);
  return(out);
}
