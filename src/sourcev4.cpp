#include <RcppArmadillo.h>   
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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
arma::vec getBetaNonMonotonicity(arma::vec zeta, arma::mat X,
                    arma::vec tau2, double sigma2) {
  
  int q = X.n_cols;
  
  //arma::vec R10(n);
  
  arma::mat InvTau2(q, q) ;
  
  arma::mat XXplusInvTau2(q, q);
  arma::mat InvXXplusInvTau2(q, q);
  
  arma::vec Mean(q);
  arma::mat Var(q, q);
  
  arma::vec out(q);
  
  //R10 = V - Beta00 - UDelta - T * Beta20;
  
  InvTau2 = getInvTau2(tau2); 
  
  XXplusInvTau2 = arma::trans(X) * X + InvTau2;
  
  inv_sympd(InvXXplusInvTau2, XXplusInvTau2);
  
  Mean = InvXXplusInvTau2 * arma::trans(X) * zeta;
  Var = sigma2 * InvXXplusInvTau2;
  
  out = rmvnormCpp(Mean, Var);
  
  return(out);
  
}

// [[Rcpp::export]]
double rtrnormCpp(double mean, double var, double lower, double upper) {
  double tmp = R::runif(0, 1);
  double sd = sqrt(var);
  double a = (lower - mean) / sd;
  double b = (upper - mean) / sd;
  
  double z1 = R::pnorm5(b, 0, 1, true, false);
  double z2 = R::pnorm5(a, 0, 1, true, false);
  
  double Z;
  double out;
  if (z1 != z2) {
    Z = z1 - z2;
    out = R::qnorm5(tmp * Z + z2, 0, 1, true, false) * sd + mean;
  } else {
      z1 = R::plnorm(exp(b), 0, 1, true, false);
      z2 = R::plnorm(exp(a), 0, 1, true, false);
      
      if (z1 != z2) {
        Z = z1 - z2;
        out = log(R::qlnorm(tmp * Z + z2, 0, 1, true, false)) * sd + mean;
      } else {
        out = R::runif(lower, upper);
      }
  }
  
  
  return(out);
}


// [[Rcpp::export]]
arma::vec getBetaMonotonicity(arma::vec zeta, 
                              arma::mat T, 
                              arma::vec tau2, 
                              double sigma2) {
  
  int p = T.n_cols;
  
  //arma::vec R20(n);
  
  arma::mat InvTau2(p, p) ;
  
  arma::mat TTplusInvTau2(p, p);
  arma::mat InvTTplusInvTau2(p, p);
  
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
  
  //R20 = V - Beta00 - UDelta - X * Beta10;
  
  InvTau2 = getInvTau2(tau2); 
  
  TTplusInvTau2 = arma::trans(T) * T + InvTau2;
  
  inv_sympd(InvTTplusInvTau2, TTplusInvTau2);
  
  Mean = InvTTplusInvTau2 * arma::trans(T) * zeta;
  Var = sigma2 * InvTTplusInvTau2;
  
  //Rcpp::Rcout << Mean << std::endl;
  //Rcpp::Rcout << Var << std::endl;
  
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
      
      //Rcpp::Rcout << "i:" << i << std::endl;
      //Rcpp::Rcout << mean << std::endl;
      //Rcpp::Rcout << var << std::endl;
      
      out(i) = rtrnormCpp(mean, var, -abs(out(i - 1)), abs(out(i - 1)));
    } else {
      mean = Mean(i);
      var = Var(i, i);
      out(i) = R::rnorm(mean, sqrt(var));
    }
    
  }
  
  //Rcpp::Rcout << "out:" << out << std::endl;
  
  return(out);
  
}

// [[Rcpp::export]]
arma::mat getUW(arma::mat U, arma::mat W) {
  int n = U.n_rows;
  int K = U.n_cols;
  int r = W.n_cols;
  arma::mat UW(n, K * r);
  
  int k = 0;
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < r; j++) {
      UW.col(k) = U.col(i) * W.col(j);
      k++;
    }
  }
  
  return(UW);
}

// [[Rcpp::export]]
arma::vec getTheta(arma::vec WRow, arma::vec delta0, arma::vec delta1, int K) {
  int r = WRow.n_elem;
  arma::vec out(K);
  arma::vec tmp(1);
  arma::vec tmpDelta1(r);
  int indx = 0;
  for (int i = 0; i < K; i++) {
    tmpDelta1 = delta1.subvec(indx, indx + r - 1);
    tmp = delta0(i) + WRow * tmpDelta1;
    out(i) = tmp(0);
    indx = indx + r;
  }
  return(out);
}

// [[Rcpp::export]]
arma::vec getEta(double zetaElem, arma::vec theta, 
                 arma::vec gamma, double sigma2) {
  arma::vec eta = arma::exp(
    - 1.0 / 2.0 / sigma2 * (theta % theta - 2.0 * theta * zetaElem) + arma::log(gamma)
  );
  return(eta);
}

// [[Rcpp::export]]
Rcpp::List getUWithW(arma::vec zetadelta, arma::mat W, int K, 
                arma::vec delta0, arma::vec delta1, double sigma2, 
                arma::vec gamma) {
  
  int n = zetadelta.n_elem;
  arma::vec theta;
  double zetaElem;
  arma::vec WRow;
  arma::vec eta;
  arma::mat prob(n, K);
  arma::mat out(n, K);
  out.zeros();
  
  arma::rowvec tmpProb;

  int i;
  
  for (i = 0; i < n; i++) {
    WRow = W.row(i);
    theta = getTheta(WRow, delta0, delta1, K);
    zetaElem = zetadelta(i);
    eta = getEta(zetaElem, theta, gamma, sigma2);
    tmpProb = arma::trans(eta / arma::accu(eta));
    prob.row(i) = tmpProb;
  }
  
  arma::mat tmpcumsumProb;
  tmpcumsumProb = arma::cumsum(prob, 1);
  double randProb;
  int cursor;
  
  for (i = 0; i < n; i++) {
    randProb = R::runif(0, 1);
    cursor = 0;
    for (int k = 0; k < K; k++) {
      if (randProb > tmpcumsumProb(i, k)) {
        cursor = cursor + 1;
      }
    }
    out(i, cursor) = 1;
  }
  
  Rcpp::List outList;
  outList = Rcpp::List::create(
    Rcpp::_["prob"] = prob,
    Rcpp::_["U"] = out
  );
  return(outList);
}

// [[Rcpp::export]]
Rcpp::List getUWithoutW(arma::vec zetadelta, int K, 
                     arma::vec delta0, double sigma2, 
                     arma::vec gamma) {
  
  //Rcpp::Rcout << "a" << std::endl;
  
  int n = zetadelta.n_elem;
  arma::vec theta;
  double zetaElem;
  arma::vec eta;
  arma::mat prob(n, K);
  arma::mat out(n, K);
  out.zeros();
  
  arma::rowvec tmpProb;

  int i;
  
  //Rcpp::Rcout << "b" << std::endl;
  
  for (i = 0; i < n; i++) {
    theta = delta0;
    zetaElem = zetadelta(i);
    eta = getEta(zetaElem, theta, gamma, sigma2);
    tmpProb = arma::trans(eta / arma::accu(eta));
    prob.row(i) = tmpProb;
  }
  
  //Rcpp::Rcout << prob << std::endl;
  
  //Rcpp::Rcout << "c" << std::endl;
  
  arma::mat tmpcumsumProb;
  tmpcumsumProb = arma::cumsum(prob, 1);
  
  //Rcpp::Rcout << tmpcumsumProb << std::endl;
  
  double randProb;
  int cursor;
  
  //Rcpp::Rcout << "d" << std::endl;
  
  for (i = 0; i < n; i++) {
    randProb = R::runif(0, 1);
    cursor = 0;
    for (int k = 0; k < K; k++) {
      if (randProb > tmpcumsumProb(i, k)) {
        cursor = cursor + 1;
      }
    }
    out(i, cursor) = 1;
  }
  
  //Rcpp::Rcout << "e" << std::endl;
  
  Rcpp::List outList;
  outList = Rcpp::List::create(
    Rcpp::_["prob"] = prob,
    Rcpp::_["U"] = out
  );
  return(outList);
}



// [[Rcpp::export]]
double rinvgammaCpp(double shape, double rate) {
  
  double scale = 1 / rate;
  double tmp = R::rgamma(shape, scale);
  double out = 1 / tmp;
  return(out);
  
}

// [[Rcpp::export]]
double getSigma2(arma::vec resi,
                 arma::vec BetaDelta,
                 arma::vec Tau2,
                 double a1, double a2) {
  
  //Rcpp::Rcout << "sigma2" << std::endl;
  
  int n = resi.n_elem;
  //int q = Beta10.n_elem;
  //int p = Beta20.n_elem;
  //int K = Delta.n_elem;
  
  int m = BetaDelta.n_elem;
  
  //arma::vec resi(n);
  
  //double shape = (n + 1 + q + p + K) / 2 + a1;
  double shape = (n + m) / 2.0 + a1;
  arma::vec tmprate;
  double rate;
  double out;
  
  arma::mat InvTau2 = getInvTau2(Tau2);
  //arma::mat InvTau210 = getInvTau2(tau210);
  //arma::mat InvTau220 = getInvTau2(tau220);
  //arma::mat InvTauDelta = getInvTau2(tauDelta);
  
  //resi = V - Beta00 - UDelta - X * Beta10 - T * Beta20;
  //tmprate = arma::accu(resi % resi) / 2 + (beta00 * beta00) / tau200 + 
  //  arma::trans(Beta10) * InvTau210 * Beta10 + arma::trans(Beta20) * InvTau220 * Beta20 +
  //  arma::trans(Delta) * InvTauDelta * Delta + 2 * a2;
  //Rcpp::Rcout << "resi:" << resi << std::endl;
  //Rcpp::Rcout << "BetaDelta:" << BetaDelta << std::endl;
  tmprate = (arma::trans(resi) * resi + arma::trans(BetaDelta) * InvTau2 * BetaDelta) / 2.0 + a2;
  //Rcpp::Rcout << "tmprate:" << tmprate << std::endl;
  rate = tmprate(0);
  
  
  //Rcpp::Rcout << "shape:" << shape << std::endl;
  //Rcpp::Rcout << "rate:" << rate << std::endl;
  
  out = rinvgammaCpp(shape, rate);
  return(out);
  
}

// [[Rcpp::export]]
double rgammaCpp(double shape, double rate) {
  double scale = 1 / rate;
  double out = R::rgamma(shape, scale);
  return(out);
}

//// [[Rcpp::export]]
//double getLambda2(double tau200, arma::vec tau210, 
//                  arma::vec tau220, arma::vec tauDelta, double b1, double b2) {
//  int q = tau210.n_elem;
//  int p = tau220.n_elem;
//  int K = tauDelta.n_elem;
//  
//  double shape = K + p + q + 1 + b1;
//  double rate = 1 / 2 * (tau200 * tau200 + arma::accu(tau210 % tau210) + arma::accu(tau220 % tau220) + 
//    arma::accu(tauDelta % tauDelta) + 2 / b2);
//  
//  double out = rgammaCpp(shape, rate);
//  return(out);
//}

// [[Rcpp::export]]
double getLambda2(arma::vec Tau2, double b1, double b2) {
  double m = Tau2.n_elem;
  
  double shape = m + b1;
  double rate = (arma::accu(Tau2 % Tau2)) / 2.0 + 1.0 / b2;
  
  //Rcpp::Rcout << "lambda2" << std::endl;
  //Rcpp::Rcout << "shape:" << shape << std::endl;
  //Rcpp::Rcout << "rate:" << rate << std::endl;
  
  double out = rgammaCpp(shape, rate);
  return(out);
}

// [[Rcpp::export]]
double getLambda2EM(arma::vec ExpectedTau2) {
  int m = ExpectedTau2.n_elem;
  double out = sqrt(2.0 * m / arma::accu(ExpectedTau2));
  return(out);
}

// [[Rcpp::export]]
double rgig_cpp(double chi, double psi, double lambda) {
  // Extract R's optim function
  Rcpp::Environment GIGrvg("package:GIGrvg"); 
  Rcpp::Function rgig = GIGrvg["rgig"];
  
  // Call the optim function from R in C++ 
  Rcpp::List out = rgig(Rcpp::_["n"]    = 1,
                                 Rcpp::_["lambda"]     = lambda,
                                 Rcpp::_["chi"] = chi,
                                 Rcpp::_["psi"] = psi);
  
  // Return estimated values
  return out(0);
}



////// //' Randomly generate a generalized inverse gaussian random variable.
////// //'
////// //' Randomly generates one draw from a generalized inverse gaussian distribution.
////// //' @param chi A positive double.
////// //' @param psi A positive double.
////// //' @param lambda A non-negative double.
////// //' @return A random draw from the generalized inverse gaussian distribution with parameters chi, psi, and lambda (double).
////// [[Rcpp::export]]
////double rgig_cpp(double chi, double psi, double lambda) {
////  double final_draw = 0;
////  double alpha = sqrt(psi / chi);  //psi = b, chi = a, lambda = p
////  double beta = sqrt(chi*psi);
////  if ((lambda > 1) || (beta > 1)) {
////    double m = (sqrt(pow(lambda - 1.0, 2) + pow(beta, 2)) + (lambda - 1.0)) / beta;
////    double a = -2.0*(lambda + 1.0) / beta - m;
////    double b = 2.0*(lambda - 1.0)*m / beta - 1.0;
////    double c = m;
////    double p = b - pow(a, 2) / 3.0;
////    double q = 2.0*pow(a, 3) / 27.0 - a*b / 3.0 + c;
////    double phi = acos(-(q / 2.0)*sqrt(-27.0 / pow(p, 3)));
////    double x_minus = sqrt(-(4.0 / 3.0)*p)*cos(phi / 3.0 + (4.0 / 3.0)*M_PI) - a / 3.0;
////    double x_plus = sqrt(-(4.0 / 3.0)*p)*cos(phi / 3.0) - a / 3.0;
////    double v_plus = sqrt(pow(m, lambda - 1.0)*exp(-(beta / 2.0)*(m + 1.0 / m)));
////    double u_minus = (x_minus - m)*sqrt(pow(x_minus, lambda - 1.0)*exp(-(beta / 2.0)*(x_minus + 1.0 / x_minus)));
////    double u_plus = (x_plus - m)*sqrt(pow(x_plus, lambda - 1.0)*exp(-(beta / 2.0)*(x_plus + 1.0 / x_plus)));
////    bool keep_looping = true;
////    double u_draw; double v_draw; double x_draw;
////    while (keep_looping) {
////      u_draw = R::runif(u_minus, u_plus);
////      v_draw = R::runif(0, v_plus);
////      x_draw = u_draw / v_draw + m;
////      if ((pow(v_draw, 2) <= pow(x_draw, lambda - 1.0)*exp(-(beta / 2.0)*(x_draw + 1.0 / x_draw))) && (x_draw > 0)) {
////        final_draw = x_draw;
////        keep_looping = false;
////      }
////    }
////  }
////  else if (lambda >= 0 && lambda <= 1 && beta >= std::min(1.0 / 2.0, (2.0 / 3.0)*sqrt(1.0 - lambda)) && beta <= 1) {
////    double m = beta / ((1.0 - lambda) + sqrt(pow(1.0 - lambda, 2) + pow(beta, 2)));
////    double x_plus = ((1.0 + lambda) + sqrt(pow(1 + lambda, 2) + pow(beta, 2))) / beta;
////    double v_plus = sqrt(pow(m, lambda - 1.0)*exp(-(beta / 2.0)*(m + 1.0 / m)));
////    double u_plus = x_plus*sqrt(pow(x_plus, lambda - 1.0)*exp(-(beta / 2.0)*(x_plus + 1.0 / x_plus)));
////    bool keep_looping = true;
////    double u_draw; double v_draw; double x_draw;
////    while (keep_looping) {
////      u_draw = R::runif(0, u_plus);
////      v_draw = R::runif(0, v_plus);
////      x_draw = u_draw / v_draw;
////      if (pow(v_draw, 2) <= pow(x_draw, lambda - 1.0)*exp(-(beta / 2.0)*(x_draw + 1.0 / x_draw))) {
////        final_draw = x_draw;
////        keep_looping = false;
////      }
////    }
////  }
////  else if (lambda >= 0 && lambda < 1 && beta > 0 && beta <= (2.0 / 3.0)*sqrt(1.0 - lambda)) {
////    double m = beta / ((1.0 - lambda) + sqrt(pow(1.0 - lambda, 2) + pow(beta, 2)));
////    double x0 = beta / (1.0 - lambda);
////    double x_star = std::max(x0, 2.0 / beta);
////    double k1 = pow(m, lambda - 1.0)*exp(-(beta / 2.0)*(m + 1.0 / m));
////    double A1 = k1*x0;
////    double A2; double k2;
////    if (x0 < 2.0 / beta) {
////      k2 = exp(-beta);
////      if (lambda == 0) {
////        A2 = k2*log(2.0 / pow(beta, 2));
////      }
////      else {
////        A2 = k2*(pow(2.0 / beta, lambda) - pow(x0, lambda)) / lambda;
////      }
////    }
////    else {
////      k2 = 0;
////      A2 = 0;
////    }
////    double k3 = pow(x_star, lambda - 1.0);
////    double A3 = 2.0*k3*exp(-x_star*beta / 2.0) / beta;
////    double A = A1 + A2 + A3;
////    bool keep_looping = true;
////    double u_draw; double v_draw; double x_draw; double h;
////    while (keep_looping) {
////      u_draw = R::runif(0, 1);
////      v_draw = R::runif(0, A);
////      if (v_draw <= A1) {
////        x_draw = x0*v_draw / A1;
////        h = k1;
////      }
////      else if (v_draw <= A1 + A2) {
////        v_draw = v_draw - A1;
////        if (lambda == 0) {
////          x_draw = beta*exp(v_draw*exp(beta));
////        }
////        else {
////          x_draw = pow(pow(x0, lambda) + v_draw*lambda / k2, 1.0 / lambda);
////        }
////        h = k2*pow(x_draw, lambda - 1.0);
////      }
////      else {
////        v_draw = v_draw - (A1 + A2);
////        x_draw = -2.0 / beta*log(exp(-x_star*beta / 2.0) - v_draw*beta / (2.0*k3));
////        h = k3*exp(-x_draw*beta / 2.0);
////      }
////      if (u_draw*h <= pow(x_draw, lambda - 1.0)*exp(-(beta / 2.0)*(x_draw + 1.0 / x_draw))) {
////        final_draw = x_draw;
////        keep_looping = false;
////      }
////    }
////  }
////  return final_draw / alpha;
////}

//// [[Rcpp::export]]
//double getTau2(double beta, double sigma2, double lambda2) {
//  //rgig_cpp(double chi, double psi, double lambda)
//  double out = rgig_cpp(beta * beta / sigma2, lambda2, 1.0 / 2.0);
//  return(out);
//}

// [[Rcpp::export]]
arma::vec getTau2(arma::vec beta, double sigma2, double lambda2) {
  //rgig_cpp(double chi, double psi, double lambda)
  int m = beta.n_elem;
  arma::vec out(m);
  for (int i = 0; i < m; i++) {
    out(i) = rgig_cpp(beta(i) * beta(i) / sigma2, lambda2, 1.0 / 2.0);
  }
  return(out);
}

// [[Rcpp::export]]
double modifiedBesselfunction2nd(double x, double nu) {
  // Extract R's optim function
  Rcpp::Environment base("package:base"); 
  Rcpp::Function besselK = base["besselK"];
  
  // Call the optim function from R in C++ 
  Rcpp::List out = besselK(Rcpp::_["x"]    = x,
                        Rcpp::_["nu"]   = nu);
  
  // Return estimated values
  return out(0);
}

// [[Rcpp::export]]
arma::vec getExpectedTau2(arma::vec beta, double sigma2, double lambda2) {
  //rgig_cpp(double chi, double psi, double lambda)
  int m = beta.n_elem;
  double chi;
  double psi;
  double lambda = 1.0 / 2.0;
  double sqrtchipsi;
  arma::vec out(m);
  for (int i = 0; i < m; i++) {
    chi = beta(i) * beta(i) / sigma2;
    psi = lambda2;
    sqrtchipsi = sqrt(chi * psi);
    out(i) = sqrt(chi) * modifiedBesselfunction2nd(sqrtchipsi, lambda + 1.0) / 
      (sqrt(psi) * modifiedBesselfunction2nd(sqrtchipsi, lambda));
  }
  return(out);
}

// [[Rcpp::export]]
arma::mat getT(arma::vec V, int p) {
  int n = V.n_elem;
  arma::mat T(n, p);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      if (i - j - 1 >= 0) {
        T(i, j) = V(i - j - 1);
      }
    }
  }
  return(T);
}


//// [[Rcpp::export]]
double obj_fun_rcpp(arma::vec beta_hat, 
                    arma::vec y, arma::mat x){
  
  arma::vec resi = y - x * beta_hat;
  double obj = arma::accu(arma::pow(resi, 2));
  return obj;
}

// [[Rcpp::export]]
arma::vec loglikelihood(arma::vec resi, double sigma2){
  
  int n = resi.n_elem;
  double pi = 3.1415926;
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    out(i) = log(1.0 / sqrt(2.0 * pi * sigma2) * 
      exp(- 1.0 / 2.0 / sigma2 * resi(i) * resi(i)));
  }
  
  return out;
}





// [[Rcpp::export]]
arma::vec optim_rcpp(arma::vec init_beta_hat,
                     arma::vec y, arma::mat x){
  
  // Extract R's optim function
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  // Call the optim function from R in C++ 
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_beta_hat,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_rcpp),
                                 Rcpp::_["method"] = "BFGS",
                                 // Pass in the other parameters as everything
                                 // is scoped environmentally
                                 Rcpp::_["y"] = y,
                                 Rcpp::_["x"] = x);
  
  // Extract out the estimated parameter values
  arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);
  
  // Return estimated values
  return out;
}


Rcpp::List fastLm(const arma::vec & y, const arma::mat & X) {
  
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y); 
  arma::colvec resid = y - X*coef; 
  
  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  
  return Rcpp::List::create(Rcpp::Named("coefs") = coef,
                            Rcpp::Named("sig2") = sig2);
}


// [[Rcpp::export]]
arma::vec idetect_rcpp(arma::vec y, int K){
  
  // Extract R's optim function
  Rcpp::Environment breakfast("package:breakfast"); 
  Rcpp::Function idetect = breakfast["sol.idetect"];
  
  // Call the optim function from R in C++ 
  Rcpp::List results = idetect(Rcpp::_["x"]    = y,
                                 // Make sure this function is not exported!
                                 Rcpp::_["thr_ic"] = 0.9,
                                 Rcpp::_["points"] = 3);
  
  // Extract out the estimated parameter values
  arma::mat results1 = Rcpp::as<arma::mat>(results["cands"]);
  //Rcpp::Rcout << results1 << std::endl;
  int ncol = results1.n_cols;
  //Rcpp::Rcout << ncol << std::endl;
  arma::vec out = results1.col(ncol - 2);
  
  // Return estimated values
  return (out.subvec(0, K - 1) - 1);
}

// [[Rcpp::export]]
arma::vec checkDim(Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> T=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> U=R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> W=R_NilValue) {
  
  int q;
  int p;
  int K;
  int r;
  
  arma::mat X_;
  arma::mat T_;
  arma::mat U_;
  arma::mat W_;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    q = X_.n_cols;
  } else {
    q = 0;
  }
  
  if (T.isNotNull()) {
    T_ = Rcpp::as<arma::mat>(T);
    p = T_.n_cols;
  } else {
    p = 0;
  }
  
  if (U.isNotNull()) {
    U_ = Rcpp::as<arma::mat>(U);
    K = U_.n_cols;
  } else {
    K = 0;
  }
  
  if (W.isNotNull()) {
    W_ = Rcpp::as<arma::mat>(W);
    r = W_.n_cols;
  } else {
    r = 0;
  }
  
  arma::vec out;
  out.zeros(4);
  out(0) = q;
  out(1) = p;
  out(2) = K;
  out(3) = r;
  
  return(out);
  
}


// [[Rcpp::export]]
Rcpp::List getGaussianPosterior(arma::vec V, arma::vec beta0, arma::vec tau2beta0,
                                double sigma2, double lambda2, 
                                Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericMatrix> T=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericMatrix> U=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericMatrix> W=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> beta1=R_NilValue, 
                                Rcpp::Nullable<Rcpp::NumericVector> beta2=R_NilValue, 
                                Rcpp::Nullable<Rcpp::NumericVector> delta0=R_NilValue, 
                                Rcpp::Nullable<Rcpp::NumericVector> delta1=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> tau2beta1=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> tau2beta2=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> tau2delta0=R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> tau2delta1=R_NilValue) {
  
  // initialize all vectors;
  
  int n = V.n_elem;
  int m = 1;
  arma::vec dim = checkDim(X, T, U, W);
  int q = dim(0);
  int p = dim(1);
  int K = dim(2);
  int r = dim(3);
  
  arma::mat X_;
  arma::mat T_;
  arma::mat U_;
  arma::mat W_;
  arma::vec beta0_ = beta0;
  arma::vec beta1_;
  arma::vec beta2_;
  arma::vec delta0_;
  arma::vec delta1_;
  arma::vec tau2beta0_ = tau2beta0;
  arma::vec tau2beta1_;
  arma::vec tau2beta2_;
  arma::vec tau2delta0_;
  arma::vec tau2delta1_;
  
  if (q > 0) {
    X_ = Rcpp::as<arma::mat>(X);
    beta1_ = Rcpp::as<arma::vec>(beta1);
    tau2beta1_ = Rcpp::as<arma::vec>(tau2beta1);
  }
  
  if (p > 0) {
    T_ = Rcpp::as<arma::mat>(T);
    beta2_ = Rcpp::as<arma::vec>(beta2);
    tau2beta2_ = Rcpp::as<arma::vec>(tau2beta2);
  } 
  
  if (K > 0) {
    U_ = Rcpp::as<arma::mat>(U);
    delta0_ = Rcpp::as<arma::vec>(delta0);
    tau2delta0_ = Rcpp::as<arma::vec>(tau2delta0);
  } 
  
  if (r > 0) {
    W_ = Rcpp::as<arma::mat>(W);
    delta1_ = Rcpp::as<arma::vec>(delta1);
    tau2delta1_ = Rcpp::as<arma::vec>(tau2delta1);
  }
  
  arma::mat UW_;
  if (K > 0) {
    if (r > 0) {
      UW_ = getUW(U_, W_);
    }
  }
  
  m = m + q + p + K + K * r;
  
  // initialize all vectors;
  
  arma::vec zeta;
  arma::vec Onebeta0;
  Onebeta0.zeros(n);
  arma::vec Xbeta1;
  Xbeta1.zeros(n);
  arma::vec Tbeta2;
  Tbeta2.zeros(n);
  arma::vec Udelta0;
  Udelta0.zeros(n);
  arma::vec UWdelta1;
  UWdelta1.zeros(n);
  arma::vec on = arma::ones(n);
  
  arma::vec betadelta;
  arma::vec tau2all;
  
  
  // update beta0;
  
  if (q > 0) {
    Xbeta1 = X_ * beta1_;
  }
  
  if (p > 0) {
    Tbeta2 = T_ * beta2_;
  }
  
  if (K > 0) {
    Udelta0 = U_ * delta0_;
    if (r > 0) {
      UWdelta1 = UW_ * delta1_;
    }
  }
  
  zeta = V - (Xbeta1 + Tbeta2 + Udelta0 + UWdelta1);
  beta0_ = getBetaNonMonotonicity(zeta, on, tau2beta0_, sigma2);
  Onebeta0 = on * beta0_;
  betadelta = beta0_;
  
  
  // update beta1;
  
  if (q > 0) {
    zeta = V - (Onebeta0 + Tbeta2 + Udelta0 + UWdelta1);
    beta1_ = getBetaNonMonotonicity(zeta, X_, tau2beta1_, sigma2);
    Xbeta1 = X_ * beta1_;
    betadelta = arma::join_cols(betadelta, beta1_);
  }
  
  
  // update beta2;
  
  
  
  if (p > 0) {
    zeta = V - (Onebeta0 + Xbeta1 + Udelta0 + UWdelta1);
    beta2_ = getBetaMonotonicity(zeta, T_, tau2beta2_, sigma2);
    Tbeta2 = T_ * beta2_;
    betadelta = arma::join_cols(betadelta, beta2_);
  }
  
  // update delta0;
  
  if (K > 0) {
    zeta = V - (Onebeta0 + Xbeta1 + Tbeta2 + UWdelta1);
    delta0_ = getBetaNonMonotonicity(zeta, U_, tau2delta0_, sigma2);
    Udelta0 = U_ * delta0_;
    betadelta = arma::join_cols(betadelta, delta0_);
    
    // update delta1;
    if (r > 0) {
      zeta = V - (Onebeta0 + Xbeta1 + Tbeta2 + Udelta0);
      delta1_ = getBetaNonMonotonicity(zeta, UW_, tau2delta1_, sigma2);
      UWdelta1 = UW_ * delta1_;
      betadelta = arma::join_cols(betadelta, delta1_);
    }
  }
  
  
  // update tau2beta0;
  
  tau2beta0_ = getTau2(beta0_, sigma2, lambda2);
  tau2all = tau2beta0_;
  
  // update tau2beta1;
  
  if (q > 0) {
    tau2beta1_ = getTau2(beta1_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta1_);
  }
  
  // update tau2beta2;
  
  if (p > 0) {
    tau2beta2_ = getTau2(beta2_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta2_);
  }
  
  // update tau2delta0;
  
  if (K > 0) {
    tau2delta0_ = getTau2(delta0_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2delta0_);
    
    // update tau2delta1;
    if(r > 0) {
      tau2delta1_ = getTau2(delta1_, sigma2, lambda2);
      tau2all = arma::join_cols(tau2all, tau2delta1_);
    }
  }
  
  
  // output;
  arma::vec fit0 = Onebeta0 + Xbeta1 + Tbeta2;
  arma::vec fit1 = Onebeta0 + Xbeta1 + Tbeta2 + Udelta0 + UWdelta1;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadelta,
    Rcpp::_["tau2all"] = tau2all,
    Rcpp::_["expectedtau2all"] = getExpectedTau2(betadelta, sigma2, lambda2),
    Rcpp::_["fit0"] = fit0,
    Rcpp::_["fit1"] = fit1,
    Rcpp::_["m"] = m,
    Rcpp::_["delta0"] = delta0_,
    Rcpp::_["delta1"] = delta1_
  );
  return(out);
}

// [[Rcpp::export]]
Rcpp::List readbetadelta(arma::vec betadelta, 
                   int q, int p, int K, int r) {
  
  int indx = 0;
  arma::vec beta0 = betadelta.subvec(indx, indx);
  indx = indx + 1;
  
  arma::vec beta1;
  if (q > 0) {
    beta1 = betadelta.subvec(indx, indx + q - 1);
    indx = indx + q;
  }
  
  arma::vec beta2;
  if (p > 0) {
    beta2 = betadelta.subvec(indx, indx + p - 1);
    indx = indx + p;
  }
  
  arma::vec delta0;
  arma::vec delta1;
  if (K > 0) {
    delta0 = betadelta.subvec(indx, indx + K - 1);
    indx = indx + K;
    if (r > 0) {
      delta1 = betadelta.subvec(indx, indx + K * r - 1);
      indx = indx + K * r;
    }
  }
  
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["beta0"] = beta0,
    Rcpp::_["beta1"] = beta1,
    Rcpp::_["beta2"] = beta2,
    Rcpp::_["delta0"] = delta0,
    Rcpp::_["delta1"] = delta1
  );
  return(out); 
}

// [[Rcpp::export]]
Rcpp::List readtau2all(arma::vec tau2all, 
                         int q, int p, int K, int r) {
  
  int indx = 0;
  arma::vec tau2beta0 = tau2all.subvec(indx, indx);
  indx = indx + 1;
  
  arma::vec tau2beta1;
  if (q > 0) {
    tau2beta1 = tau2all.subvec(indx, indx + q - 1);
    indx = indx + q;
  }
  
  arma::vec tau2beta2;
  if (p > 0) {
    tau2beta2 = tau2all.subvec(indx, indx + p - 1);
    indx = indx + p;
  }
  
  arma::vec tau2delta0;
  arma::vec tau2delta1;
  if (K > 0) {
    tau2delta0 = tau2all.subvec(indx, indx + K - 1);
    indx = indx + K;
    if (r > 0) {
      tau2delta1 = tau2all.subvec(indx, indx + K * r - 1);
      indx = indx + K *r;
    }
  }
  
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["tau2beta0"] = tau2beta0,
    Rcpp::_["tau2beta1"] = tau2beta1,
    Rcpp::_["tau2beta2"] = tau2beta2,
    Rcpp::_["tau2delta0"] = tau2delta0,
    Rcpp::_["tau2delta1"] = tau2delta1
  );
  return(out); 
}

// [[Rcpp::export]]
Rcpp::List initializeGaussianPosterior(
    arma::vec V, double lambda2, 
    Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> T=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> U=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> W=R_NilValue){
  
  // initialize all vectors;
  
  int n = V.n_elem;
  double sigma2;
  int m = 1;
  arma::vec dim = checkDim(X, T, U, W);
  int q = dim(0);
  int p = dim(1);
  int K = dim(2);
  int r = dim(3);
  
  arma::mat X_;
  arma::mat T_;
  arma::mat U_;
  arma::mat W_;
  arma::vec beta0_;
  arma::vec beta1_;
  arma::vec beta2_;
  arma::vec delta0_;
  arma::vec delta1_;
  arma::vec tau2beta0_;
  arma::vec tau2beta1_;
  arma::vec tau2beta2_;
  arma::vec tau2delta0_;
  arma::vec tau2delta1_;
  
  arma::mat XTUUW_;
  XTUUW_.ones(n, 1);
  
  arma::mat XT_;
  XT_.ones(n, 1);
  
  if (q > 0) {
    X_ = Rcpp::as<arma::mat>(X);
    XTUUW_ = arma::join_rows(XTUUW_, X_);
    XT_ = arma::join_rows(XT_, X_);
  }
  
  if (p > 0) {
    T_ = Rcpp::as<arma::mat>(T);
    XTUUW_ = arma::join_rows(XTUUW_, T_);
    XT_ = arma::join_rows(XT_, T_);
  } 
  
  if (K > 0) {
    U_ = Rcpp::as<arma::mat>(U);
    XTUUW_ = arma::join_rows(XTUUW_, U_);
  } 
  
  if (r > 0) {
    W_ = Rcpp::as<arma::mat>(W);
  }
  
  arma::mat UW_;
  if (K > 0) {
    if (r > 0) {
      UW_ = getUW(U_, W_);
      XTUUW_ = arma::join_rows(XTUUW_, UW_);
    }
  }
  
  m = m + q + p + K + K * r;
  
  // initialize beta and delta;
  
  arma::vec betadelta;
  arma::vec tau2all;
  
  betadelta.zeros(m);
  betadelta(0) = arma::mean(V);
  betadelta = optim_rcpp(betadelta, V, XTUUW_);
  //Rcpp::Rcout << betadelta << std::endl;
  sigma2 = obj_fun_rcpp(betadelta, V, XTUUW_) / n;
  
  Rcpp::List betadeltaList = readbetadelta(betadelta, q, p, K, r);
  
  beta0_ = Rcpp::as<arma::vec>(betadeltaList["beta0"]);
  beta1_ = Rcpp::as<arma::vec>(betadeltaList["beta1"]);
  beta2_ = Rcpp::as<arma::vec>(betadeltaList["beta2"]);
  delta0_ = Rcpp::as<arma::vec>(betadeltaList["delta0"]);
  delta1_ = Rcpp::as<arma::vec>(betadeltaList["delta1"]);
  
  // initialize tau2;
  
  // update tau2beta0;
  
  tau2beta0_ = getTau2(beta0_, sigma2, lambda2);
  tau2all = tau2beta0_;
  
  // update tau2beta1;
  
  if (q > 0) {
    tau2beta1_ = getTau2(beta1_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta1_);
  }
  
  // update tau2beta2;
  
  if (p > 0) {
    tau2beta2_ = getTau2(beta2_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta2_);
  }
  
  // update tau2delta0;
  
  if (K > 0) {
    tau2delta0_ = getTau2(delta0_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2delta0_);
    
    // update tau2delta1;
    if(r > 0) {
      tau2delta1_ = getTau2(delta1_, sigma2, lambda2);
      tau2all = arma::join_cols(tau2all, tau2delta1_);
    }
  }
  
  // output;
  arma::vec fit0 = XT_ * betadelta.subvec(0, 1 + q + p - 1);
  arma::vec fit1 = XTUUW_ * betadelta;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadelta,
    Rcpp::_["tau2all"] = tau2all,
    Rcpp::_["expectedtau2all"] = getExpectedTau2(betadelta, sigma2, lambda2),
    Rcpp::_["sigma2"] = sigma2,
    Rcpp::_["fit0"] = fit0,
    Rcpp::_["fit1"] = fit1,
    Rcpp::_["m"] = m,
    Rcpp::_["delta0"] = delta0_,
    Rcpp::_["delta1"] = delta1_
  );
  return(out);
  
};

// [[Rcpp::export]]
arma::mat getU(arma::vec Uvec, int n, int K) {
  
  arma::mat U(n, K);
  U.zeros();
  
  arma::vec Uvec_ = arma::sort(Uvec);
  
  int startidx;
  int endidx;
  
  for (int i = 0; i < K; i++) {
    
    if ((0 < i) & (i < (K - 1))) {
      startidx = Uvec_(i - 1) + 1;
      endidx = Uvec_(i);
    } else if (i == 0) {
      startidx = 0;
      endidx = Uvec_(0);
    } else if (i == (K - 1)) {
      startidx = Uvec_(i - 1) + 1;
      endidx = n - 1;
    }
    U.submat(startidx, i, endidx, i).ones();
  }
  
  return(U);
  
}


// [[Rcpp::export]]
arma::rowvec getTvec(arma::vec V, int p, int i) {
  
  arma::rowvec Tvec(p);
  
    for (int j = 0; j < p; j++) {
      if (i - j - 1 >= 0) {
        Tvec(j) = V(i - j - 1);
      }
    }
  
  return(Tvec);
  
}


// [[Rcpp::export]]
double dzinfbinom(double Yi, double V1i, double V2i, double psi) {
  
  double omega = 1 / (1 + exp(-V2i));
  double mu = exp(V1i);
  double out = 0.0;
  
  if (Yi == 0.0) {
    out = out * omega;
  }
  
  out = out + (1 - omega) * (R::gammafn(Yi + psi) / R::gammafn(Yi + 1) / R::gammafn(psi)) * 
    pow(mu / (mu + psi), Yi) * pow(psi/ (mu + psi), psi);
  
  return(out);
}

// [[Rcpp::export]]
double rV1i(double Yi, double V1i, double V2i, double psi, double V1ihat, double sigma21, int burnin) {
  
  double newV1 = V1i;
  double tmpV1;
  double tmppb;
  double a;
  
  for (int i = 0; i < burnin; i++) {
    tmpV1 = R::rnorm(newV1, sqrt(sigma21));
    tmppb = R::runif(0.0, 1.0);
    a = dzinfbinom(Yi, tmpV1, V2i, psi) * exp(- 1.0 / 2.0 / sigma21 * pow(tmpV1 - V1ihat, 2)) / 
      dzinfbinom(Yi, newV1, V2i, psi) * exp(- 1.0 / 2.0 / sigma21 * pow(newV1 - V1ihat, 2));
    if (tmppb < a) {
      newV1 = tmpV1;
    } 
  }
  
  return(newV1);
}

// [[Rcpp::export]]
arma::vec getV1(arma::vec Y, arma::vec V1, arma::vec V2, double psi, arma::vec V1hat, double sigma21, int burnin) {
  
  int n = Y.n_elem;
  arma::vec newV1(n);
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (int i = 0; i < n; i++) {
    //Rcpp::Rcout << i << std::endl;
    newV1(i) = rV1i(Y(i), V1(i), V2(i), psi, V1hat(i), sigma21, burnin);
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  return(newV1);
  
}

// [[Rcpp::export]]
double rV2i(double Yi, double V1i, double V2i, double psi, double V2ihat, double sigma22, int burnin) {
  
  double newV2 = V2i;
  double tmpV2;
  double tmppb;
  double a;
  
  for (int i = 0; i < burnin; i++) {
    tmpV2 = R::rnorm(newV2, sqrt(sigma22));
    tmppb = R::runif(0.0, 1.0);
    a = dzinfbinom(Yi, V1i, tmpV2, psi) * exp(- 1.0 / 2.0 / sigma22 * pow(tmpV2 - V2ihat, 2)) / 
      dzinfbinom(Yi, V1i, newV2, psi) * exp(- 1.0 / 2.0 / sigma22 * pow(newV2 - V2ihat, 2));
    if (tmppb < a) {
      newV2 = tmpV2;
    } 
  }
  
  return(newV2);
}

// [[Rcpp::export]]
arma::vec getV2(arma::vec Y, arma::vec V1, arma::vec V2, double psi, arma::vec V2hat, double sigma22, int burnin) {
  
  int n = Y.n_elem;
  arma::vec newV2(n);
  
  for (int i = 0; i < n; i++) {
    newV2(i) = rV2i(Y(i), V1(i), V2(i), psi, V2hat(i), sigma22, burnin);
  }
  
  return(newV2);
  
}

// [[Rcpp::export]]
arma::vec getGammaParm(double mean, double var) {
  arma::vec out(2);
  out(0) = pow(mean, 2) / var;
  out(1) = var / mean;
  return(out);
}

// [[Rcpp::export]]
double sumdiflogdzinbinom(arma::vec Y, arma::vec V1, arma::vec V2, double psi, double oldpsi) {
  int n = Y.n_elem;
  double tmp = 0.0;
  
  for (int i = 0; i < n; i++) {
    tmp = tmp + log(dzinfbinom(Y(i), V1(i), V2(i), psi)) - log(dzinfbinom(Y(i), V1(i), V2(i), oldpsi));
  }
  return(tmp);
}

// [[Rcpp::export]]
double getPsi(arma::vec Y, arma::vec V1, arma::vec V2, double psi, int burnin, double d1, double d2) {
  
  double a;
  
  double newpsi = psi;
  arma::vec tmpparsnew;
  arma::vec tmpparsold;
  double tmppsi;
  double tmppb;
  
  for (int i = 0; i < burnin; i++) {
    tmpparsold = getGammaParm(newpsi, newpsi);
    tmppsi = R::rgamma(tmpparsold(0), tmpparsold(1));
    tmpparsnew = getGammaParm(tmppsi, tmppsi);
    tmppb = R::runif(0.0, 1.0);
    a = sumdiflogdzinbinom(Y, V1, V2, tmppsi, newpsi);
    a = a + log(pow(tmppsi, d1 - 1) * exp(- tmppsi / d2)) - log(pow(newpsi, d1 - 1) * exp(- newpsi / d2));
    a = exp(a) * (R::dgamma(newpsi, tmpparsnew(0), tmpparsnew(1), false) / 
      R::dgamma(tmppsi, tmpparsold(0), tmpparsold(1), false));
    if (tmppb < a) {
      newpsi = tmppsi;
    }
  }
  return(newpsi);
}