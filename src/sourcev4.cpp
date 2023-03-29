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
        out = R::runif(lower, upper); //// Taylor series first order
  }
  
  
  return(out);
}


// [[Rcpp::export]]
double dtrnormCpp(double x, double mean, double var, double lower, double upper) {
  double sd = sqrt(var);
  double a = (lower - mean) / sd;
  double b = (upper - mean) / sd;
  
  double z1 = R::pnorm5(b, 0, 1, true, false);
  double z2 = R::pnorm5(a, 0, 1, true, false);
  
  double Z;
  double out;
  
  if ((lower <= x) & (x <= upper)) {
    
    if (z1 != z2) {
      Z = z1 - z2;
      out = R::dnorm4(x, mean, sd, false) / Z;
    } else {
      out = R::dunif(x, lower, upper, false); //// Taylor series first order
    }
    
  } else {
    out = 0;
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
Rcpp::List getUWithoutW_MC(arma::vec zetadelta, int K, 
                     arma::vec delta0, double sigma2, 
                     arma::mat Gamma, arma::vec gamma) {
  
  //Rcpp::Rcout << "a" << std::endl;
  
  int n = zetadelta.n_elem;
  arma::vec theta;
  double zetaElem;
  arma::vec eta;
  arma::mat prob(n, K);
  arma::mat out(n, K);
  out.zeros();
  
  arma::vec tmpGamma; 
  
  arma::rowvec tmpProb;

  arma::mat tmpcumsumProb;
  
  double randProb;
  int cursor;
  
  arma::rowvec tmpout; 
  
  int i;
  int j;
  //Rcpp::Rcout << "b" << std::endl;
  
  for (i = 0; i < n; i++) {
    
    theta = delta0;
    zetaElem = zetadelta(i);
    
    if (i == 0) {
      tmpGamma = gamma;
    } else {
      tmpGamma = Gamma * arma::trans(out.row(i - 1));
    }
    
    eta = getEta(zetaElem, theta, tmpGamma, sigma2);
    tmpProb = arma::trans(eta / arma::accu(eta));
    prob.row(i) = tmpProb;
    tmpcumsumProb = arma::cumsum(prob.row(i));
    
    randProb = R::runif(0, 1);
    cursor = 0;
    for (int k = 0; k < K; k++) {
      if (randProb > tmpcumsumProb(k)) {
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
  Rcpp::List out = rgig(Rcpp::_["n"] = 1,
                        Rcpp::_["lambda"]  = lambda,
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
arma::mat IsolatedShift(int T) {
  arma::mat out(T, T);
  out.eye();
  return(out);
} 

// [[Rcpp::export]]
arma::mat SustainedShift(int T) {
  
  arma::mat out(T, T - 2);
  out.zeros();
  
  int i = 0;
  int j = 0;
  for (j = 0; j < (T - 2); j++) {
    for (i = 0; i < T; i++) {
      if (i > j) {
        out(i, j) = 1.0;
      }
    }
  }
  
  return(out);
} 

// [[Rcpp::export]]
arma::mat GradualShift(int T) {
  
  arma::mat out(T, T - 2);
  out.zeros();
  
  int i = 0;
  int j = 0;
  for (j = 0; j < (T - 2); j++) {
    for (i = 0; i < T; i++) {
      if (i > j) {
        out(i, j) = (i - j) * 1.0;
      }
    }
  }
  
  return(out);
} 

// [[Rcpp::export]]
Rcpp::List getGaussianPosteriorCM(arma::vec Y, arma::vec beta0, arma::vec tau2beta0,
                                double sigma2, double lambda2, 
                                arma::mat X,
                                arma::mat V,
                                arma::vec beta1, 
                                arma::vec beta2, 
                                arma::vec tau2beta1,
                                arma::vec tau2beta2,
                                int q, int p) {
  
  // initialize all vectors;
  
  
  
  int n = Y.n_elem;
  int m = 1;
  
  //Rcpp::Rcout << 0.1 << std::endl;
  
  arma::vec Y_ = Y; 
  arma::mat X_ = X;
  arma::mat V_ = V;
  
  arma::vec beta0_ = beta0;
  arma::vec beta1_ = beta1;
  arma::vec beta2_ = beta2;

  arma::vec tau2beta0_ = tau2beta0;
  arma::vec tau2beta1_ = tau2beta1;
  arma::vec tau2beta2_ = tau2beta2;

  //Rcpp::Rcout << 0.2 << std::endl;
  
  m = m + q + p;
  
  // initialize all vectors;
  
  arma::vec zeta;
  arma::vec Onebeta0;
  Onebeta0.zeros(n);
  arma::vec Xbeta1;
  Xbeta1.zeros(n);
  arma::vec Vbeta2;
  Vbeta2.zeros(n);

  arma::vec on = arma::ones(n);
  
  arma::vec betadelta;
  arma::vec tau2all;
  
  // update beta0;
  
  if (q > 0) {
    Xbeta1 = X_ * beta1_;
  }
  
  
  if (p > 0) {
    Vbeta2 = V_ * beta2_;
  }
  
  
  //Rcpp::Rcout << "U_" << U_ << std::endl;
  //Rcpp::Rcout << "delta0_" << delta0_ << std::endl;
  
  //cpp::Rcout << "Xbeta1" << Xbeta1 << std::endl;
  //cpp::Rcout << "Tbeta2" << Tbeta2 << std::endl;
  //cpp::Rcout << "Udelta0" << Udelta0 << std::endl;
  //cpp::Rcout << "UWdelta1" << UWdelta1 << std::endl;
  
  
  zeta = Y_ - (Xbeta1 + Vbeta2);
  beta0_ = getBetaNonMonotonicity(zeta, on, tau2beta0_, sigma2);
  Onebeta0 = on * beta0_;
  betadelta = beta0_;
  
  
  // update beta1;
  
  if (q > 0) {
    zeta = Y_ - (Onebeta0 + Vbeta2);
    beta1_ = getBetaNonMonotonicity(zeta, X_, tau2beta1_, sigma2);
    Xbeta1 = X_ * beta1_;
    betadelta = arma::join_cols(betadelta, beta1_);
  }
  
  
  // update beta2;
  
  if (p > 0) {
    zeta = Y_ - (Onebeta0 + Xbeta1);
    beta2_ = getBetaMonotonicity(zeta, V_, tau2beta2_, sigma2);
    Vbeta2 = V_ * beta2_;
    betadelta = arma::join_cols(betadelta, beta2_);
  }
  
  
  // update tau2beta0;
  
  //Rcpp::Rcout << beta0_ << std::endl;
  //Rcpp::Rcout << sigma2 << std::endl;
  //Rcpp::Rcout << lambda2 << std::endl;
  
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
  
  
  // output;
  arma::vec fit0 = Onebeta0 + Vbeta2;
  arma::vec fit1 = Onebeta0 + Xbeta1 + Vbeta2;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadelta,
    Rcpp::_["tau2all"] = tau2all,
    Rcpp::_["expectedtau2all"] = getExpectedTau2(betadelta, sigma2, lambda2),
    Rcpp::_["fit0"] = fit0,
    Rcpp::_["fit1"] = fit1,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p
  );
  return(out);
}


// [[Rcpp::export]]
Rcpp::List getGaussianPosterior(arma::vec V, arma::vec beta0, arma::vec tau2beta0,
                                double sigma2, double lambda2, 
                                arma::mat X,
                                arma::mat T,
                                arma::mat U,
                                arma::mat W,
                                arma::vec beta1, 
                                arma::vec beta2, 
                                arma::vec delta0, 
                                arma::vec delta1,
                                arma::vec tau2beta1,
                                arma::vec tau2beta2,
                                arma::vec tau2delta0,
                                arma::vec tau2delta1,
                                int q, int p, int K, int r) {
  
  // initialize all vectors;
  
  
  
  int n = V.n_elem;
  int m = 1;
  
  //Rcpp::Rcout << 0.1 << std::endl;
  
  arma::mat X_ = X;
  arma::mat T_ = T;
  arma::mat U_ = U;
  arma::mat W_ = W;
  arma::vec beta0_ = beta0;
  arma::vec beta1_ = beta1;
  arma::vec beta2_ = beta2;
  arma::vec delta0_ = delta0;
  arma::vec delta1_ = delta1;
  arma::vec tau2beta0_ = tau2beta0;
  arma::vec tau2beta1_ = tau2beta1;
  arma::vec tau2beta2_ = tau2beta2;
  arma::vec tau2delta0_ = tau2delta0;
  arma::vec tau2delta1_ = tau2delta1;
  
  //Rcpp::Rcout << 0.14 << std::endl;
  
  arma::mat UW_;
  if (K > 0) {
    if (r > 0) {
      UW_ = getUW(U_, W_);
    }
  }
  
  //Rcpp::Rcout << 0.2 << std::endl;
  
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
  
  
  //Rcpp::Rcout << 0.3 << std::endl;
  
  //Rcpp::Rcout << "U_" << U_ << std::endl;
  //Rcpp::Rcout << "delta0_" << delta0_ << std::endl;
  
  //cpp::Rcout << "Xbeta1" << Xbeta1 << std::endl;
  //cpp::Rcout << "Tbeta2" << Tbeta2 << std::endl;
  //cpp::Rcout << "Udelta0" << Udelta0 << std::endl;
  //cpp::Rcout << "UWdelta1" << UWdelta1 << std::endl;
  
  
  zeta = V - (Xbeta1 + Tbeta2 + Udelta0 + UWdelta1);
  beta0_ = getBetaNonMonotonicity(zeta, on, tau2beta0_, sigma2);
  Onebeta0 = on * beta0_;
  betadelta = beta0_;
  
  
  //Rcpp::Rcout << 1 << std::endl;
  
  // update beta1;
  
  if (q > 0) {
    zeta = V - (Onebeta0 + Tbeta2 + Udelta0 + UWdelta1);
    beta1_ = getBetaNonMonotonicity(zeta, X_, tau2beta1_, sigma2);
    Xbeta1 = X_ * beta1_;
    betadelta = arma::join_cols(betadelta, beta1_);
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  // update beta2;
  
  
  
  if (p > 0) {
    zeta = V - (Onebeta0 + Xbeta1 + Udelta0 + UWdelta1);
    beta2_ = getBetaMonotonicity(zeta, T_, tau2beta2_, sigma2);
    Tbeta2 = T_ * beta2_;
    betadelta = arma::join_cols(betadelta, beta2_);
  }
  
  //Rcpp::Rcout << 3 << std::endl;
  
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
  
  //Rcpp::Rcout << 4 << std::endl;
  
  // update tau2beta0;
  
  //Rcpp::Rcout << beta0_ << std::endl;
  //Rcpp::Rcout << sigma2 << std::endl;
  //Rcpp::Rcout << lambda2 << std::endl;
  
  tau2beta0_ = getTau2(beta0_, sigma2, lambda2);
  tau2all = tau2beta0_;
  
  //Rcpp::Rcout << 5 << std::endl;
  
  // update tau2beta1;
  
  if (q > 0) {
    tau2beta1_ = getTau2(beta1_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta1_);
  }
  
  //Rcpp::Rcout << 6 << std::endl;
  
  // update tau2beta2;
  
  if (p > 0) {
    tau2beta2_ = getTau2(beta2_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta2_);
  }
  
  //Rcpp::Rcout << 7 << std::endl;
  
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
  
  //Rcpp::Rcout << 8 << std::endl;
  
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
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K"] = K,
    Rcpp::_["r"] = r
  );
  return(out);
}

// [[Rcpp::export]]
Rcpp::List getGaussianPosterior_Isolated_Sustained(
  arma::vec V, 
  arma::vec beta0, arma::vec tau2beta0,
  double sigma2, double lambda2, 
  arma::mat X,
  arma::mat T,
  arma::mat U1,
  arma::mat U2,
  arma::mat W,
  arma::vec beta1, 
  arma::vec beta2, 
  arma::vec delta10, 
  arma::vec delta11,
  arma::vec delta20, 
  arma::vec delta21,
  arma::vec tau2beta1,
  arma::vec tau2beta2,
  arma::vec tau2delta10,
  arma::vec tau2delta11,
  arma::vec tau2delta20,
  arma::vec tau2delta21,
  int q, int p, int K1, int K2, int r) {
  
  // initialize all vectors;
  
  
  
  int n = V.n_elem;
  int m = 1;
  
  //Rcpp::Rcout << 0.1 << std::endl;
  
  arma::mat X_ = X;
  arma::mat T_ = T;
  arma::mat U1_ = U1;
  arma::mat U2_ = U2;
  arma::mat W_ = W;
  arma::vec beta0_ = beta0;
  arma::vec beta1_ = beta1;
  arma::vec beta2_ = beta2;
  arma::vec delta10_ = delta10;
  arma::vec delta11_ = delta11;
  arma::vec delta20_ = delta20;
  arma::vec delta21_ = delta21;
  arma::vec tau2beta0_ = tau2beta0;
  arma::vec tau2beta1_ = tau2beta1;
  arma::vec tau2beta2_ = tau2beta2;
  arma::vec tau2delta10_ = tau2delta10;
  arma::vec tau2delta11_ = tau2delta11;
  arma::vec tau2delta20_ = tau2delta20;
  arma::vec tau2delta21_ = tau2delta21;
  //Rcpp::Rcout << 0.14 << std::endl;
  
  arma::mat UW1_;
  if (K1 > 0) {
    if (r > 0) {
      UW1_ = getUW(U1_, W_);
    }
  }
  
  arma::mat UW2_;
  if (K2 > 0) {
    if (r > 0) {
      UW2_ = getUW(U2_, W_);
    }
  }
  
  //Rcpp::Rcout << 0.2 << std::endl;
  
  m = m + q + p + K1 + K1 * r+ K2 + K2 * r;
  
  // initialize all vectors;
  
  arma::vec zeta;
  arma::vec Onebeta0;
  Onebeta0.zeros(n);
  arma::vec Xbeta1;
  Xbeta1.zeros(n);
  arma::vec Tbeta2;
  Tbeta2.zeros(n);
  arma::vec Udelta10;
  Udelta10.zeros(n);
  arma::vec UWdelta11;
  UWdelta11.zeros(n);
  arma::vec Udelta20;
  Udelta20.zeros(n);
  arma::vec UWdelta21;
  UWdelta21.zeros(n);
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
  
  if (K1 > 0) {
    Udelta10 = U1_ * delta10_;
    if (r > 0) {
      UWdelta11 = UW1_ * delta11_;
    }
  }
  
  if (K2 > 0) {
    Udelta20 = U2_ * delta20_;
    if (r > 0) {
      UWdelta21 = UW2_ * delta21_;
    }
  }
  
  //Rcpp::Rcout << 0.3 << std::endl;
  
  //Rcpp::Rcout << "U_" << U_ << std::endl;
  //Rcpp::Rcout << "delta0_" << delta0_ << std::endl;
  
  //cpp::Rcout << "Xbeta1" << Xbeta1 << std::endl;
  //cpp::Rcout << "Tbeta2" << Tbeta2 << std::endl;
  //cpp::Rcout << "Udelta0" << Udelta0 << std::endl;
  //cpp::Rcout << "UWdelta1" << UWdelta1 << std::endl;
  
  
  zeta = V - (Xbeta1 + Tbeta2 + Udelta10 + UWdelta11 + Udelta20 + UWdelta21);
  beta0_ = getBetaNonMonotonicity(zeta, on, tau2beta0_, sigma2);
  Onebeta0 = on * beta0_;
  betadelta = beta0_;
  
  
  //Rcpp::Rcout << 1 << std::endl;
  
  // update beta1;
  
  if (q > 0) {
    zeta = V - (Onebeta0 + Tbeta2 + Udelta10 + UWdelta11 + Udelta20 + UWdelta21);
    beta1_ = getBetaNonMonotonicity(zeta, X_, tau2beta1_, sigma2);
    Xbeta1 = X_ * beta1_;
    betadelta = arma::join_cols(betadelta, beta1_);
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  // update beta2;
  
  
  
  if (p > 0) {
    zeta = V - (Onebeta0 + Xbeta1 + Udelta10 + UWdelta11 + Udelta20 + UWdelta21);
    beta2_ = getBetaMonotonicity(zeta, T_, tau2beta2_, sigma2);
    Tbeta2 = T_ * beta2_;
    betadelta = arma::join_cols(betadelta, beta2_);
  }
  
  //Rcpp::Rcout << 3 << std::endl;
  
  // update delta10;
  
  if (K1 > 0) {
    zeta = V - (Onebeta0 + Xbeta1 + Tbeta2 + UWdelta11+ Udelta20 + UWdelta21);
    delta10_ = getBetaNonMonotonicity(zeta, U1_, tau2delta10_, sigma2);
    Udelta10 = U1_ * delta10_;
    betadelta = arma::join_cols(betadelta, delta10_);
    
    // update delta11;
    if (r > 0) {
      zeta = V - (Onebeta0 + Xbeta1 + Tbeta2 + Udelta10+ Udelta20 + UWdelta21);
      delta11_ = getBetaNonMonotonicity(zeta, UW1_, tau2delta11_, sigma2);
      UWdelta11 = UW1_ * delta11_;
      betadelta = arma::join_cols(betadelta, delta11_);
    }
  }
  
  // update delta0;
  
  if (K2 > 0) {
    zeta = V - (Onebeta0 + Xbeta1 + Tbeta2 + Udelta10 + UWdelta11 + UWdelta21);
    delta20_ = getBetaNonMonotonicity(zeta, U2_, tau2delta20_, sigma2);
    Udelta20 = U2_ * delta20_;
    betadelta = arma::join_cols(betadelta, delta20_);
    
    // update delta1;
    if (r > 0) {
      zeta = V - (Onebeta0 + Xbeta1 + Tbeta2 + Udelta10 + UWdelta11 + Udelta20);
      delta21_ = getBetaNonMonotonicity(zeta, UW2_, tau2delta21_, sigma2);
      UWdelta21 = UW2_ * delta21_;
      betadelta = arma::join_cols(betadelta, delta21_);
    }
  }
  
  //Rcpp::Rcout << 4 << std::endl;
  
  // update tau2beta0;
  
  //Rcpp::Rcout << beta0_ << std::endl;
  //Rcpp::Rcout << sigma2 << std::endl;
  //Rcpp::Rcout << lambda2 << std::endl;
  
  tau2beta0_ = getTau2(beta0_, sigma2, lambda2);
  tau2all = tau2beta0_;
  
  //Rcpp::Rcout << 5 << std::endl;
  
  // update tau2beta1;
  
  if (q > 0) {
    tau2beta1_ = getTau2(beta1_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta1_);
  }
  
  //Rcpp::Rcout << 6 << std::endl;
  
  // update tau2beta2;
  
  if (p > 0) {
    tau2beta2_ = getTau2(beta2_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2beta2_);
  }
  
  //Rcpp::Rcout << 7 << std::endl;
  
  // update tau2delta10;
  
  if (K1 > 0) {
    tau2delta10_ = getTau2(delta10_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2delta10_);
    
    // update tau2delta11;
    if(r > 0) {
      tau2delta11_ = getTau2(delta11_, sigma2, lambda2);
      tau2all = arma::join_cols(tau2all, tau2delta11_);
    }
  }
  
  // update tau2delta20;
  
  if (K2 > 0) {
    tau2delta20_ = getTau2(delta20_, sigma2, lambda2);
    tau2all = arma::join_cols(tau2all, tau2delta20_);
    
    // update tau2delta21;
    if(r > 0) {
      tau2delta21_ = getTau2(delta21_, sigma2, lambda2);
      tau2all = arma::join_cols(tau2all, tau2delta21_);
    }
  }
  
  //Rcpp::Rcout << 8 << std::endl;
  
  // output;
  arma::vec fit0 = Onebeta0 + Xbeta1 + Tbeta2;
  arma::vec fit1 = Onebeta0 + Xbeta1 + Tbeta2 + Udelta10 + UWdelta11 + Udelta20 + UWdelta21;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadelta,
    Rcpp::_["tau2all"] = tau2all,
    Rcpp::_["expectedtau2all"] = getExpectedTau2(betadelta, sigma2, lambda2),
    Rcpp::_["fit0"] = fit0,
    Rcpp::_["fit1"] = fit1,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K1"] = K1,
    Rcpp::_["K2"] = K2,
    Rcpp::_["r"] = r
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
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K"] = K,
    Rcpp::_["r"] = r
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
    out = omega;
  }
  
  out = out + (1.0 - omega) * (R::gammafn(Yi + psi) / R::gammafn(Yi + 1.0) / R::gammafn(psi)) * 
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
      (dzinfbinom(Yi, newV1, V2i, psi) * exp(- 1.0 / 2.0 / sigma21 * pow(newV1 - V1ihat, 2)));
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
      (dzinfbinom(Yi, V1i, newV2, psi) * exp(- 1.0 / 2.0 / sigma22 * pow(newV2 - V2ihat, 2)));
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
double proddiflogdzinbinom(arma::vec Y, arma::vec V1, arma::vec V2, double psi, double oldpsi) {
  int n = Y.n_elem;
  double tmp = 0.0;
  
  for (int i = 0; i < n; i++) {
    tmp = tmp + log(dzinfbinom(Y(i), V1(i), V2(i), psi)) - log(dzinfbinom(Y(i), V1(i), V2(i), oldpsi));
    //tmp = tmp * (dzinfbinom(Y(i), V1(i), V2(i), psi) / dzinfbinom(Y(i), V1(i), V2(i), oldpsi)) ;
  }
  return(tmp);
}

// [[Rcpp::export]]
double getPsi(arma::vec Y, arma::vec V1, arma::vec V2, double psi, int burnin, double c1, double c2) {
  
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
    //a = proddiflogdzinbinom(Y, V1, V2, tmppsi, newpsi);
    a = exp(proddiflogdzinbinom(Y, V1, V2, tmppsi, newpsi));
    //Rcpp::Rcout << a << std::endl;
    
    //a = a + log(pow(tmppsi, d1 - 1) * exp(- tmppsi / d2)) - log(pow(newpsi, d1 - 1) * exp(- newpsi / d2));
    a = a * pow(tmppsi, c1 - 1) * exp(- tmppsi / c2) / pow(newpsi, c1 - 1) / exp(- newpsi / c2);
    //a = exp(a) * (R::dgamma(newpsi, tmpparsnew(0), tmpparsnew(1), false) / 
    //  R::dgamma(tmppsi, tmpparsold(0), tmpparsold(1), false));
    a = a * (R::dgamma(newpsi, tmpparsnew(0), tmpparsnew(1), false) / 
      R::dgamma(tmppsi, tmpparsold(0), tmpparsold(1), false));
    if (tmppb < a) {
      newpsi = tmppsi;
    }
  }
  return(newpsi);
}

// [[Rcpp::export]]
double getOmega(arma::vec Y, arma::vec V2, double psi, double omega, int burnin) {
  
  double a;
  
  double newomega = omega;
  arma::vec tmpparsnew;
  arma::vec tmpparsold;
  double tmpomega;
  double tmppb;
  
  int n = V2.n_elem;
  arma::vec V1(n);
  
  V1.fill(log(omega / (1 - omega)));
  
  for (int i = 0; i < burnin; i++) {
    tmpomega = R::runif(0.0, 1.0);
    tmppb = R::runif(0.0, 1.0);
    
    a = exp(proddiflogdzinbinom(Y, V1, V2, psi, psi));
  
    if (tmppb < a) {
      newomega = tmpomega;
    }
  }
  return(newomega);
}

// [[Rcpp::export]]
arma::mat getCubeMean(arma::cube x) {
  arma::mat out = arma::mean(x, 2);
  return(out);
}

// [[Rcpp::export]]
arma::mat getUMaxProb(arma::cube prob) {
  arma::mat probmean = getCubeMean(prob);
  int n = probmean.n_rows;
  int K = probmean.n_cols;
  arma::ucolvec maxprob = arma::index_max(probmean, 1);
  arma::mat out(n, K);
  out.zeros();
  
  for (int i = 0; i < n; i++) {
    out(i, maxprob(i)) = 1;
  } 
  
  return(out);
}

// [[Rcpp::export]]
Rcpp::List getPosteriorLayer2Shift(arma::vec V, arma::vec beta0, arma::vec tau2beta0, 
                              double sigma2, double lambda2, 
                              bool updatelambda2 = true, int burnin = 50, int nsim = 100,
                              Rcpp::Nullable<Rcpp::NumericVector> gamma=R_NilValue,
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
 
  int n = V.n_elem;
 
  arma::mat betadeltaout;
  arma::mat tau2allout;
  arma::mat expectedtau2allout;
  arma::vec sigma2out;
  arma::vec lambda2out;
  arma::mat fit0out;
  arma::mat fit1out;
    
  arma::vec tmpbetadelta;
  arma::vec tmptau2all;
  arma::vec tmpexpectedtau2all;
  arma::vec tmpsigma2;
  arma::vec tmpfit0;
  arma::vec tmpfit1;
  arma::vec tmpdelta0;
  arma::vec tmpdelta1;
  
  arma::vec dim = checkDim(X, T, U, W);
  
  Rcpp::List tmpbetadeltalist;
  Rcpp::List tmptau2alllist;
  
  int q = dim(0);
  int p = dim(1);
  int K = dim(2);
  int r = dim(3);
  
  int m = 1 + q + p + K + K * r;
  
  arma::mat X_;
  arma::mat T_;
  arma::mat U_;
  arma::mat W_;
  arma::mat prob_;
  
  arma::vec beta0_ = beta0;
  arma::vec tau2beta0_ = tau2beta0;
  
  arma::vec beta1_;
  arma::vec tau2beta1_;
  arma::vec beta2_;
  arma::vec tau2beta2_;
  arma::vec delta0_;
  arma::vec tau2delta0_;
  arma::vec delta1_;
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
    
    if (r > 0) {
      W_ = Rcpp::as<arma::mat>(W);
      delta1_ = Rcpp::as<arma::vec>(delta1);
      tau2delta1_ = Rcpp::as<arma::vec>(tau2delta1);
    }
  }
  
  arma::vec sigma2_(1);
  sigma2_(0) = sigma2;
  double itersigma2 = sigma2;
  
  arma::vec lambda2_(1);
  lambda2_(0) = lambda2;
  double iterlambda2 = lambda2;
  
  ////////////////////////
  
  arma::cube Uout;
  arma::cube probout;
  arma::vec zetadelta;
  
  arma::vec gamma_;
  
  if (K > 0) {
    Uout.zeros(n, K, nsim);
    probout.zeros(n, K, nsim);
    if (gamma.isNotNull()) {
      gamma_ = Rcpp::as<arma::vec>(gamma);
    } else {
      gamma_.zeros(K);
      gamma_.fill(1.0/K);
    }
  }
  
  
  arma::vec fit0_;
  arma::vec fit1_;
  
  arma::vec resi;
  
  Rcpp::List tmplist;
  
  Rcpp::List tmpUlist;
  
  int cnt = 0;
  
  for (int sim = 0; sim < (burnin + nsim); sim++) {
    
    tmplist = getGaussianPosterior(V, beta0_, tau2beta0_,
                                   itersigma2, iterlambda2, 
                                   X_, T_, U_, W_, 
                                   beta1_, beta2_, delta0_, delta1_, 
                                   tau2beta1_, tau2beta2_, 
                                   tau2delta0_, tau2delta1_, q, p, K, r);
    
    tmpbetadelta = Rcpp::as<arma::vec>(tmplist["betadelta"]);
    tmptau2all = Rcpp::as<arma::vec>(tmplist["tau2all"]);
      
    tmpbetadeltalist = readbetadelta(tmpbetadelta, q, p, K, r);
    tmptau2alllist = readtau2all(tmptau2all, q, p, K, r);
      
    tmpbetadeltalist = readbetadelta(tmpbetadelta, q, p, K, r);
      
    beta0_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta0"]);
    tau2beta0_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta0"]);
    
    if (q > 0) {
      beta1_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta1"]);
      tau2beta1_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta1"]);
    }
    
    if (p > 0) {
      beta2_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta2"]);
      tau2beta2_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta2"]);
    }
    
    if (K > 0) {
      delta0_ = Rcpp::as<arma::vec>(tmpbetadeltalist["delta0"]);
      tau2delta0_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2delta0"]);
      if (r > 0) {
        delta1_ = Rcpp::as<arma::vec>(tmpbetadeltalist["delta1"]);
        tau2delta1_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2delta1"]);
      }
    }
    
    fit0_ = Rcpp::as<arma::vec>(tmplist["fit0"]);
    fit1_ = Rcpp::as<arma::vec>(tmplist["fit1"]);
    
    resi = V - fit1_;
    itersigma2 = arma::accu(arma::pow(resi, 2)) / n;
      
    tmpexpectedtau2all = Rcpp::as<arma::vec>(tmplist["expectedtau2all"]);
      
    zetadelta = V - fit0_;
      
    if (K > 0) {
      if (r == 0) {
        tmpUlist = getUWithoutW(zetadelta, K, delta0_, itersigma2, gamma_);
      } else {
        tmpUlist = getUWithW(zetadelta, W_, K, 
                  delta0_, delta1_, 
                  itersigma2, gamma_);
      }
      U_ = Rcpp::as<arma::mat>(tmpUlist["U"]);
      prob_ = Rcpp::as<arma::mat>(tmpUlist["prob"]);
    }
      
    if (updatelambda2 == true) {
      iterlambda2 = getLambda2EM(tmpexpectedtau2all);
    }
      
    if (sim >= burnin) {
      
      sigma2_(0) = itersigma2;
      lambda2_(0) = iterlambda2;
      
      if (cnt > 0) {
        betadeltaout = arma::join_cols(betadeltaout, arma::trans(tmpbetadelta));
        tau2allout = arma::join_cols(tau2allout, arma::trans(tmptau2all));
        expectedtau2allout = arma::join_cols(expectedtau2allout, arma::trans(tmpexpectedtau2all));
        sigma2out = arma::join_cols(sigma2out, sigma2_);
        lambda2out = arma::join_cols(lambda2out, lambda2_);
        fit0out = arma::join_cols(fit0out, arma::trans(fit0_));
        fit1out = arma::join_cols(fit1out, arma::trans(fit1_));
      } else {
        betadeltaout = arma::trans(tmpbetadelta);
        tau2allout = arma::trans(tmptau2all);
        expectedtau2allout = arma::trans(tmpexpectedtau2all);
        sigma2out = sigma2_;
        lambda2out = lambda2_;
        fit0out = arma::trans(fit0_);
        fit1out = arma::trans(fit1_);
      }
      
      Uout.slice(cnt) = U_;
      probout.slice(cnt) = prob_;
      cnt++;
    }
    
  }
 
 
 // output;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadeltaout,
    Rcpp::_["tau2all"] = tau2allout,
    Rcpp::_["sigma2"] = sigma2out,
    Rcpp::_["lambda2"] = lambda2out,
    Rcpp::_["fit0"] = fit0out,
    Rcpp::_["fit1"] = fit1out,
    Rcpp::_["U"] = Uout,
    Rcpp::_["prob"] = probout,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K"] = K,
    Rcpp::_["r"] = r
  );
  return(out);
  
}




// [[Rcpp::export]]
Rcpp::List getPosteriorLayer2NoShift(arma::vec V, arma::vec beta0, arma::vec tau2beta0, 
                                   double sigma2, double lambda2, 
                                   bool updatelambda2 = true, int burnin = 50, int nsim = 100,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> T=R_NilValue,
                                   Rcpp::Nullable<Rcpp::NumericVector> beta1=R_NilValue, 
                                   Rcpp::Nullable<Rcpp::NumericVector> beta2=R_NilValue, 
                                   Rcpp::Nullable<Rcpp::NumericVector> tau2beta1=R_NilValue,
                                   Rcpp::Nullable<Rcpp::NumericVector> tau2beta2=R_NilValue) {
  
  int n = V.n_elem;
  
  arma::mat betadeltaout;
  arma::mat tau2allout;
  arma::mat expectedtau2allout;
  arma::vec sigma2out;
  arma::vec lambda2out;
  arma::mat fit0out;
  arma::mat fit1out;
  
  arma::vec tmpbetadelta;
  arma::vec tmptau2all;
  arma::vec tmpexpectedtau2all;
  arma::vec tmpsigma2;
  arma::vec tmpfit0;
  arma::vec tmpfit1;
  
  arma::mat X_;
  arma::mat T_;
  
  //Rcpp::Rcout << 0.1 << std::endl;
  
  int q;
  int p;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    q = X_.n_cols;
  } else {
    q = 0;
  }
  
  //Rcpp::Rcout << q << std::endl;
  
  //Rcpp::Rcout << 0.2 << std::endl;
  
  if (T.isNotNull()) {
    T_ = Rcpp::as<arma::mat>(T);
    p = T_.n_cols;
  } else {
    p = 0;
  }
  
  //Rcpp::Rcout << p << std::endl;
  
  //Rcpp::Rcout << 0.3 << std::endl;
  
  Rcpp::List tmpbetadeltalist;
  Rcpp::List tmptau2alllist;
  
  
  int m = 1 + q + p;
  
  arma::vec beta0_ = beta0;
  arma::vec tau2beta0_ = tau2beta0;
  
  arma::vec beta1_;
  arma::vec tau2beta1_;
  arma::vec beta2_;
  arma::vec tau2beta2_;
  
  //Rcpp::Rcout << 0.4 << std::endl;
  
  if (q > 0) {
    beta1_ = Rcpp::as<arma::vec>(beta1);
    tau2beta1_ = Rcpp::as<arma::vec>(tau2beta1);
    //Rcpp::Rcout << 0.41 << std::endl;
  }
  
  if (p > 0) {
    beta2_ = Rcpp::as<arma::vec>(beta2);
    tau2beta2_ = Rcpp::as<arma::vec>(tau2beta2);
    //Rcpp::Rcout << 0.42 << std::endl;
  }
  
  //Rcpp::Rcout << 0.5 << std::endl;
  
  arma::vec sigma2_(1);
  sigma2_(0) = sigma2;
  double itersigma2 = sigma2;
  
  arma::vec lambda2_(1);
  lambda2_(0) = lambda2;
  double iterlambda2 = lambda2;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  ////////////////////////
  
  arma::mat U_;
  arma::mat W_;
  arma::vec delta0_;
  arma::vec delta1_;
  arma::vec tau2delta0_;
  arma::vec tau2delta1_;
  
  arma::vec fit0_;
  arma::vec fit1_;
  
  arma::vec resi;
  
  Rcpp::List tmplist;
  
  Rcpp::List tmpUlist;
  
  int cnt = 0;
  
  //Rcpp::Rcout << 2 << std::endl;
  
  for (int sim = 0; sim < (burnin + nsim); sim++) {
    
    tmplist = getGaussianPosterior(V, beta0_, tau2beta0_,
                                   itersigma2, iterlambda2, 
                                   X_, T_, U_, W_, 
                                   beta1_, beta2_, delta0_, delta1_, 
                                   tau2beta1_, tau2beta2_, 
                                   tau2delta0_, tau2delta1_, q, p, 0, 0);
    
    //Rcpp::Rcout << 3 << std::endl;
    
    tmpbetadelta = Rcpp::as<arma::vec>(tmplist["betadelta"]);
    tmptau2all = Rcpp::as<arma::vec>(tmplist["tau2all"]);
    
    tmpbetadeltalist = readbetadelta(tmpbetadelta, q, p, 0, 0);
    tmptau2alllist = readtau2all(tmptau2all, q, p, 0, 0);
    
    beta0_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta0"]);
    tau2beta0_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta0"]);
    
    if (q > 0) {
      beta1_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta1"]);
      tau2beta1_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta1"]);
    }
    
    if (p > 0) {
      beta2_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta2"]);
      tau2beta2_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta2"]);
    }
    
    //Rcpp::Rcout << 4 << std::endl;
    
    fit0_ = Rcpp::as<arma::vec>(tmplist["fit0"]);
    fit1_ = Rcpp::as<arma::vec>(tmplist["fit1"]);
    
    resi = V - fit1_;
    itersigma2 = arma::accu(arma::pow(resi, 2)) / n;
    
    tmpexpectedtau2all = Rcpp::as<arma::vec>(tmplist["expectedtau2all"]);
    
    
    if (updatelambda2 == true) {
      iterlambda2 = getLambda2EM(tmpexpectedtau2all);
    }
    
    //Rcpp::Rcout << 5 << std::endl;
    
    if (sim >= burnin) {
      
      sigma2_(0) = itersigma2;
      lambda2_(0) = iterlambda2;
      
      if (cnt > 0) {
        betadeltaout = arma::join_cols(betadeltaout, arma::trans(tmpbetadelta));
        tau2allout = arma::join_cols(tau2allout, arma::trans(tmptau2all));
        expectedtau2allout = arma::join_cols(expectedtau2allout, arma::trans(tmpexpectedtau2all));
        sigma2out = arma::join_cols(sigma2out, sigma2_);
        lambda2out = arma::join_cols(lambda2out, lambda2_);
        fit0out = arma::join_cols(fit0out, arma::trans(fit0_));
        fit1out = arma::join_cols(fit1out, arma::trans(fit1_));
      } else {
        betadeltaout = arma::trans(tmpbetadelta);
        tau2allout = arma::trans(tmptau2all);
        expectedtau2allout = arma::trans(tmpexpectedtau2all);
        sigma2out = sigma2_;
        lambda2out = lambda2_;
        fit0out = arma::trans(fit0_);
        fit1out = arma::trans(fit1_);
      }
      
      cnt++;
    }
    
    //Rcpp::Rcout << 6 << std::endl;
    
  }
  
  
  // output;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadeltaout,
    Rcpp::_["tau2all"] = tau2allout,
    Rcpp::_["sigma2"] = sigma2out,
    Rcpp::_["lambda2"] = lambda2out,
    Rcpp::_["fit0"] = fit0out,
    Rcpp::_["fit1"] = fit1out,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K"] = 0,
    Rcpp::_["r"] = 0
  );
  return(out);
  
}




// [[Rcpp::export]]
Rcpp::List getPosteriorLayer1NBShift(arma::vec Y, arma::vec V, arma::vec Vhat, 
                                   arma::vec beta0, arma::vec tau2beta0, 
                                   double sigma2, double psi, double lambda2, int p,
                                   bool updatelambda2 = true, bool updatepsi = true,
                                   int c1 = 1, int c2 = 1,
                                   int burnin1 = 50, int burnin2 = 50, int nsim = 100,
                                   Rcpp::Nullable<Rcpp::NumericVector> gamma=R_NilValue,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
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
  
  int n = Y.n_elem;
  
  arma::mat betadeltaout;
  arma::mat tau2allout;
  arma::mat expectedtau2allout;
  arma::vec sigma2out;
  arma::vec psiout;
  arma::vec lambda2out;
  arma::mat fit0out;
  arma::mat fit1out;
  arma::mat Vout;
  
  
  arma::vec tmpbetadelta;
  arma::vec tmptau2all;
  arma::vec tmpexpectedtau2all;
  arma::vec tmpsigma2;
  arma::vec tmpfit0;
  arma::vec tmpfit1;
  arma::vec tmpdelta0;
  arma::vec tmpdelta1;
  
  arma::vec beta1_;
  arma::vec tau2beta1_;
  arma::vec beta2_;
  arma::vec tau2beta2_;
  arma::vec delta0_;
  arma::vec tau2delta0_;
  arma::vec delta1_;
  arma::vec tau2delta1_;
  
  arma::mat X_;
  arma::mat T_;
  arma::mat U_;
  arma::mat W_;
  arma::mat prob_;
  
  int q = 0;
  int K = 0;
  int r = 0;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    beta1_ = Rcpp::as<arma::vec>(beta1);
    tau2beta1_ = Rcpp::as<arma::vec>(tau2beta1);
    q = X_.n_cols;
    //Rcpp::Rcout << 0.41 << std::endl;
  }
  
  if (p > 0) {
    beta2_ = Rcpp::as<arma::vec>(beta2);
    tau2beta2_ = Rcpp::as<arma::vec>(tau2beta2);
  }
  
  if (U.isNotNull()) {
    U_ = Rcpp::as<arma::mat>(U);
    delta0_ = Rcpp::as<arma::vec>(delta0);
    tau2delta0_ = Rcpp::as<arma::vec>(tau2delta0);
    K = U_.n_cols;
    if (W.isNotNull()) {
      W_ = Rcpp::as<arma::mat>(W);
      delta1_ = Rcpp::as<arma::vec>(delta1);
      tau2delta1_ = Rcpp::as<arma::vec>(tau2delta1);
      r = W_.n_cols;
    } 
  }
  
  
  Rcpp::List tmpbetadeltalist;
  Rcpp::List tmptau2alllist;
  
  int m = 1 + q + p + K + K * r;
  
  
  
  arma::vec beta0_ = beta0;
  arma::vec tau2beta0_ = tau2beta0;
  
  
  arma::vec sigma2_(1);
  sigma2_(0) = sigma2;
  double itersigma2 = sigma2;
  
  arma::vec psi_(1);
  psi_(0) = psi;
  double iterpsi = psi;
  
  arma::vec lambda2_(1);
  lambda2_(0) = lambda2;
  double iterlambda2 = lambda2;

  ////////////////////////
  
  arma::cube Uout;
  arma::cube probout;
  
  
  arma::vec zetadelta;
  
  arma::vec gamma_;
  
  if (K > 0) {
    Uout.zeros(n, K, nsim);
    probout.zeros(n, K, nsim);
    if (gamma.isNotNull()) {
      gamma_ = Rcpp::as<arma::vec>(gamma);
    } else {
      gamma_.zeros(K);
      gamma_.fill(1.0/K);
    }
  }
  
  
  arma::vec fit0_;
  arma::vec fit1_ = Vhat;
  
  arma::vec resi;
  
  Rcpp::List tmplist;
  
  Rcpp::List tmpUlist;
  
  arma::vec V1 = V;
  arma::vec V2(n);
  V2.fill(-arma::datum::inf);
  
  int cnt = 0;
  
  for (int sim1 = 0; sim1 < (burnin1 + nsim); sim1++) {
    
    V1 = getV1(Y, V1, V2, iterpsi, fit1_, itersigma2, burnin1);
    T_ = getT(V1, p);
    
    
    for (int sim2 = 0; sim2 < (burnin2 + 1); sim2++) {
      
      tmplist = getGaussianPosterior(V1, beta0_, tau2beta0_,
                                     itersigma2, iterlambda2, 
                                     X_, T_, U_, W_, 
                                     beta1_, beta2_, delta0_, delta1_, 
                                     tau2beta1_, tau2beta2_, 
                                     tau2delta0_, tau2delta1_, q, p, K, r);
      
      tmpbetadelta = Rcpp::as<arma::vec>(tmplist["betadelta"]);
      tmptau2all = Rcpp::as<arma::vec>(tmplist["tau2all"]);
      
      tmpbetadeltalist = readbetadelta(tmpbetadelta, q, p, K, r);
      tmptau2alllist = readtau2all(tmptau2all, q, p, K, r);
      
      beta0_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta0"]);
      tau2beta0_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta0"]);
      
      if (q > 0) {
        beta1_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta1"]);
        tau2beta1_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta1"]);
      }
      
      if (p > 0) {
        beta2_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta2"]);
        tau2beta2_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta2"]);
      }
      
      if (K > 0) {
        delta0_ = Rcpp::as<arma::vec>(tmpbetadeltalist["delta0"]);
        tau2delta0_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2delta0"]);
        if (r > 0) {
          delta1_ = Rcpp::as<arma::vec>(tmpbetadeltalist["delta1"]);
          tau2delta1_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2delta1"]);
        }
      }
      
      fit0_ = Rcpp::as<arma::vec>(tmplist["fit0"]);
      fit1_ = Rcpp::as<arma::vec>(tmplist["fit1"]);
      
      resi = V1 - fit1_;
      itersigma2 = arma::accu(arma::pow(resi, 2)) / n;
      
      tmpexpectedtau2all = Rcpp::as<arma::vec>(tmplist["expectedtau2all"]);
      
      zetadelta = V1 - fit0_;
      
      if (K > 0) {
        if (r == 0) {
          tmpUlist = getUWithoutW(zetadelta, K, delta0_, itersigma2, gamma_);
        } else {
          tmpUlist = getUWithW(zetadelta, W_, K, 
                               delta0_, delta1_, 
                               itersigma2, gamma_);
        }
        U_ = Rcpp::as<arma::mat>(tmpUlist["U"]);
        prob_ = Rcpp::as<arma::mat>(tmpUlist["prob"]);
      }
      
    }
    
    
    if (updatepsi == true) {
      iterpsi = getPsi(Y, V1, V2, iterpsi, burnin1, c1, c2);
    }
    
    
    if (updatelambda2 == true) {
      iterlambda2 = getLambda2EM(tmpexpectedtau2all);
    }
    
    if (sim1 >= burnin1) {
      
      sigma2_(0) = itersigma2;
      lambda2_(0) = iterlambda2;
      psi_(0) = iterpsi;
      
      if (cnt > 0) {
        betadeltaout = arma::join_cols(betadeltaout, arma::trans(tmpbetadelta));
        tau2allout = arma::join_cols(tau2allout, arma::trans(tmptau2all));
        expectedtau2allout = arma::join_cols(expectedtau2allout, arma::trans(tmpexpectedtau2all));
        sigma2out = arma::join_cols(sigma2out, sigma2_);
        psiout = arma::join_cols(psiout, psi_);
        lambda2out = arma::join_cols(lambda2out, lambda2_);
        fit0out = arma::join_cols(fit0out, arma::trans(fit0_));
        fit1out = arma::join_cols(fit1out, arma::trans(fit1_));
        Vout = arma::join_cols(Vout, arma::trans(V1));
      } else {
        betadeltaout = arma::trans(tmpbetadelta);
        tau2allout = arma::trans(tmptau2all);
        expectedtau2allout = arma::trans(tmpexpectedtau2all);
        sigma2out = sigma2_;
        psiout = psi_;
        lambda2out = lambda2_;
        fit0out = arma::trans(fit0_);
        fit1out = arma::trans(fit1_);
        Vout = arma::trans(V1);
      }
      
      Uout.slice(cnt) = U_;
      probout.slice(cnt) = prob_;
      cnt++;
    }
    
  }
  
  
  // output;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadeltaout,
    Rcpp::_["tau2all"] = tau2allout,
    Rcpp::_["psi"] = psiout,
    Rcpp::_["sigma2"] = sigma2out,
    Rcpp::_["lambda2"] = lambda2out,
    Rcpp::_["fit0"] = fit0out,
    Rcpp::_["fit1"] = fit1out,
    Rcpp::_["V"] = Vout,
    Rcpp::_["U"] = Uout,
    Rcpp::_["prob"] = probout,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K"] = K,
    Rcpp::_["r"] = r
  );
  return(out);
  
}





// [[Rcpp::export]]
Rcpp::List getPosteriorLayer1NBNoShift(arma::vec Y, arma::vec V, arma::vec Vhat, 
                                     arma::vec beta0, arma::vec tau2beta0, 
                                     double sigma2, double psi, double lambda2, int p,
                                     bool updatelambda2 = true, bool updatepsi = true,
                                     int c1 = 1, int c2 = 1,
                                     int burnin1 = 50, int burnin2 = 50, int nsim = 100,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> beta1=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> beta2=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2beta1=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2beta2=R_NilValue) {
  
  int n = Y.n_elem;
  
  arma::mat betadeltaout;
  arma::mat tau2allout;
  arma::mat expectedtau2allout;
  arma::vec sigma2out;
  arma::vec psiout;
  arma::vec lambda2out;
  arma::mat fit0out;
  arma::mat fit1out;
  arma::mat Vout;
  
  
  arma::vec tmpbetadelta;
  arma::vec tmptau2all;
  arma::vec tmpexpectedtau2all;
  arma::vec tmpsigma2;
  arma::vec tmpfit0;
  arma::vec tmpfit1;
  
  arma::vec beta1_;
  arma::vec tau2beta1_;
  arma::vec beta2_;
  arma::vec tau2beta2_;
  
  arma::mat X_;
  arma::mat T_;
  arma::mat prob_;
  
  int q = 0;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    beta1_ = Rcpp::as<arma::vec>(beta1);
    tau2beta1_ = Rcpp::as<arma::vec>(tau2beta1);
    q = X_.n_cols;
    //Rcpp::Rcout << 0.41 << std::endl;
  }
  
  if (p > 0) {
    beta2_ = Rcpp::as<arma::vec>(beta2);
    tau2beta2_ = Rcpp::as<arma::vec>(tau2beta2);
  }
  
  
  Rcpp::List tmpbetadeltalist;
  Rcpp::List tmptau2alllist;
  
  int m = 1 + q + p;
  
  arma::vec beta0_ = beta0;
  arma::vec tau2beta0_ = tau2beta0;
  
  
  arma::vec sigma2_(1);
  sigma2_(0) = sigma2;
  double itersigma2 = sigma2;
  
  arma::vec psi_(1);
  psi_(0) = psi;
  double iterpsi = psi;
  
  arma::vec lambda2_(1);
  lambda2_(0) = lambda2;
  double iterlambda2 = lambda2;
  
  ////////////////////////
  
  arma::mat U_;
  arma::mat W_;
  arma::vec delta0_;
  arma::vec delta1_;
  arma::vec tau2delta0_;
  arma::vec tau2delta1_;
  
  arma::vec fit0_;
  arma::vec fit1_ = Vhat;
  
  arma::vec resi;
  
  Rcpp::List tmplist;
  
  Rcpp::List tmpUlist;
  
  arma::vec V1 = V;
  arma::vec V2(n);
  V2.fill(-arma::datum::inf);
  
  int cnt = 0;
  
  for (int sim1 = 0; sim1 < (burnin1 + nsim); sim1++) {
    
    V1 = getV1(Y, V1, V2, iterpsi, fit1_, itersigma2, burnin2);
    T_ = getT(V1, p);
    
    //Rcpp::Rcout << "T_" << T_ << std::endl;
    
    for (int sim2 = 0; sim2 < (burnin2 + 1); sim2++) {
      
      //Rcpp::Rcout << "beta0_" << beta0_ << std::endl;
      //Rcpp::Rcout << "tau2beta0_" << tau2beta0_ << std::endl;
      //Rcpp::Rcout << "itersigma2" << itersigma2 << std::endl;
      //Rcpp::Rcout << "iterlambda2" << iterlambda2 << std::endl;
      //Rcpp::Rcout << "beta1_" << beta1_ << std::endl;
      //Rcpp::Rcout << "beta2_" << beta2_ << std::endl;
      //Rcpp::Rcout << "delta0_" << delta0_ << std::endl;
      //Rcpp::Rcout << "delta1_" << delta1_ << std::endl;
      //Rcpp::Rcout << "tau2beta1_" << tau2beta1_ << std::endl;
      //Rcpp::Rcout << "tau2beta2_" << tau2beta2_ << std::endl;
      //Rcpp::Rcout << "tau2delta0_" << tau2delta0_ << std::endl;
      //Rcpp::Rcout << "tau2delta1_" << tau2delta1_ << std::endl;
      
      tmplist = getGaussianPosterior(V1, beta0_, tau2beta0_,
                                     itersigma2, iterlambda2, 
                                     X_, T_, U_, W_, 
                                     beta1_, beta2_, delta0_, delta1_, 
                                     tau2beta1_, tau2beta2_, 
                                     tau2delta0_, tau2delta1_, q, p, 0, 0);
      
      
      
      tmpbetadelta = Rcpp::as<arma::vec>(tmplist["betadelta"]);
      tmptau2all = Rcpp::as<arma::vec>(tmplist["tau2all"]);
      
      //Rcpp::Rcout << "tmpbetadelta" << tmplist << std::endl;
      //Rcpp::Rcout << "tmptau2all" << tmplist << std::endl;
      
      tmpbetadeltalist = readbetadelta(tmpbetadelta, q, p, 0, 0);
      tmptau2alllist = readtau2all(tmptau2all, q, p, 0, 0);
      
      beta0_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta0"]);
      tau2beta0_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta0"]);
      
      if (q > 0) {
        beta1_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta1"]);
        tau2beta1_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta1"]);
      }
      
      if (p > 0) {
        beta2_ = Rcpp::as<arma::vec>(tmpbetadeltalist["beta2"]);
        tau2beta2_ = Rcpp::as<arma::vec>(tmptau2alllist["tau2beta2"]);
      }
      
      
      fit0_ = Rcpp::as<arma::vec>(tmplist["fit0"]);
      fit1_ = Rcpp::as<arma::vec>(tmplist["fit1"]);
      
      resi = V1 - fit1_;
      itersigma2 = arma::accu(arma::pow(resi, 2)) / n;
      
      tmpexpectedtau2all = Rcpp::as<arma::vec>(tmplist["expectedtau2all"]);
      
    }
    
    
    if (updatepsi == true) {
      iterpsi = getPsi(Y, V1, V2, iterpsi, burnin2, c1, c2);
    }
    
    
    if (updatelambda2 == true) {
      iterlambda2 = getLambda2EM(tmpexpectedtau2all);
    }
    
    if (sim1 >= burnin1) {
      
      sigma2_(0) = itersigma2;
      lambda2_(0) = iterlambda2;
      psi_(0) = iterpsi;
      
      if (cnt > 0) {
        betadeltaout = arma::join_cols(betadeltaout, arma::trans(tmpbetadelta));
        tau2allout = arma::join_cols(tau2allout, arma::trans(tmptau2all));
        expectedtau2allout = arma::join_cols(expectedtau2allout, arma::trans(tmpexpectedtau2all));
        sigma2out = arma::join_cols(sigma2out, sigma2_);
        psiout = arma::join_cols(psiout, psi_);
        lambda2out = arma::join_cols(lambda2out, lambda2_);
        fit0out = arma::join_cols(fit0out, arma::trans(fit0_));
        fit1out = arma::join_cols(fit1out, arma::trans(fit1_));
        Vout = arma::join_cols(Vout, arma::trans(V1));
      } else {
        betadeltaout = arma::trans(tmpbetadelta);
        tau2allout = arma::trans(tmptau2all);
        expectedtau2allout = arma::trans(tmpexpectedtau2all);
        sigma2out = sigma2_;
        psiout = psi_;
        lambda2out = lambda2_;
        fit0out = arma::trans(fit0_);
        fit1out = arma::trans(fit1_);
        Vout = arma::trans(V1);
      }
      
      cnt++;
    }
    
  }
  
  
  // output;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta"] = betadeltaout,
    Rcpp::_["tau2all"] = tau2allout,
    Rcpp::_["psi"] = psiout,
    Rcpp::_["sigma2"] = sigma2out,
    Rcpp::_["lambda2"] = lambda2out,
    Rcpp::_["fit0"] = fit0out,
    Rcpp::_["fit1"] = fit1out,
    Rcpp::_["V"] = Vout,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K"] = 0,
    Rcpp::_["r"] = 0
  );
  return(out);
  
}





// [[Rcpp::export]]
Rcpp::List getPosteriorLayer1ZinfNBNoShift(arma::vec Y, 
                                       arma::vec V1, arma::vec V1hat, 
                                       arma::vec V2, arma::vec V2hat, 
                                       arma::vec beta01, arma::vec tau2beta01, 
                                       arma::vec beta02, arma::vec tau2beta02, 
                                       double sigma21, double sigma22, double psi, double lambda2, int p,
                                       bool updatelambda2 = true, bool updatepsi = true,
                                       int c1 = 1, int c2 = 1,
                                       int burnin1 = 50, int burnin2 = 50, int nsim = 100,
                                       Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericVector> beta11=R_NilValue, 
                                       Rcpp::Nullable<Rcpp::NumericVector> beta21=R_NilValue, 
                                       Rcpp::Nullable<Rcpp::NumericVector> tau2beta11=R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericVector> tau2beta21=R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericVector> beta12=R_NilValue, 
                                       Rcpp::Nullable<Rcpp::NumericVector> beta22=R_NilValue, 
                                       Rcpp::Nullable<Rcpp::NumericVector> tau2beta12=R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericVector> tau2beta22=R_NilValue) {
  
  int n = Y.n_elem;
  
  arma::mat betadelta1out;
  arma::mat tau2all1out;
  arma::mat expectedtau2all1out;
  arma::vec sigma21out;
  arma::mat fit01out;
  arma::mat fit11out;
  arma::mat V1out;
  
  arma::mat betadelta2out;
  arma::mat tau2all2out;
  arma::mat expectedtau2all2out;
  arma::vec sigma22out;
  arma::mat fit02out;
  arma::mat fit12out;
  arma::mat V2out;
  
  arma::vec psiout;
  arma::vec lambda2out;
  
  arma::vec tmpbetadelta1;
  arma::vec tmptau2all1;
  arma::vec tmpexpectedtau2all1;
  arma::vec tmpsigma21;
  arma::vec tmpfit01;
  arma::vec tmpfit11;
  
  arma::vec tmpbetadelta2;
  arma::vec tmptau2all2;
  arma::vec tmpexpectedtau2all2;
  arma::vec tmpsigma22;
  arma::vec tmpfit02;
  arma::vec tmpfit12;
  
  arma::vec beta11_;
  arma::vec tau2beta11_;
  arma::vec beta21_;
  arma::vec tau2beta21_;
  
  arma::vec beta12_;
  arma::vec tau2beta12_;
  arma::vec beta22_;
  arma::vec tau2beta22_;
  
  arma::mat X_;
  arma::mat T1_;
  arma::mat T2_;
  
  int q = 0;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    beta11_ = Rcpp::as<arma::vec>(beta11);
    tau2beta11_ = Rcpp::as<arma::vec>(tau2beta11);
    beta12_ = Rcpp::as<arma::vec>(beta12);
    tau2beta12_ = Rcpp::as<arma::vec>(tau2beta12);
    q = X_.n_cols;
    //Rcpp::Rcout << 0.41 << std::endl;
  }
  
  if (p > 0) {
    beta21_ = Rcpp::as<arma::vec>(beta21);
    tau2beta21_ = Rcpp::as<arma::vec>(tau2beta22);
    beta22_ = Rcpp::as<arma::vec>(beta22);
    tau2beta22_ = Rcpp::as<arma::vec>(tau2beta22);
  }
  
  
  Rcpp::List tmpbetadelta1list;
  Rcpp::List tmptau2all1list;
  Rcpp::List tmpbetadelta2list;
  Rcpp::List tmptau2all2list;
  
  int m = 1 + q + p;
  
  arma::vec beta01_ = beta01;
  arma::vec tau2beta01_ = tau2beta01;
  arma::vec beta02_ = beta02;
  arma::vec tau2beta02_ = tau2beta02;
  
  arma::vec sigma21_(1);
  sigma21_(0) = sigma21;
  double itersigma21 = sigma21;
  
  arma::vec sigma22_(1);
  sigma22_(0) = sigma22;
  double itersigma22 = sigma22;
  
  arma::vec psi_(1);
  psi_(0) = psi;
  double iterpsi = psi;
  
  arma::vec lambda2_(1);
  lambda2_(0) = lambda2;
  double iterlambda2 = lambda2;
  
  ////////////////////////
  
  arma::mat U1_;
  arma::mat U2_;
  arma::mat W_;
  arma::vec delta01_;
  arma::vec delta11_;
  arma::vec tau2delta01_;
  arma::vec tau2delta11_;
  arma::vec delta02_;
  arma::vec delta12_;
  arma::vec tau2delta02_;
  arma::vec tau2delta12_;
  
  arma::vec fit01_;
  arma::vec fit11_ = V1hat;
  
  arma::vec fit02_;
  arma::vec fit12_ = V2hat;
  
  arma::vec resi1;
  arma::vec resi2;
  
  Rcpp::List tmp1list;
  Rcpp::List tmpU1list;
  Rcpp::List tmp2list;
  Rcpp::List tmpU2list;
  
  arma::vec V1_ = V1;
  arma::vec V2_ = V2;
  
  int cnt = 0;
  
  for (int sim1 = 0; sim1 < (burnin1 + nsim); sim1++) {
    
    //// update V1 and T1
    V1_ = getV1(Y, V1_, V2_, iterpsi, fit11_, itersigma21, burnin1);
    T1_ = getT(V1_, p);
    
    for (int sim2 = 0; sim2 < (burnin2 + 1); sim2++) {
      
      //// update beta, delta and tau2
      
      tmp1list = getGaussianPosterior(V1_, beta01_, tau2beta01_,
                                     itersigma21, iterlambda2, 
                                     X_, T1_, U1_, W_, 
                                     beta11_, beta21_, delta01_, delta11_, 
                                     tau2beta11_, tau2beta21_, 
                                     tau2delta01_, tau2delta11_, q, p, 0, 0);
      
      tmpbetadelta1 = Rcpp::as<arma::vec>(tmp1list["betadelta"]);
      tmptau2all1 = Rcpp::as<arma::vec>(tmp1list["tau2all"]);
      
      tmpbetadelta1list = readbetadelta(tmpbetadelta1, q, p, 0, 0);
      tmptau2all1list = readtau2all(tmptau2all1, q, p, 0, 0);
      
      beta01_ = Rcpp::as<arma::vec>(tmpbetadelta1list["beta0"]);
      tau2beta01_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2beta0"]);
      
      if (q > 0) {
        beta11_ = Rcpp::as<arma::vec>(tmpbetadelta1list["beta1"]);
        tau2beta11_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2beta1"]);
      }
      
      if (p > 0) {
        beta21_ = Rcpp::as<arma::vec>(tmpbetadelta1list["beta2"]);
        tau2beta21_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2beta2"]);
      }
      
      //// update fit
      
      fit01_ = Rcpp::as<arma::vec>(tmp1list["fit0"]);
      fit11_ = Rcpp::as<arma::vec>(tmp1list["fit1"]);
      
      //// update sigma2
      
      resi1 = V1_ - fit11_;
      itersigma21 = arma::accu(arma::pow(resi1, 2)) / n;
      
      tmpexpectedtau2all1 = Rcpp::as<arma::vec>(tmp1list["expectedtau2all"]);
      
    }
    
    //// update V2 and T2
    V2_ = getV2(Y, V1_, V2_, iterpsi, fit12_, itersigma22, burnin1);
    T2_ = getT(V2_, p);
    
    for (int sim2 = 0; sim2 < (burnin2 + 1); sim2++) {
      
      //// update beta, delta and tau2
      
      tmp2list = getGaussianPosterior(V2_, beta02_, tau2beta02_,
                                      itersigma22, iterlambda2, 
                                      X_, T2_, U2_, W_, 
                                      beta12_, beta22_, delta02_, delta12_, 
                                      tau2beta12_, tau2beta22_, 
                                      tau2delta02_, tau2delta12_, q, p, 0, 0);
      
      tmpbetadelta2 = Rcpp::as<arma::vec>(tmp2list["betadelta"]);
      tmptau2all2 = Rcpp::as<arma::vec>(tmp2list["tau2all"]);
      
      tmpbetadelta2list = readbetadelta(tmpbetadelta2, q, p, 0, 0);
      tmptau2all2list = readtau2all(tmptau2all2, q, p, 0, 0);
      
      beta02_ = Rcpp::as<arma::vec>(tmpbetadelta2list["beta0"]);
      tau2beta02_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2beta0"]);
      
      if (q > 0) {
        beta12_ = Rcpp::as<arma::vec>(tmpbetadelta2list["beta1"]);
        tau2beta12_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2beta1"]);
      }
      
      if (p > 0) {
        beta22_ = Rcpp::as<arma::vec>(tmpbetadelta2list["beta2"]);
        tau2beta22_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2beta2"]);
      }
      
      //// update fit
      fit02_ = Rcpp::as<arma::vec>(tmp2list["fit0"]);
      fit12_ = Rcpp::as<arma::vec>(tmp2list["fit1"]);
      
      
      //// update sigma2
      resi2 = V2_ - fit12_;
      itersigma22 = arma::accu(arma::pow(resi2, 2)) / n;
      
      tmpexpectedtau2all2 = Rcpp::as<arma::vec>(tmp1list["expectedtau2all"]);
      
    }
    
    
    //// update psi
    if (updatepsi == true) {
      iterpsi = getPsi(Y, V1_, V2_, iterpsi, burnin1, c1, c2);
    }
    
    //// update lambda2
    if (updatelambda2 == true) {
      iterlambda2 = getLambda2EM(arma::join_cols(tmpexpectedtau2all1, tmpexpectedtau2all2));
    }
    
    if (sim1 >= burnin1) {
      
      sigma21_(0) = itersigma21;
      sigma22_(0) = itersigma22;
      lambda2_(0) = iterlambda2;
      psi_(0) = iterpsi;
      
      if (cnt > 0) {
        betadelta1out = arma::join_cols(betadelta1out, arma::trans(tmpbetadelta1));
        tau2all1out = arma::join_cols(tau2all1out, arma::trans(tmptau2all1));
        expectedtau2all1out = arma::join_cols(expectedtau2all1out, arma::trans(tmpexpectedtau2all1));
        sigma21out = arma::join_cols(sigma21out, sigma21_);
        fit01out = arma::join_cols(fit01out, arma::trans(fit01_));
        fit11out = arma::join_cols(fit11out, arma::trans(fit11_));
        V1out = arma::join_cols(V1out, arma::trans(V1_));
        
        betadelta2out = arma::join_cols(betadelta2out, arma::trans(tmpbetadelta2));
        tau2all2out = arma::join_cols(tau2all2out, arma::trans(tmptau2all2));
        expectedtau2all2out = arma::join_cols(expectedtau2all2out, arma::trans(tmpexpectedtau2all2));
        sigma22out = arma::join_cols(sigma22out, sigma22_);
        fit02out = arma::join_cols(fit02out, arma::trans(fit02_));
        fit12out = arma::join_cols(fit12out, arma::trans(fit12_));
        V2out = arma::join_cols(V2out, arma::trans(V2_));
        
        psiout = arma::join_cols(psiout, psi_);
        lambda2out = arma::join_cols(lambda2out, lambda2_);
      } else {
        betadelta1out = arma::trans(tmpbetadelta1);
        tau2all1out = arma::trans(tmptau2all1);
        expectedtau2all1out = arma::trans(tmpexpectedtau2all1);
        sigma21out = sigma21_;
        fit01out = arma::trans(fit01_);
        fit11out = arma::trans(fit11_);
        V1out = arma::trans(V1_);
        
        betadelta2out = arma::trans(tmpbetadelta2);
        tau2all2out = arma::trans(tmptau2all2);
        expectedtau2all2out = arma::trans(tmpexpectedtau2all2);
        sigma22out = sigma22_;
        fit02out = arma::trans(fit02_);
        fit12out = arma::trans(fit12_);
        V2out = arma::trans(V2_);
        
        psiout = psi_;
        lambda2out = lambda2_;
      }
      
      cnt++;
    }
    
  }
  
  
  // output;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["betadelta1"] = betadelta1out,
    Rcpp::_["tau2all1"] = tau2all1out,
    Rcpp::_["sigma21"] = sigma21out,
    Rcpp::_["fit01"] = fit01out,
    Rcpp::_["fit11"] = fit11out,
    Rcpp::_["V1"] = V1out,
    Rcpp::_["betadelta2"] = betadelta2out,
    Rcpp::_["tau2all2"] = tau2all2out,
    Rcpp::_["sigma22"] = sigma22out,
    Rcpp::_["fit02"] = fit02out,
    Rcpp::_["fit12"] = fit12out,
    Rcpp::_["V2"] = V2out,
    Rcpp::_["psi"] = psiout,
    Rcpp::_["lambda2"] = lambda2out,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p,
    Rcpp::_["K"] = 0,
    Rcpp::_["r"] = 0
  );
  return(out);
  
}



template <typename T>
inline void set_item_impl( List& target, int i, const T& obj, CharacterVector& names, traits::true_type ){
  target[i] = obj.object ;
  names[i] = obj.name ;
}

template <typename T>
inline void set_item_impl( List& target, int i, const T& obj, CharacterVector&, traits::false_type ){
  target[i] = obj ;
}

template <typename T>
inline void set_item( List& target, int i, const T& obj, CharacterVector& names){
  set_item_impl( target, i, obj, names, typename traits::is_named<T>::type() ) ;
}

class SequentialInserter{
public:
  SequentialInserter( List& target_, CharacterVector& names_ ) : target(target_), names(names_), index(0){}
  
  class SequentialInserterProxy{
  public:
    SequentialInserterProxy( SequentialInserter& parent_, const char* name_ ) : parent(parent_), name(name_){}
    
    template <typename T>
    void operator=( const T& obj){
      parent.set( name, obj ) ;
    }  
    
  private:
    SequentialInserter& parent ;
    const char* name ;
  } ;
  
  inline SequentialInserterProxy operator[]( const char* name){
    return SequentialInserterProxy( *this, name );    
  }
  
  template <typename T>
  void set( const char* name, const T& obj){
    target[index] = obj ;
    names[index] = name ;
    index++ ;
  }
  
private:
  List& target ;
  CharacterVector& names;
  int index ;
  
} ;




// [[Rcpp::export]]
Rcpp::List getPosteriorLayer1ZinfNBShift(arma::vec Y,
                                     arma::vec V1, arma::vec V1hat, 
                                     arma::vec V2, arma::vec V2hat, 
                                     arma::vec beta01, arma::vec tau2beta01, 
                                     arma::vec beta02, arma::vec tau2beta02, 
                                     double sigma21, double sigma22, 
                                     double psi, double lambda2, int p,
                                     bool updatelambda2 = true, bool updatepsi = true,
                                     int c1 = 1, int c2 = 1,
                                     int burnin1 = 50, int burnin2 = 50, int nsim = 100,
                                     Rcpp::Nullable<Rcpp::NumericVector> gamma1=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> gamma2=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> U1=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> U2=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> W=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> beta11=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> beta21=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2beta11=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2beta21=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> beta12=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> beta22=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2beta12=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2beta22=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> delta01=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> delta11=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> delta02=R_NilValue, 
                                     Rcpp::Nullable<Rcpp::NumericVector> delta12=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2delta01=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2delta11=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2delta02=R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> tau2delta12=R_NilValue) {
  
  int n = Y.n_elem;
  
  arma::mat betadelta1out;
  arma::mat tau2all1out;
  arma::mat expectedtau2all1out;
  arma::vec sigma21out;
  arma::mat fit01out;
  arma::mat fit11out;
  arma::mat V1out;
  
  arma::mat betadelta2out;
  arma::mat tau2all2out;
  arma::mat expectedtau2all2out;
  arma::vec sigma22out;
  arma::mat fit02out;
  arma::mat fit12out;
  arma::mat V2out;
  
  arma::vec psiout;
  arma::vec lambda2out;
  
  arma::vec tmpbetadelta1;
  arma::vec tmptau2all1;
  arma::vec tmpexpectedtau2all1;
  arma::vec tmpsigma21;
  arma::vec tmpfit01;
  arma::vec tmpfit11;
  arma::vec tmpdelta01;
  arma::vec tmpdelta11;
  
  arma::vec tmpbetadelta2;
  arma::vec tmptau2all2;
  arma::vec tmpexpectedtau2all2;
  arma::vec tmpsigma22;
  arma::vec tmpfit02;
  arma::vec tmpfit12;
  arma::vec tmpdelta02;
  arma::vec tmpdelta12;
  
  arma::vec beta11_;
  arma::vec tau2beta11_;
  arma::vec beta21_;
  arma::vec tau2beta21_;
  arma::vec delta01_;
  arma::vec tau2delta01_;
  arma::vec delta11_;
  arma::vec tau2delta11_;
  
  arma::vec beta12_;
  arma::vec tau2beta12_;
  arma::vec beta22_;
  arma::vec tau2beta22_;
  arma::vec delta02_;
  arma::vec tau2delta02_;
  arma::vec delta12_;
  arma::vec tau2delta12_;
  
  arma::mat X_;
  arma::mat T1_;
  arma::mat U1_;
  arma::mat T2_;
  arma::mat U2_;
  arma::mat W_;
  arma::mat prob1_;
  arma::mat prob2_;
  
  int q = 0;
  int K = 0;
  int r = 0;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    beta11_ = Rcpp::as<arma::vec>(beta11);
    tau2beta11_ = Rcpp::as<arma::vec>(tau2beta11);
    beta12_ = Rcpp::as<arma::vec>(beta12);
    tau2beta12_ = Rcpp::as<arma::vec>(tau2beta12);
    q = X_.n_cols;
    //Rcpp::Rcout << 0.41 << std::endl;
  }
  
  if (p > 0) {
    beta21_ = Rcpp::as<arma::vec>(beta21);
    tau2beta21_ = Rcpp::as<arma::vec>(tau2beta21);
    beta22_ = Rcpp::as<arma::vec>(beta22);
    tau2beta22_ = Rcpp::as<arma::vec>(tau2beta22);
  }
  
  int K1;
  int K2;
  
  if (U1.isNotNull()) {
    U1_ = Rcpp::as<arma::mat>(U1);
    delta01_ = Rcpp::as<arma::vec>(delta01);
    tau2delta01_ = Rcpp::as<arma::vec>(tau2delta01);
    K1 = U1_.n_cols;
  }
  
  if (U2.isNotNull()) {
    U2_ = Rcpp::as<arma::mat>(U2);
    delta02_ = Rcpp::as<arma::vec>(delta02);
    tau2delta02_ = Rcpp::as<arma::vec>(tau2delta02);
    K2 = U2_.n_cols;
  }
  
  if (K1 > K2) {
    K = K2;
  } else {
    K = K1;
  }
  
  if (K > 0) {
    if (W.isNotNull()) {
      W_ = Rcpp::as<arma::mat>(W);
      delta11_ = Rcpp::as<arma::vec>(delta11);
      tau2delta11_ = Rcpp::as<arma::vec>(tau2delta11);
      delta12_ = Rcpp::as<arma::vec>(delta12);
      tau2delta12_ = Rcpp::as<arma::vec>(tau2delta12);
      r = W_.n_cols;
    } 
  }
  
  
  Rcpp::List tmpbetadelta1list;
  Rcpp::List tmptau2all1list;
  
  Rcpp::List tmpbetadelta2list;
  Rcpp::List tmptau2all2list;
  
  int m = 1 + q + p + K + K * r;
  
  arma::vec beta01_ = beta01;
  arma::vec tau2beta01_ = tau2beta01;
  
  arma::vec beta02_ = beta02;
  arma::vec tau2beta02_ = tau2beta02;
  
  arma::vec sigma21_(1);
  sigma21_(0) = sigma21;
  double itersigma21 = sigma21;
  
  arma::vec sigma22_(1);
  sigma22_(0) = sigma22;
  double itersigma22 = sigma22;
  
  arma::vec psi_(1);
  psi_(0) = psi;
  double iterpsi = psi;
  
  arma::vec lambda2_(1);
  lambda2_(0) = lambda2;
  double iterlambda2 = lambda2;
  
  ////////////////////////
  
  arma::cube U1out;
  arma::cube prob1out;
  arma::vec zetadelta1;
  arma::vec gamma1_;
  
  arma::cube U2out;
  arma::cube prob2out;
  arma::vec zetadelta2;
  arma::vec gamma2_;
  
  if (K > 0) {
    
    U1out.zeros(n, K, nsim);
    prob1out.zeros(n, K, nsim);
    U2out.zeros(n, K, nsim);
    prob2out.zeros(n, K, nsim);
    
    if (gamma1.isNotNull()) {
      gamma1_ = Rcpp::as<arma::vec>(gamma1);
    } else {
      gamma1_.zeros(K);
      gamma1_.fill(1.0/K);
    }
    
    if (gamma2.isNotNull()) {
      gamma2_ = Rcpp::as<arma::vec>(gamma2);
    } else {
      gamma2_.zeros(K);
      gamma2_.fill(1.0/K);
    }
  }
  
  arma::vec fit01_;
  arma::vec fit11_ = V1hat;
  arma::vec resi1;
  
  arma::vec fit02_;
  arma::vec fit12_ = V2hat;
  arma::vec resi2;
  
  Rcpp::List tmp1list;
  Rcpp::List tmpU1list;
  
  Rcpp::List tmp2list;
  Rcpp::List tmpU2list;
  
  arma::vec V1_ = V1;
  arma::vec V2_ = V2;
  
  int cnt = 0;
  
  for (int sim1 = 0; sim1 < (burnin1 + nsim); sim1++) {
    
    V1_ = getV1(Y, V1_, V2_, iterpsi, fit11_, itersigma21, burnin1);
    T1_ = getT(V1_, p);
    
    
    for (int sim2 = 0; sim2 < (burnin2 + 1); sim2++) {
      
      tmp1list = getGaussianPosterior(V1_, beta01_, tau2beta01_,
                                     itersigma21, iterlambda2, 
                                     X_, T1_, U1_, W_, 
                                     beta11_, beta21_, delta01_, delta11_, 
                                     tau2beta11_, tau2beta21_, 
                                     tau2delta01_, tau2delta11_, q, p, K, r);
      
      tmpbetadelta1 = Rcpp::as<arma::vec>(tmp1list["betadelta"]);
      tmptau2all1 = Rcpp::as<arma::vec>(tmp1list["tau2all"]);
      
      tmpbetadelta1list = readbetadelta(tmpbetadelta1, q, p, K, r);
      tmptau2all1list = readtau2all(tmptau2all1, q, p, K, r);
      
      beta01_ = Rcpp::as<arma::vec>(tmpbetadelta1list["beta0"]);
      tau2beta01_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2beta0"]);
      
      if (q > 0) {
        beta11_ = Rcpp::as<arma::vec>(tmpbetadelta1list["beta1"]);
        tau2beta11_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2beta1"]);
      }
      
      if (p > 0) {
        beta21_ = Rcpp::as<arma::vec>(tmpbetadelta1list["beta2"]);
        tau2beta21_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2beta2"]);
      }
      
      if (K > 0) {
        delta01_ = Rcpp::as<arma::vec>(tmpbetadelta1list["delta0"]);
        tau2delta01_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2delta0"]);
        if (r > 0) {
          delta11_ = Rcpp::as<arma::vec>(tmpbetadelta1list["delta1"]);
          tau2delta11_ = Rcpp::as<arma::vec>(tmptau2all1list["tau2delta1"]);
        }
      }
      
      fit01_ = Rcpp::as<arma::vec>(tmp1list["fit0"]);
      fit11_ = Rcpp::as<arma::vec>(tmp1list["fit1"]);
      
      resi1 = V1_ - fit11_;
      itersigma21 = arma::accu(arma::pow(resi1, 2)) / n;
      
      tmpexpectedtau2all1 = Rcpp::as<arma::vec>(tmp1list["expectedtau2all"]);
      
      zetadelta1 = V1_ - fit01_;
      
      if (K > 0) {
        if (r == 0) {
          tmpU1list = getUWithoutW(zetadelta1, K, delta01_, itersigma21, gamma1_);
        } else {
          tmpU1list = getUWithW(zetadelta1, W_, K, 
                               delta01_, delta11_, 
                               itersigma21, gamma1_);
        }
        U1_ = Rcpp::as<arma::mat>(tmpU1list["U"]);
        prob1_ = Rcpp::as<arma::mat>(tmpU1list["prob"]);
      }
      
    }
    
    ///////////////
    
    V2_ = getV2(Y, V1_, V2_, iterpsi, fit12_, itersigma22, burnin1);
    T2_ = getT(V2_, p);
    
    for (int sim2 = 0; sim2 < (burnin2 + 1); sim2++) {
      
      tmp2list = getGaussianPosterior(V2_, beta02_, tau2beta02_,
                                     itersigma22, iterlambda2, 
                                     X_, T2_, U2_, W_, 
                                     beta12_, beta22_, delta02_, delta12_, 
                                     tau2beta12_, tau2beta22_, 
                                     tau2delta02_, tau2delta12_, q, p, K, r);
      
      tmpbetadelta2 = Rcpp::as<arma::vec>(tmp2list["betadelta"]);
      tmptau2all2 = Rcpp::as<arma::vec>(tmp2list["tau2all"]);
      
      tmpbetadelta2list = readbetadelta(tmpbetadelta2, q, p, K, r);
      tmptau2all2list = readtau2all(tmptau2all2, q, p, K, r);
      
      beta02_ = Rcpp::as<arma::vec>(tmpbetadelta2list["beta0"]);
      tau2beta02_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2beta0"]);
      
      if (q > 0) {
        beta12_ = Rcpp::as<arma::vec>(tmpbetadelta2list["beta1"]);
        tau2beta12_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2beta1"]);
      }
      
      if (p > 0) {
        beta22_ = Rcpp::as<arma::vec>(tmpbetadelta2list["beta2"]);
        tau2beta22_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2beta2"]);
      }
      
      if (K > 0) {
        delta02_ = Rcpp::as<arma::vec>(tmpbetadelta2list["delta0"]);
        tau2delta02_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2delta0"]);
        if (r > 0) {
          delta12_ = Rcpp::as<arma::vec>(tmpbetadelta2list["delta1"]);
          tau2delta12_ = Rcpp::as<arma::vec>(tmptau2all2list["tau2delta1"]);
        }
      }
      
      fit02_ = Rcpp::as<arma::vec>(tmp2list["fit0"]);
      fit12_ = Rcpp::as<arma::vec>(tmp2list["fit1"]);
      
      resi2 = V2_ - fit12_;
      itersigma22 = arma::accu(arma::pow(resi2, 2)) / n;
      
      tmpexpectedtau2all2 = Rcpp::as<arma::vec>(tmp2list["expectedtau2all"]);
      
      zetadelta2 = V2_ - fit02_;
      
      if (K > 0) {
        if (r == 0) {
          tmpU2list = getUWithoutW(zetadelta2, K, delta02_, itersigma22, gamma2_);
        } else {
          tmpU2list = getUWithW(zetadelta2, W_, K, 
                               delta02_, delta12_, 
                               itersigma22, gamma2_);
        }
        U2_ = Rcpp::as<arma::mat>(tmpU2list["U"]);
        prob2_ = Rcpp::as<arma::mat>(tmpU2list["prob"]);
      }
      
    }
    
    ///////////////
    
    if (updatepsi == true) {
      iterpsi = getPsi(Y, V1_, V2_, iterpsi, burnin1, c1, c2);
    }
    
    if (updatelambda2 == true) {
      iterlambda2 = getLambda2EM(arma::join_cols(tmpexpectedtau2all1, tmpexpectedtau2all2));
    }
    
    if (sim1 >= burnin1) {
      
      sigma21_(0) = itersigma21;
      sigma22_(0) = itersigma22;
      lambda2_(0) = iterlambda2;
      psi_(0) = iterpsi;
      
      if (cnt > 0) {
        betadelta1out = arma::join_cols(betadelta1out, arma::trans(tmpbetadelta1));
        tau2all1out = arma::join_cols(tau2all1out, arma::trans(tmptau2all1));
        expectedtau2all1out = arma::join_cols(expectedtau2all1out, arma::trans(tmpexpectedtau2all1));
        sigma21out = arma::join_cols(sigma21out, sigma21_);
        fit01out = arma::join_cols(fit01out, arma::trans(fit01_));
        fit11out = arma::join_cols(fit11out, arma::trans(fit11_));
        V1out = arma::join_cols(V1out, arma::trans(V1_));
        
        betadelta2out = arma::join_cols(betadelta2out, arma::trans(tmpbetadelta2));
        tau2all2out = arma::join_cols(tau2all2out, arma::trans(tmptau2all2));
        expectedtau2all2out = arma::join_cols(expectedtau2all2out, arma::trans(tmpexpectedtau2all2));
        sigma22out = arma::join_cols(sigma22out, sigma22_);
        fit02out = arma::join_cols(fit02out, arma::trans(fit02_));
        fit12out = arma::join_cols(fit12out, arma::trans(fit12_));
        V2out = arma::join_cols(V2out, arma::trans(V2_));
        
        psiout = arma::join_cols(psiout, psi_);
        lambda2out = arma::join_cols(lambda2out, lambda2_);
      } else {
        betadelta1out = arma::trans(tmpbetadelta1);
        tau2all1out = arma::trans(tmptau2all1);
        expectedtau2all1out = arma::trans(tmpexpectedtau2all1);
        sigma21out = sigma21_;
        fit01out = arma::trans(fit01_);
        fit11out = arma::trans(fit11_);
        V1out = arma::trans(V1_);
        
        betadelta2out = arma::trans(tmpbetadelta2);
        tau2all2out = arma::trans(tmptau2all2);
        expectedtau2all2out = arma::trans(tmpexpectedtau2all2);
        sigma22out = sigma22_;
        fit02out = arma::trans(fit02_);
        fit12out = arma::trans(fit12_);
        V2out = arma::trans(V2_);
        
        psiout = psi_;
        lambda2out = lambda2_;
      }
      
      U1out.slice(cnt) = U1_;
      prob1out.slice(cnt) = prob1_;
      
      U2out.slice(cnt) = U2_;
      prob2out.slice(cnt) = prob2_;
      
      cnt++;
    }
    
  }
  
  
  // output;
  Rcpp::List out(23);
  Rcpp::CharacterVector names(23);
  SequentialInserter inserter(out, names);
  
  inserter["betadelta1"] = betadelta1out;
  inserter["tau2all1"] = tau2all1out;
  inserter["sigma21"] = sigma21out;
  inserter["fit01"] = fit01out;
  inserter["fit11"] = fit11out;
  inserter["V1"] = V1out;
  inserter["U1"] = U1out;
  inserter["prob1"] = prob1out;
  inserter["betadelta2"] = betadelta2out;
  inserter["tau2all2"] = tau2all2out;
  inserter["sigma22"] = sigma22out;
  inserter["fit02"] = fit02out;
  inserter["fit12"] = fit12out;
  inserter["V2"] = V2out;
  inserter["U2"] = U2out;
  inserter["prob2"] = prob2out;
  inserter["psi"] = psiout;
  inserter["lambda2"] = lambda2out;
  inserter["m"] = m;
  inserter["q"] = q;
  inserter["p"] = p;
  inserter["K"] = K;
  inserter["r"] = r;
  out.names() = names;
  
  //out = Rcpp::List::create(
  //  Rcpp::_["betadelta1"] = betadelta1out,
  //  Rcpp::_["tau2all1"] = tau2all1out,
  //  Rcpp::_["sigma21"] = sigma21out,
  //  Rcpp::_["fit01"] = fit01out,
  //  Rcpp::_["fit11"] = fit11out,
  //  Rcpp::_["V1"] = V1out,
  //  Rcpp::_["U1"] = U1out,
  //  Rcpp::_["prob1"] = prob1out,
  //  Rcpp::_["betadelta2"] = betadelta2out,
  //  Rcpp::_["tau2all2"] = tau2all2out,
  //  Rcpp::_["sigma22"] = sigma22out,
  //  Rcpp::_["fit02"] = fit02out,
  //  Rcpp::_["fit12"] = fit12out,
  //  Rcpp::_["V2"] = V2out,
  //  Rcpp::_["U2"] = U2out,
  //  Rcpp::_["prob2"] = prob2out,
  //  Rcpp::_["psi"] = psiout,
  //  Rcpp::_["lambda2"] = lambda2out,
  //  Rcpp::_["Yfit0"] = Yfit0,
  //  Rcpp::_["Yfit1"] = Yfit1,
  //  Rcpp::_["m"] = m,
  //  Rcpp::_["q"] = q,
  //  Rcpp::_["p"] = p,
  //  Rcpp::_["K"] = K,
  //  Rcpp::_["r"] = r
  //);
  return(out);
  
}





//// [[Rcpp::export]]
//double quad(arma::vec beta, arma::vec tau2, double sigma2) {
//  int n = beta.n_elem;
//  arma::mat Tau2(n, n);
//  Tau2 = arma::diagmat(tau2) * sigma2;
//  arma::vec tmp = arma::trans(beta) * Tau2 * beta;
//  double out = tmp(0);
//  return(out);
//}
//
//// [[Rcpp::export]]
//double quadzinfnb(arma::vec beta1, arma::vec beta2, 
//            arma::vec tau21, arma::vec tau22, 
//            double sigma21, double sigma22) {
//  int n1 = beta1.n_elem;
//  int n2 = beta2.n_elem;
//  arma::mat Tau21(n1, n1);
//  arma::mat Tau22(n2, n2);
//  Tau21 = arma::diagmat(tau21) / sigma21;
//  Tau22 = arma::diagmat(tau22) / sigma22;
//  arma::vec tmp = arma::trans(beta1) * Tau21 * beta1  + 
//    arma::trans(beta2) * Tau22 * beta2;
//  double out = tmp(0);
//  return(out);
//}

//// [[Rcpp::export]]
//double KL(double mu1, double mu2, double sigma21, double sigma22) {
//  
//  double KL = log(sqrt(sigma22) / sqrt(sigma21)) + (sigma21 + pow(mu1 - mu2, 2)) / 2.0 / sigma22 - 1.0 / 2.0;
//  return(KL);
//  
//}

// [[Rcpp::export]]
arma::vec loglikelihoodLayer3(arma::vec V, arma::vec beta0, arma::vec tau2beta0,
                           double sigma2, double lambda2, 
                           arma::vec gamma,
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
  
  arma::vec fit(n);
  
  // update beta0;
  
  Onebeta0 = on * beta0_;
  
  fit = Onebeta0;
  
  if (q > 0) {
    Xbeta1 = X_ * beta1_;
    fit = fit + Xbeta1;
  }
  
  if (p > 0) {
    Tbeta2 = T_ * beta2_;
    fit = fit + Tbeta2;
  }
  
  if (K > 0) {
    Udelta0 = U_ * delta0_;
    fit = fit + Udelta0;
    if (r > 0) {
      UWdelta1 = UW_ * delta1_;
      fit = fit + UWdelta1;
    }
  }
  
  //
  
  arma::vec tmp(n);
  arma::rowvec Uvec;
  arma::vec llU(1);
  
  for (int i = 0; i < n; i++) {
    llU(0) = 0.0;
    if (K > 0) {
      Uvec = U_.row(i);
      llU = Uvec * log(gamma);
    }
    tmp(i) = log(R::dnorm4(V(i), fit(i), sqrt(sigma2), false)) + 
      llU(0);
  }
  
  return(tmp);
}

// [[Rcpp::export]]
double loglikelihoodLayer2(arma::vec Y, arma::vec V, double psi) {
  
  int n = Y.n_elem;
  double tmp = 0;
  double mu;
  double pb;
  
  for (int i = 0; i < n; i++) {
    mu = exp(V(i));
    pb = psi / (mu + psi);
    tmp = tmp + log(R::dnbinom(Y(i), psi, pb, false));
  }
  
  return(tmp);
}

// [[Rcpp::export]]
double loglikelihoodLayer1(arma::vec Y, arma::vec V1, arma::vec V2, double psi) {
  
  int n = Y.n_elem;
  double tmp = 0;

  for (int i = 0; i < n; i++) {
    tmp = tmp + log(dzinfbinom(Y(i), V1(i), V2(i), psi));
  }
  
  return(tmp);
}




//
//// [[Rcpp::export]]
//double loglikelihoodRatioLayer31(arma::vec V, arma::vec fit0, arma::vec fit1, 
//                                double sigma2, arma::mat U, arma::vec gamma) {
//  
//  int n = V.n_elem;
//  double tmp = 0;
//  arma::rowvec Uvec;
//
//  
//  for (int i = 0; i < n; i++) {
//
//    tmp = tmp + (-2.0) * (log(R::dnorm4(V(i), fit0(i), sqrt(sigma2), false))
//                            - log(R::dnorm4(V(i), fit1(i), sqrt(sigma2), false)));
//  }
//  
//  return(tmp);
//}