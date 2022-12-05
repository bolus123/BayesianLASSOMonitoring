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
  double Z = R::pnorm5(b, 0, 1, true, false) - R::pnorm5(a, 0, 1, true, false);
  double out;
  out = R::qnorm5(tmp * Z + R::pnorm5(a, 0, 1, true, false), 0, 1, true, false) * sd + mean;
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
arma::mat getUW(arma::mat U, arma::mat W) {
  int K = U.n_cols;
  int r = W.n_cols;
  arma::mat UW(n, K * r);
  
  k = 0;
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < r; j++) {
      UW.col(k) = U.col(i) * W.col(j);
      k++;
    }
  }
  
  return(UW);
}

// [[Rcpp::export]]
arma::vec getTheta(arma::mat W, double delta0, arma::vec delta1) {
  arma::vec out(1);
  out = W * delta1;
  out(0) = out(0) + delta0;
  return(out);
}

// [[Rcpp::export]]
arma::vec getEta(arma::vec zeta, arma::vec theta, 
                 arma::vec gamma, double sigma2) {
  arma::vec eta = arma::exp(
    - 1.0 / 2.0 / sigma2 * (theta % theta + 2 * zeta * theta)
  );
}

// [[Rcpp::export]]
Rcpp::List getU(arma::vec RU, 
                arma::vec Delta, double sigma2, arma::vec gamma) {
  
  int n = RU.n_elem;
  int K = gamma.n_elem;
  
  //arma::vec R(n);
  //arma::vec tmpR;
  
  arma::mat theta(n, K);
  arma::mat prob(n, K);
  double randProb;
  double cursor;
  arma::mat tmpprobcumsum;
  arma::mat U(n, K);
  U.zeros();
  int i;
  
  for (i = 0; i < n; i++) {
    theta.row(i) = arma::trans(arma::exp(- 1.0 / 2.0 / sigma2 * 
      (Delta % Delta - 2.0 * Delta * RU(i)) + arma::log(gamma)));
    prob.row(i) = theta.row(i) / arma::accu(theta.row(i));
  }
  
  tmpprobcumsum = arma::cumsum(prob, 1);
  
  for (i = 0; i < n; i++) {
    randProb = R::runif(0, 1);
    cursor = 0;
    for (int k = 0; k < K; k++) {
      if (randProb > tmpprobcumsum(i, k)) {
        cursor = cursor + 1;
      }
    }
    U(i, cursor) = 1;
  }
  
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["theta"] = theta,
    Rcpp::_["U"] = U
  );
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
double getSigma2(arma::vec resi,
                 arma::vec BetaDelta,
                 arma::vec Tau2,
                 double a1, double a2) {
  
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
  tmprate = (arma::trans(resi) * resi + arma::trans(BetaDelta) * InvTau2 * BetaDelta) / 2.0 + a2;
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
  
  double out = rgammaCpp(shape, rate);
  return(out);
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
  //rgig_cpp(double chi, double psi, double lambda)
  double out = rgig_cpp(beta * beta / sigma2, lambda2, 1.0 / 2.0);
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
                            Rcpp::Named("sig2")       = sig2);
}

// [[Rcpp::export]]
Rcpp::List getPosteriorLayer2(
    arma::vec V, int p, int K, double a1, double a2, double b1, double b2,
    arma::vec gamma,
    int nsim, int burnin, int maxsim, bool simlambda = true, 
    Rcpp::Nullable<double> lambda=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> X=R_NilValue) {
  
  int n = V.n_elem;
  arma::mat fit0(nsim, n);
  arma::mat fit1(nsim, n);
  arma::mat T_;
  
  if (p > 0) {
    T_ = getT(V, p);
  }
  
  arma::mat X_;
  int q;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    q = X_.n_cols;
  } else {
    q = 0;
  }
  
  int m = 1 + p + q + K;
  
  //Rcpp::Rcout << q << std::endl;
  
  ////////////////////////////////////////
  
  double lambda2;
  arma::vec lambdaout(nsim);
  double lambda_;
  
  if (lambda.isNotNull()) {
    lambda_ = Rcpp::as<double>(lambda);
    lambda2 = lambda_ * lambda_;
  } else {
    lambda2 = 1000; //temporary value. This needs to be modified using others.
  }
  
  ////////////////////////////////////////
  
  int i;
  Rcpp::List Utheta;
  arma::mat U;
  arma::mat theta;
  arma::cube Uout;
  arma::cube thetaout;
  
  if (K > 0) {
    Uout.zeros(n, K, nsim);
    thetaout.zeros(n, K, nsim);
  }
  
  int tmpId;
  
  if (K > 0) {
    U.zeros(n, K);
    for (i = 0; i < n; i++) {
      tmpId = round(R::runif(-0.5, K - 1 + 0.5));
      U(i, tmpId) = 1;
    }
  }
  
  arma::mat SU = getSUMat(U);
  
  arma::mat XTSU_;
  arma::mat tmp(n, 1);
  tmp.ones();
  
  //Rcpp::Rcout << tmp.n_elem << std::endl;
  
  int indx;
  
  int complete = 0;
  
  XTSU_ = tmp;
  
  //Rcpp::Rcout << XT_ << std::endl;
  
  if (q > 0) {
    XTSU_ = arma::join_rows(XTSU_, X_);
  }
  
  if (p > 0) {
    XTSU_ = arma::join_rows(XTSU_, T_);
  }
  
  if (K > 0) {
    XTSU_ = arma::join_rows(XTSU_, SU);
  }
  
  //Rcpp::Rcout << XT_ << std::endl;
  
  double beta00 = arma::accu(V) / n;
  arma::vec Beta00(n);
  arma::vec Beta10;
  arma::vec Beta20;
  arma::vec Delta;
  
  arma::vec beta00out(nsim);
  arma::mat Beta10out;
  if (q > 0) {
    Beta10out.zeros(nsim, q);
  }
  
  arma::mat Beta20out;
  if (p > 0) {
    Beta20out.zeros(nsim, p);
  }
  
  arma::mat Deltaout;
  if (K > 0) {
    Deltaout.zeros(nsim, K);
    Delta.zeros(K);
  }
  
  arma::vec Init_BetaDeltaVec(m);
  Init_BetaDeltaVec.zeros();
  Init_BetaDeltaVec[0] = beta00;
  arma::vec tmp_BetaDeltaVec(m);
  double sigma2;
  arma::vec sigma2out(nsim);
  
  tmp_BetaDeltaVec = optim_rcpp(Init_BetaDeltaVec, V, XTSU_);
  sigma2 = obj_fun_rcpp(tmp_BetaDeltaVec, V, XTSU_) / n;
  
  beta00 = tmp_BetaDeltaVec(0);
  indx = 1;
  if (q > 0) {
    Beta10 = tmp_BetaDeltaVec.subvec(indx, indx + q - 1);
    indx = indx + q;
  }
  
  if (p > 0) {
    Beta20 = tmp_BetaDeltaVec.subvec(indx, indx + p - 1);
    indx = indx + p;
  }
  
  if (K > 0) {
    Delta = tmp_BetaDeltaVec.subvec(indx, indx + K - 1);
  }
  
  ////////////////////////////////////////
  
  
  //if (K > 0) {
  //  UDelta = getUDelta(U, Delta, n);
  //}
  
  arma::vec Tau200(1);
  arma::vec Tau210;
  arma::vec Tau220;
  arma::vec Tau2Delta;
  arma::vec Tau2;
  
  arma::vec R00(n);
  arma::vec R10(n);
  arma::vec R20(n);
  arma::vec RDelta(n);
  arma::vec RU(n);
  
  arma::vec BetaDelta;
  arma::vec tmpbeta00(1);
  arma::vec resi(n);
  
  arma::vec SUDelta;
  arma::vec XBeta10;
  arma::vec TBeta20;
  
  int cnt = 0;
  
  
  Beta00 = getBeta00(beta00, n);
  SUDelta = SU * Delta;
  XBeta10 = X_ * Beta10;
  TBeta20 = T_ * Beta20;
  
  if (maxsim >= nsim) {
    for (int sim = 0; sim < (maxsim + burnin); sim++) {
      
      Rcpp::Rcout << "sim:" << sim << std::endl;
      Rcpp::Rcout << "cnt:" << cnt << std::endl;
      Rcpp::Rcout << "sigma2:" << sigma2 << std::endl;
      Rcpp::Rcout << "lambda2:" << lambda2 << std::endl;
      
      Tau200(0) = getTau2(beta00, sigma2, lambda2);
      Tau2 = Tau200;
      
      if (q > 0) {
        Tau210.zeros(q);
        for (i = 0; i < q; i++) {
          Tau210(i) = getTau2(Beta10(i), sigma2, lambda2);
        }
        Tau2 = arma::join_cols(Tau2, Tau210);
      }
      
      //Rcpp::Rcout << Tau2 << std::endl;
      
      if (p > 0) {
        Tau220.zeros(p);
        for (i = 0; i < p; i++) {
          Tau220(i) = getTau2(Beta20(i), sigma2, lambda2);
        }
        Tau2 = arma::join_cols(Tau2, Tau220);
      }
      
      if (K > 0) {
        Tau2Delta.zeros(K);
        for (i = 0; i < K; i++) {
          Tau2Delta(i) = getTau2(Delta(i), sigma2, lambda2);
        }
        Tau2 = arma::join_cols(Tau2, Tau2Delta);
      }
      
      //Rcpp::Rcout << Tau2 << std::endl;
      
      if (simlambda == true) {
        lambda2 = getLambda2(Tau2, b1, b2);
      }
      
      ////////////////////////////////////////
      
      R00 = V;
      
      if (K > 0) {
        R00 = R00 - SUDelta;
      }
      
      if (q > 0) {
        R00 = R00 - XBeta10;
      }
      
      if (p > 0) {
        R00 = R00 - TBeta20;
      }
      
      beta00 = getbeta00(R00, Tau200(0), sigma2);
      Beta00 = getBeta00(beta00, n);
      
      ////////////////////////////////////////
      
      if (q > 0) {
        R10 = V - Beta00;
        
        if (K > 0) {
          R10 = R10 - SUDelta;
        }
        
        if (p > 0) {
          R10 = R10 - TBeta20;
        }
        
        Beta10 = getBeta10(R10, X_, Tau210, sigma2);
      }
      
      XBeta10 = X_ * Beta10;
      
      ////////////////////////////////////////
      
      if (p > 0) {
        R20 = V - Beta00;
        
        if (K > 0) {
          R20 = R20 - SUDelta;
        }
        
        if (q > 0) {
          R20 = R20 - XBeta10;
        }
        
        Beta20 = getBeta20(R20, T_, Tau220, sigma2);
      }
      
      TBeta20 = T_ * Beta20;
      
      ////////////////////////////////////////
      
      if (K > 0) {
        RDelta = V - Beta00;
        
        if (q > 0) {
          RDelta = RDelta - XBeta10;
        }
        
        if (p > 0) {
          RDelta = RDelta - TBeta20;
        }
        
        Delta = getDelta(RDelta, SU, Tau2Delta, sigma2);
        
        SUDelta = SU * Delta;
        
        RU = RDelta;
        
        for (i = 0; i < n; i++) {
          if (i > 0) {
            RU(i) = RU(i) - SUDelta(i - 1);
          }
        }
        
        Utheta = getU(RU, Delta, sigma2, gamma);
        U = Rcpp::as<arma::mat>(Utheta["U"]);
        theta = Rcpp::as<arma::mat>(Utheta["theta"]);
        //U = getU(RU, Delta, sigma2, gamma);
        
        SU = getSUMat(U);
        
        SUDelta = SU * Delta;
      }
      
      ////////////////////////////////////////
      
      tmpbeta00(0) = beta00;
      BetaDelta = tmpbeta00;
      resi = V - Beta00;
      if (q > 0) {
        resi = resi - XBeta10;
        BetaDelta = arma::join_cols(BetaDelta, Beta10);
      }
      
      if (p > 0) {
        resi = resi - TBeta20;
        BetaDelta = arma::join_cols(BetaDelta, Beta20);
      }
      
      if (K > 0) {
        resi = resi - SUDelta;
        BetaDelta = arma::join_cols(BetaDelta, Delta);
      }
      
      sigma2 = getSigma2(resi, BetaDelta, Tau2, a1, a2);
      
      ////////////////////////////////////////
      
      if (cnt == nsim) {
        complete = 1;
        break;
      }
      
      ////////////////////////////////////////
      
      if (sim >= burnin) {
        beta00out(cnt) = beta00;
        
        if (q > 0) {
          Beta10out.row(cnt) = arma::trans(Beta10);
        }
        
        if (p > 0) {
          Beta20out.row(cnt) = arma::trans(Beta20);
        }
        
        if (K > 0) {
          Deltaout.row(cnt) = arma::trans(Delta);
          Uout.slice(cnt) = U;
          thetaout.slice(cnt) = theta;
        }
        
        sigma2out(cnt) = sigma2;
        lambdaout(cnt) = sqrt(lambda2);
        
        fit1.row(cnt) = arma::trans(Beta00);
        fit0.row(cnt) = fit1.row(cnt);
        
        if (q > 0) {
          fit1.row(cnt) = fit1.row(cnt) + arma::trans(XBeta10);
          fit0.row(cnt) = fit1.row(cnt);
        }
        
        if (p > 0) {
          fit1.row(cnt) = fit1.row(cnt) + arma::trans(TBeta20);
          fit0.row(cnt) = fit1.row(cnt);
        }
        
        if (K > 0) {
          fit1.row(cnt) = fit1.row(cnt) + arma::trans(SUDelta);
        }
        
        cnt++;
      }
      
    }
  } else {
    Rcpp::Rcout << "The maximum number of simulations must be greater than the number of requested simulations" << std::endl;
  }
  
  if (complete == 0) {
    Rcpp::Rcout << "The simulation cannot be done" << std::endl;
  }
  
  Rcpp::List out;
  
  return out = Rcpp::List::create(
    Rcpp::_["beta00"] = beta00out,
    Rcpp::_["Beta10"] = Beta10out,
    Rcpp::_["Beta20"] = Beta20out,
    Rcpp::_["Delta"] = Deltaout,
    Rcpp::_["sigma2"] = sigma2out,
    Rcpp::_["lambda"] = lambdaout,
    Rcpp::_["U"] = Uout,
    Rcpp::_["theta"] = thetaout,
    Rcpp::_["fit0"] = fit0,
    Rcpp::_["fit1"] = fit1
  );
  
}
