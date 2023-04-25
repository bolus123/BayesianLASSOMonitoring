// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat getInvTau2(arma::vec Tau2) {
  
  int q = Tau2.n_elem;
  arma::mat out(q, q);
  out.zeros();
  
  for (int i = 0; i < q; i++) {
    
    out(i, i) = 1 / Tau2(i);
    
  }
  
  return(out);
  
}

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

arma::vec getEta(double zetaElem, arma::vec theta, 
                 arma::vec gamma, double sigma2) {
  arma::vec eta = arma::exp(
    - 1.0 / 2.0 / sigma2 * (theta % theta - 2.0 * theta * zetaElem) + arma::log(gamma)
  );
  return(eta);
}

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

double getLambda2EM(arma::vec ExpectedTau2) {
  int m = ExpectedTau2.n_elem;
  double out = sqrt(2.0 * m / arma::accu(ExpectedTau2));
  return(out);
}

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


arma::vec getTau2(arma::vec beta, double sigma2, double lambda2) {
  //rgig_cpp(double chi, double psi, double lambda)
  int m = beta.n_elem;
  arma::vec out(m);
  for (int i = 0; i < m; i++) {
    out(i) = rgig_cpp(beta(i) * beta(i) / sigma2, lambda2, 1.0 / 2.0);
  }
  return(out);
}

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

double obj_fun_rcpp(arma::vec beta_hat, 
                    arma::vec y, arma::mat x){
  
  arma::vec resi = y - x * beta_hat;
  double obj = arma::accu(arma::pow(resi, 2));
  return obj;
}

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

Rcpp::List getGaussianPosteriorCM(arma::vec Y, arma::mat V, arma::mat X,
                                  arma::vec beta0, arma::vec beta1, arma::vec beta2, 
                                  arma::vec tau02, arma::vec tau12, arma::vec tau22,
                                  double sigma2, double lambda2, 
                                  int q, int p) {
  
  // initialize all vectors;
  
  int n = Y.n_elem;
  int m = 1;
  
  //Rcpp::Rcout << 0.1 << std::endl;
  
  arma::vec Y_ = Y; 
  arma::mat V_ = V;
  arma::mat X_ = X;
  
  arma::vec beta0_ = beta0;
  arma::vec beta1_ = beta1;
  arma::vec beta2_ = beta2;
  
  arma::vec tau02_ = tau02;
  arma::vec tau12_ = tau12;
  arma::vec tau22_ = tau22;
  
  //Rcpp::Rcout << 0.2 << std::endl;
  
  m = m + q + p;
  
  // initialize all vectors;
  
  arma::vec zeta;
  arma::vec Onebeta0;
  Onebeta0.zeros(n);
  
  arma::vec Vbeta1;
  Vbeta1.zeros(n);
  
  arma::vec Xbeta2;
  Xbeta2.zeros(n);
  
  arma::vec on = arma::ones(n);
  
  arma::vec beta;
  arma::vec tau2;
  
  // update beta0;
  
  if (q > 0) {
    Vbeta1 = V_ * beta1_;
  }
  
  if (p > 0) {
    Xbeta2 = X_ * beta2_;
  }
  
  zeta = Y_ - (Xbeta2 + Vbeta1);
  beta0_ = getBetaNonMonotonicity(zeta, on, tau02_, sigma2);
  Onebeta0 = on * beta0_;
  beta = beta0_;
  
  // update beta1;
  
  if (q > 0) {
    zeta = Y_ - (Onebeta0 + Xbeta2);
    beta1_ = getBetaMonotonicity(zeta, V_, tau12_, sigma2);
    Vbeta1 = V_ * beta1_;
    beta = arma::join_cols(beta, beta1_);
  }
  
  // update beta2;
  
  if (p > 0) {
    zeta = Y_ - (Onebeta0 + Vbeta1);
    beta2_ = getBetaNonMonotonicity(zeta, X_, tau22_, sigma2);
    Xbeta2 = X_ * beta2_;
    beta = arma::join_cols(beta, beta2_);
  }
  
  tau02_ = getTau2(beta0_, sigma2, lambda2);
  tau2 = tau02_;
  
  // update tau12;
  
  if (q > 0) {
    tau12_ = getTau2(beta1_, sigma2, lambda2);
    tau2 = arma::join_cols(tau2, tau12_);
  }
  
  // update tau22;
  
  if (p > 0) {
    tau22_ = getTau2(beta2_, sigma2, lambda2);
    tau2 = arma::join_cols(tau2, tau22_);
  }
  
  // output;
  arma::vec fit0 = Onebeta0 + Vbeta1;
  arma::vec fit1 = Onebeta0 + Xbeta2 + Vbeta1;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["beta"] = beta,
    Rcpp::_["tau2"] = tau2,
    Rcpp::_["expectedTau2"] = getExpectedTau2(beta, sigma2, lambda2),
    Rcpp::_["fit0"] = fit0,
    Rcpp::_["fit1"] = fit1,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q,
    Rcpp::_["p"] = p
  );
  return(out);
}

Rcpp::List getGaussianPosteriorCMH0(arma::vec Y, arma::mat V, 
                                  arma::vec beta0, arma::vec beta1, 
                                  arma::vec tau02, arma::vec tau12,
                                  double sigma2, double lambda2, 
                                  int q) {
  
  // initialize all vectors;
  
  int n = Y.n_elem;
  int m = 1;
  
  //Rcpp::Rcout << 0.1 << std::endl;
  
  arma::vec Y_ = Y; 
  arma::mat V_ = V;
  
  arma::vec beta0_ = beta0;
  arma::vec beta1_ = beta1;
  
  arma::vec tau02_ = tau02;
  arma::vec tau12_ = tau12;
  
  //Rcpp::Rcout << 0.2 << std::endl;
  
  m = m + q;
  
  // initialize all vectors;
  
  arma::vec zeta;
  arma::vec Onebeta0;
  Onebeta0.zeros(n);
  
  arma::vec Vbeta1;
  Vbeta1.zeros(n);
  
  arma::vec on = arma::ones(n);
  
  arma::vec beta;
  arma::vec tau2;
  
  // update beta0;
  
  if (q > 0) {
    Vbeta1 = V_ * beta1_;
  }
  
  zeta = Y_ - (Vbeta1);
  beta0_ = getBetaNonMonotonicity(zeta, on, tau02_, sigma2);
  Onebeta0 = on * beta0_;
  beta = beta0_;
  
  // update beta1;
  
  if (q > 0) {
    zeta = Y_ - (Onebeta0);
    beta1_ = getBetaMonotonicity(zeta, V_, tau12_, sigma2);
    Vbeta1 = V_ * beta1_;
    beta = arma::join_cols(beta, beta1_);
  }
  
  tau02_ = getTau2(beta0_, sigma2, lambda2);
  tau2 = tau02_;
  
  // update tau12;
  
  if (q > 0) {
    tau12_ = getTau2(beta1_, sigma2, lambda2);
    tau2 = arma::join_cols(tau2, tau12_);
  }
  
  
  // output;
  arma::vec fit = Onebeta0 + Vbeta1;
  Rcpp::List out;
  out = Rcpp::List::create(
    Rcpp::_["beta"] = beta,
    Rcpp::_["tau2"] = tau2,
    Rcpp::_["expectedTau2"] = getExpectedTau2(beta, sigma2, lambda2),
    Rcpp::_["fit"] = fit,
    Rcpp::_["m"] = m,
    Rcpp::_["q"] = q
  );
  return(out);
}


// [[Rcpp::export]]
Rcpp::List getPosterior(arma::vec Y, arma::mat V, arma::mat X, double lambda2, 
                        arma::vec beta0, arma::vec beta1, arma::vec beta2, 
                        int burnin, int nsim) {
  
  int T = Y.n_elem;
  
  int q = V.n_cols;
  int p = X.n_cols;
  
  double lambda2_ = lambda2;
  
  //initialize output
  Rcpp::List out;
  
  arma::mat betaout(nsim, 1 + p + q);
  arma::vec sigma2out(nsim); 
  arma::vec lambda2out(nsim);
  
  //initialize subjects
  Rcpp::List m1;
  arma::vec fit1;
  
  arma::vec beta(1 + p + q);
  beta(0) = beta0(0);
  beta.subvec(1, q) = beta1;
  beta.subvec((q + 1), (q + p)) = beta2;
  
  arma::vec tau2(1 + p + q); 
  arma::vec tau02(1); 
  arma::vec tau12(q);
  arma::vec tau22(p);
  
  arma::vec expectedTau2(1 + p + q); 
  
  double sigma2 = arma::var(Y - arma::ones(T) * beta0 - V * beta1 - X * beta2);
  
  //
  
  
  int i = 0;
  int k = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    tau2 = getTau2(beta, sigma2, lambda2_);
    
    tau02 = tau2(0);
    tau12 = tau2.subvec(1, q);
    tau22 = tau2.subvec(q + 1, q + p);
    
    beta0 = beta(0);
    beta1 = beta.subvec(1, q);
    beta2 = beta.subvec(q + 1, q + p);
    
    m1 = getGaussianPosteriorCM(Y, V, X,
                                beta0, beta1, beta2, 
                                tau02, tau12, tau22,
                                sigma2, lambda2_, 
                                q, p);
    
    beta = Rcpp::as<arma::vec>(Rcpp::wrap(m1["beta"]));
    tau2 = Rcpp::as<arma::vec>(Rcpp::wrap(m1["tau2"]));
    
    expectedTau2 = Rcpp::as<arma::vec>(Rcpp::wrap(m1["expectedTau2"]));
    
    fit1 = Rcpp::as<arma::vec>(Rcpp::wrap(m1["fit1"]));
    
    sigma2 = arma::var(Y - fit1);
    lambda2_ = getLambda2EM(expectedTau2);
    
    if (i >= burnin) {
      betaout.row(i - burnin) = arma::trans(beta);
      sigma2out(i - burnin) = sigma2;
      lambda2out(i - burnin) = lambda2_;
    }
    
  }
  
  out = Rcpp::List::create(
    Rcpp::_["beta"] = betaout,
    Rcpp::_["sigma2"] = sigma2out,
    Rcpp::_["lambda2"] = lambda2out
  );
  
  return(out);

}

// [[Rcpp::export]]
Rcpp::List getPosteriorH0(arma::vec Y, arma::mat V, double lambda2, 
                        arma::vec beta0, arma::vec beta1, 
                        int burnin, int nsim) {
  
  int T = Y.n_elem;
  
  int q = V.n_cols;
  
  double lambda2_ = lambda2;
  
  //initialize output
  Rcpp::List out;
  
  arma::mat betaout(nsim, 1 + q);
  arma::vec sigma2out(nsim); 
  arma::vec lambda2out(nsim);
  
  //initialize subjects
  Rcpp::List m1;
  arma::vec fit1;
  
  arma::vec beta(1 + q);
  beta(0) = beta0(0);
  beta.subvec(1, q) = beta1;
  
  arma::vec tau2(1 + q); 
  arma::vec tau02(1); 
  arma::vec tau12(q);
  
  arma::vec expectedTau2(1 + q); 
  
  double sigma2 = arma::var(Y - arma::ones(T) * beta0 - V * beta1);
  
  //
  
  
  int i = 0;
  int k = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    tau2 = getTau2(beta, sigma2, lambda2_);
    
    tau02 = tau2(0);
    tau12 = tau2.subvec(1, q);
    
    beta0 = beta(0);
    beta1 = beta.subvec(1, q);
    
    arma::mat tmpX(1, 1);
    tmpX.zeros();
    arma::vec tmpbeta2(1);
    tmpbeta2.zeros();
    
    m1 = getGaussianPosteriorCM(Y, V, tmpX,
                                beta0, beta1, tmpbeta2, 
                                tau02, tau12, tmpbeta2,
                                sigma2, lambda2_, 
                                q, 0);
    
    beta = Rcpp::as<arma::vec>(Rcpp::wrap(m1["beta"]));
    tau2 = Rcpp::as<arma::vec>(Rcpp::wrap(m1["tau2"]));
    
    expectedTau2 = Rcpp::as<arma::vec>(Rcpp::wrap(m1["expectedTau2"]));
    
    fit1 = Rcpp::as<arma::vec>(Rcpp::wrap(m1["fit1"]));
    
    sigma2 = arma::var(Y - fit1);
    lambda2_ = getLambda2EM(expectedTau2);
    
    if (i >= burnin) {
      betaout.row(i - burnin) = arma::trans(beta);
      sigma2out(i - burnin) = sigma2;
      lambda2out(i - burnin) = lambda2_;
    }
    
  }
  
  out = Rcpp::List::create(
    Rcpp::_["beta"] = betaout,
    Rcpp::_["sigma2"] = sigma2out,
    Rcpp::_["lambda2"] = lambda2out
  );
  
  return(out);
  
}


arma::vec getPvalue(arma::mat beta2, Rcpp::String side)  {
  
  int nbeta = beta2.n_rows;
  int p = beta2.n_cols;
  arma::vec pvalue(p);
  double tmp1 = 0;
  double tmp2 = 0;
  
  int i = 0;
  int j = 0;
  
  if (side == "two-sided") {
    for (i = 0; i < p; i++) {
      tmp1 = 0;
      tmp2 = 0;
      for (j = 0; j < nbeta; j++) {
        if (beta2(j, i) >= 0) {
          tmp1 = tmp1 + 1;
        }
        if (beta2(j, i) <= 0) {
          tmp2 = tmp2 + 1;
        }
      }
      tmp1 = tmp1 / nbeta;
      tmp2 = tmp2 / nbeta;
      if (tmp1 > tmp2) {
        pvalue(i) = 2 * (1 - tmp1);
      } else {
        pvalue(i) = 2 * (1 - tmp2);
      }
    }
  } else if (side == "upper-sided") {
    for (i = 0; i < p; i++) {
      tmp1 = 0;
      for (j = 0; j < nbeta; j++) {
        if (beta2(j, i) <= 0) {
          tmp1 = tmp1 + 1;
        }
      }
      tmp1 = tmp1 / nbeta;
      pvalue(i) = 1 - tmp1;
    }
  } else if (side == "lower-sided") {
    for (i = 0; i < p; i++) {
      tmp1 = 0;
      for (j = 0; j < nbeta; j++) {
        if (beta2(j, i) >= 0) {
          tmp1 = tmp1 + 1;
        }
      }
      tmp1 = tmp1 / nbeta;
      pvalue(i) = 1 - tmp1;
    }
  }
  
  
  return(pvalue);
}

// [[Rcpp::export]]
arma::mat BenjaminiHochberg(double FDR, arma::mat beta2, Rcpp::String side)  {
  
  int n = beta2.n_cols;
  arma::vec pvalue = getPvalue(beta2, side);
  arma::uvec idx = arma::sort_index(pvalue);

  int i;
  arma::vec rank(n);
  arma::vec lim(n); 
  arma::vec sig(n);
  
  for (i = 0; i < n; i++) {
    rank(idx(i)) = i + 1.0;
  }
  
  lim = rank / n * FDR;
  
  for (i = 0; i < n; i++) {
    if (pvalue(i) < lim(i)) {
      sig(i) = 1;
    } else {
      sig(i) = 0; 
    }
  }
  
  arma::mat out(n, 4);
  
  out.col(0) = pvalue;
  out.col(1) = rank;
  out.col(2) = lim;
  out.col(3) = sig;
  
  return(out);
  
}
