// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
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

// [[Rcpp::export]]
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



double cdfKernelEstimationGaussian(double x, arma::vec beta, double bandwidth)  {
  int n = beta.n_elem;
  double p = 0; 
  int i;
  
  for (i = 0; i < n; i++) {
    p = p + R::pnorm5(x, beta(i), bandwidth, 1, 0);
  }
  
  return(p / n);
}

arma::vec getPvalueKernelEstimationGaussian(arma::mat beta2, 
                                            Rcpp::String side, double bandwidth)  {
  
  int nbeta = beta2.n_rows;
  int p = beta2.n_cols;
  arma::vec pvalue(p);
  double tmp1 = 0;
  double tmp2 = 0;
  
  int i = 0;
  int j = 0;
  
  if (side == "two-sided") {
    for (i = 0; i < p; i++) {
      
      tmp2 = cdfKernelEstimationGaussian(0.0, beta2.col(i), bandwidth);
      tmp1 = 1.0 - tmp2;
      
      if (tmp1 > tmp2) {
        pvalue(i) = 2 * (1 - tmp1);
      } else {
        pvalue(i) = 2 * (1 - tmp2);
      }
    }
  } else if (side == "upper-sided") {
    for (i = 0; i < p; i++) {
      tmp1 = cdfKernelEstimationGaussian(0.0, beta2.col(i), bandwidth);
      pvalue(i) = 1 - tmp1;
    }
  } else if (side == "lower-sided") {
    for (i = 0; i < p; i++) {
      tmp2 = cdfKernelEstimationGaussian(0.0, beta2.col(i), bandwidth);
      //tmp1 = 1.0 - tmp2;
      pvalue(i) = tmp2;
    }
  }
  
  
  return(pvalue);
}


// [[Rcpp::export]]
arma::mat BenjaminiHochberg(double FDR, arma::mat beta2, Rcpp::String side, 
                            int KernelSmoothing, double bandwidth)  {
  
  int n = beta2.n_cols;
  arma::vec pvalue(n);
  if (KernelSmoothing == 0) {
    pvalue = getPvalue(beta2, side);
  } else {
    pvalue = getPvalueKernelEstimationGaussian(
      beta2, side, bandwidth);
  }
  
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

// [[Rcpp::export]]
arma::mat BonferroniCorrection(double FAP, arma::mat beta2, Rcpp::String side, 
                               int KernelSmoothing, double bandwidth)  {
  
  int n = beta2.n_cols;
  arma::vec pvalue(n);
  if (KernelSmoothing == 0) {
    pvalue = getPvalue(beta2, side);
  } else {
    pvalue = getPvalueKernelEstimationGaussian(
      beta2, side, bandwidth);
  }
  
  int i;
  arma::vec lim(n); 
  lim.fill(FAP / n);
  arma::vec sig(n);
  
  for (i = 0; i < n; i++) {
    if (pvalue(i) < lim(i)) {
      sig(i) = 1;
    } else {
      sig(i) = 0; 
    }
  }
  
  arma::mat out(n, 4);
  out.zeros();
  
  out.col(0) = pvalue;
  out.col(2) = lim;
  out.col(3) = sig;
  
  return(out);
  
}


