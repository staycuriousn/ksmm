// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
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
  
  // [[Rcpp::depends(RcppArmadillo)]] 
  // [[Rcpp::export]]
List smo_smm(ListMatrix K, List init_sd, arma::mat W, NumericVector init_a, NumericVector y, double C, int dim_x, int dim_y, double epsilon, int maxit) {
  
  int n = init_sd.length();
  
  double tb;
  int n1;
  int n2;
  double maxDiff;
  arma::mat AA(dim_x, dim_y);
  AA.fill(0.0);
  int S;
  double tempSum1;
  double tempSum2;
  
  NumericVector a = clone(init_a);
  // NumericVector a_old(a.length());
  NumericVector a_old = clone(init_a);
  double a1_new;
  double a2_new;
  double U;
  double V;
  double b1;
  double b2;
  double b = 0.0;
  List sd(n);
  
  
  for (int i = 0; i < n; ++i) {
    arma::mat sd_mat_tmp = init_sd[i];
    sd[i] = sd_mat_tmp;
  }
  
  for (int l = 0; l < maxit; ++l) {
    
    IntegerVector ns_vec(0);
    double f_norm2 = accu(W % W);
    
    for (int i = 0; i < n; ++i) {
      arma::mat sd_mat = sd[i];
      double mul_sd = accu(sd_mat % W);
      
      tb = y[i] * (mul_sd / std::sqrt(f_norm2) + b);
      
      if ((std::abs(tb - 1.0) < epsilon) & (a[i] <= C) & (a[i] >= 0.0)) {
        
      } else if ((tb > 1.0) & (std::abs(a[i] - 0.0) < epsilon)) {
        
      } else if ((tb < 1.0) & (std::abs(a[i] - C) < epsilon)) {
        
      } else {
        ns_vec.push_back(i);
      }
    }
    if (ns_vec.length() < 1) break;
    // Rprintf("%i \n", ns_vec.length());
    
    n2 = RcppArmadillo::sample(ns_vec, 1, FALSE)[0];
    // Rprintf("%i \n", n1);
    if (ns_vec.length() <= 2) {
      IntegerVector s_vec = Rcpp::Range(0, n - 1);
      IntegerVector n1_temp = s_vec[s_vec != n2];
      n1 = RcppArmadillo::sample(n1_temp, 1, FALSE)[0];
    } else {
      IntegerVector n1_temp = ns_vec[ns_vec != n2];
      n1 = RcppArmadillo::sample(n1_temp, 1, FALSE)[0];
    }
    
    // IntegerVector s_vec = Rcpp::Range(0, n - 1);
    // IntegerVector n1_temp = s_vec[s_vec != n2];
    // n1 = RcppArmadillo::sample(n1_temp, 1, FALSE)[0];
    
    
    arma::mat AW(dim_x, dim_y);
    AW.fill(0.0);
    
    for (int j = 0; j < n; ++j) {
      arma::mat K1 = K(n1, j);
      arma::mat K2 = K(n2, j);
      //arma::mat K1 = K(0, j);
      //arma::mat K2 = K(1, j);
      arma::mat diff_mat = (K2 - K1);
      //diff_mat.attr("dim") = Dimension(K1.nrow(), K1.ncol());
      //NumericMatrix tmp_coef(1, 1, a[j] * y[j] * y[n2]);
      double coef = (a[j] * y[j] * y[n2]);
      AW = AW + coef * diff_mat;
    }
    arma::mat K22 = K(n2, n2);
    arma::mat K11 = K(n1, n1);
    arma::mat K12 = K(n1, n2);
    arma::mat K21 = K(n2, n1);
    AA = K22 + K11 - K12 - K21; 
    
    S = y[n1] * y[n2];
    tempSum1 = 1.0 - S - accu(AW % W) / std::sqrt(f_norm2);
    tempSum2 = -1.0 * (f_norm2 * (accu(AA % W) + accu(AW % (AW + trans(AW)))) - 2 * std::pow(accu(AW % W), 2)) / std::pow(f_norm2, 1.5);
    
    maxDiff = tempSum1 / ( tempSum2 + 1e-6);
    
    if (std::abs(maxDiff) < 1e-8) {
      continue;
    }
    
    // if (abs(maxDiff) >= C) {
      //   maxDiff = C / 10.0;
      // }
    // Rprintf("%e \n", accu(AW));
    // Rprintf("%e \n", accu(AA % W));
    // Rprintf("%e \n", accu(AW % (AW + trans(AW))));
    // Rprintf("%e \n", pow(accu(AW % W), 2));
    // 
      // Rprintf("%e \n", -1.0 * (f_norm2 * (accu(AA % W) + accu(AW % (AW + trans(AW)))) - 2 * pow(accu(AW % W), 2)));
    // Rprintf("%e \n", pow(f_norm2, 1.5));
    // Rprintf("%e \n", maxDiff);
    a_old[n1] = a[n1];
    a_old[n2] = a[n2];
    
    a2_new = a[n2] - maxDiff;
    
    if (S == -1) {
      U = std::max(0.0, a_old[n2] - a_old[n1]);
      V = std::min(C, C - a_old[n1] + a_old[n2]);
    } else {
      U = std::max(0.0, a_old[n2] + a_old[n1] - C);
      V = std::min(C, a_old[n1] + a_old[n2]);
    }
    
    if (a2_new > V) {
      a2_new = V;
    }
    if (a2_new < U) {
      a2_new = U;
    }
    
    a1_new = a_old[n1] + S * (a_old[n2] - a2_new);
    a[n1] = a1_new;
    a[n2] = a2_new;
    
    maxDiff = -(a2_new - a_old[n2]);
    
    W = W + std::pow(maxDiff, 2) * AA - maxDiff * (AW + trans(AW));
    AW = AW - maxDiff * AA;
    
    for (int k = 0; k < n; ++k) {
      arma::mat K1 = K(n1, k);
      arma::mat K2 = K(n2, k);
      arma::mat diff_mat = (K2 - K1);
      double coef_new = (y[n2] * (a2_new - a_old[n2]));
      arma::mat sd_mat = sd[k];
      sd[k] = sd_mat + coef_new * diff_mat;
    }
    
    arma::mat b1_mat(dim_x, dim_y);
    b1_mat.fill(0.0);
    arma::mat b2_mat(dim_x, dim_y);
    b2_mat.fill(0.0);
    
    for (int i = 0; i < n; ++i) {
      arma::mat K1 = K(i, n1);
      arma::mat K2 = K(i, n2);
      double coef = a[i] * y[i];
      b1_mat = b1_mat + coef * K1;
      b2_mat = b2_mat + coef * K2;
    }
    
    b1 = y[n1] - accu(b1_mat % W) / std::sqrt(f_norm2);
    b2 = y[n2] - accu(b2_mat % W) / std::sqrt(f_norm2);
    
    if ((0.0 < a1_new) & (a1_new < C)) {
      b = b1;
    } else if ((0.0 < a2_new) & (a2_new < C)) {
      b = b2;
    } else {
      b = (b1 + b2) / 2;
    }
  }
  List out = List::create(a, b, W);
  return out;
}


