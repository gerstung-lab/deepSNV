#include <iostream>
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;


#define CN 8
static double cof[CN]={
  2.5066282746310005,
  1.0000000000190015,
  76.18009172947146,
  -86.50532032941677,
  24.01409824083091,
  -1.231739572450155,
  0.1208650973866179e-2,
  -0.5395239384953e-5,
};


  // [[Rcpp::export]]
double c_log_gamma(double x) {
  double y,ser,*co;
  int j;

    ser=cof[1]; y=x; co=cof+2;
    for(j=2;j<CN;j++) {
      y+=1.; ser+=(*co)/y; co++;
    }

      y=x+5.5;
      y-=(x+0.5)*log(y);
      return(-y+log(cof[0]*ser/x));
}

// [[Rcpp::export]]
double c_log_beta(double x, double y) {
  double res;
  res = c_log_gamma(x) + c_log_gamma(y) - c_log_gamma(x+y);
  return res;
}

// [[Rcpp::export]]
double c_logbb(int x, int n, double mu, double disp) {
  double res;
  res = c_log_beta(x+mu,n-x-mu+disp) - c_log_beta(mu, disp - mu);
  return res;
}

// [[Rcpp::export]]
NumericVector c_bayesian_factor_part(IntegerVector x, IntegerVector n, NumericVector mu, NumericVector disp) {
  int i=0;
  int dim=x.size();
  NumericVector res(dim);
      for (i=0; i<dim; i++) {
            res[i] = c_logbb(x[i],n[i],mu[i],disp[i]);
  }
  return res;
}

// [[Rcpp::export]]
NumericVector c_bayesian_factor(IntegerVector x_fw, IntegerVector X_fw, IntegerVector x_bw, IntegerVector X_bw,
                                IntegerVector n_fw, IntegerVector N_fw, IntegerVector n_bw, IntegerVector N_bw,
                                NumericVector nu0_fw, NumericVector nu0_bw, NumericVector nu1s_both, NumericVector nu1c_fw, NumericVector nu1c_bw,
                                NumericVector rdisp) {
  int i=0;
  int dim=x_fw.size();
  NumericVector res(dim);

  for (i=0; i<dim; i++) {

    res[i] = c_logbb(x_fw[i], n_fw[i], nu0_fw[i], rdisp[i]) +
      c_logbb(X_fw[i], N_fw[i], nu0_fw[i], rdisp[i]) +
      c_logbb(x_bw[i], n_bw[i], nu0_bw[i], rdisp[i]) +
      c_logbb(X_bw[i], N_bw[i], nu0_bw[i], rdisp[i]) -
      c_logbb(x_fw[i], n_fw[i], nu1s_both[i], rdisp[i]) -
      c_logbb(X_fw[i], N_fw[i], nu1c_fw[i], rdisp[i]) -
      c_logbb(x_bw[i], n_bw[i], nu1s_both[i], rdisp[i]) -
      c_logbb(X_bw[i], N_bw[i], nu1c_bw[i], rdisp[i]);

  }
  return res;
}
