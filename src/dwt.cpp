
#include "dwt.h"

// We use QMF, Haar, Select filter, etc.
#include "wv_filters.h"

// We use reverse vec
#include "armadillo_manipulations.h"

/* --------------------- Start DWT and MODWT Functions --------------------- */


//' @title Discrete Wavelet Transform
//' @description Calculation of the coefficients for the discrete wavelet transformation.
//' @param x           A \code{vector} with dimensions \eqn{N\times 1}{N x 1}.
//' @param filter_name A \code{string} indicating the filter.
//' @param nlevels     An \code{integer}, \eqn{J}, indicating the level of the decomposition.
//' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
//' @details
//' Performs a level J decomposition of the time series using the pyramid algorithm
//' @author James Balamuta and Justin Lee
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name, unsigned int nlevels) {
  
  unsigned int N = x.n_elem;

  unsigned int J;
  if (nlevels == floor(log2(N))){
    J = nlevels-1; 
  }else{
    J = nlevels; 
  }
  
  unsigned int tau = pow(2.0,double(J));

  if(double(N)/double(tau) != floor(double(N)/double(tau))){
    Rcpp::stop("The supplied sample size ('x') must be divisible by 2^(nlevels). Either truncate or expand the number of samples.");
  }
  if(tau > N){
    Rcpp::stop("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");
  }

  arma::field<arma::vec> filter_info = select_filter(filter_name);

  int L = arma::as_scalar(filter_info(0));
  arma::vec h = filter_info(1); //check the pulls
  arma::vec g = filter_info(2);

  arma::field<arma::vec> y(J);

  for(unsigned int j = 0; j < J; j++) {

    unsigned int M = N/pow(2.0,double(j));
    unsigned int M_over_2 = M/2;

    arma::vec Wj(M_over_2);
    arma::vec Vj(M_over_2);

    for(unsigned t = 0; t < M_over_2; t++) {

      int u = 2*t + 1;

      double Wjt = h(0)*x(u);
      double Vjt = g(0)*x(u);

      for(int n = 1; n < L; n++){
        u -= 1;
        if(u < 0){
          u = M - 1;
        }
        Wjt += h(n)*x(u);
        Vjt += g(n)*x(u);
      }

      Wj[t] = Wjt;
      Vj[t] = Vjt;
    }

    y(j) = Wj;
    x = Vj;
  }
  
  return y;
}

//' @title Maximum Overlap Discrete Wavelet Transform
//' @description
//' Calculation of the coefficients for the discrete wavelet transformation
//' @inheritParams dwt_cpp
//' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
//' @keywords internal
//' @details
//' Performs a level J decomposition of the time series using the pyramid algorithm.
//' Use this implementation to supply custom parameters instead of modwt(x),
//' which serves as a wrapper function.
//' @author James Balamuta and Justin Lee
//' @keywords internal
// [[Rcpp::export]]
arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name, unsigned int nlevels){
  
  unsigned int N = x.n_elem;
  
  unsigned int J = nlevels;
  
  unsigned int tau = pow(2.0,double(J));
  
  if(tau > N) Rcpp::stop("The number of levels [ 2^(nlevels) ] exceeds sample size ('x'). Supply a lower number of levels.");
  
  arma::field<arma::vec> filter_info = select_filter(filter_name);
  
  int L = arma::as_scalar(filter_info(0));
  arma::vec ht = filter_info(1);
  arma::vec gt = filter_info(2);
  
  // modwt transform
  double transform_factor = sqrt(2.0);
  ht /= transform_factor;
  gt /= transform_factor;
  
  arma::field<arma::vec> y(J);
  
  arma::vec Wj(N);
  arma::vec Vj(N);
  
  // changing index throughout the loop 
  int index = 0; 
  
  for(unsigned int j = 0; j < J; j++) {
    for(unsigned int t = 2*j+1; t < N; t++) {
      
      int k = t; 
      
      double Wjt = ht(0)*x(k);
      double Vjt = gt(0)*x(k);
      
      for(int n = 1; n < L; n++){
        k -= pow(2.0,double(j));
        if(k < 0){
          k += N;
        }
        Wjt += ht(n)*x(k);
        Vjt += gt(n)*x(k);
      }
      
      Wj[t] = Wjt;
      Vj[t] = Vjt;
    }
    index = 2*index+1;
    y(j) = Wj.rows(index,N-1);
    x = Vj;
  }
  
  return y;
}

/* -------------------------------- END DWT and MODWT Functions ---------------------- */
