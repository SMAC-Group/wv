
#include "wv_cov.h"
#include "dwt.h"

// [[Rcpp::export]]
arma::mat compute_cov_cpp(arma::field<arma::vec> coef1, arma::field<arma::vec> coef2, arma::vec variance, arma::vec lower, arma::vec upper){

  unsigned int J = coef1.n_elem; 
  
  arma::mat V(J, 4); //1st column, covariance
                     //2nd & 3rd column, confidence interval 
  
  V.col(1) = variance; 
  
  for(unsigned int i = 0; i < J; i++){
    V(i, 0) = mean(coef1(i) % coef2(i));
    V(i, 2) = V(i, 0) + lower(i); 
    V(i, 3) = V(i, 0) + upper(i);
  }
  
  return V; 
}