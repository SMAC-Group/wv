
#ifndef WV_COV
#define WV_COV

#include <RcppArmadillo.h>

arma::mat compute_cov_cpp(arma::field<arma::vec> coef1, arma::field<arma::vec> coef2, arma::vec variance, arma::vec lower, arma::vec upper);

#endif
