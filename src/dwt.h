#ifndef DWT
#define DWT

#include <RcppArmadillo.h>

arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name = "haar", unsigned int nlevels = 4);

arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name = "haar", unsigned int nlevels = 4);

#endif
