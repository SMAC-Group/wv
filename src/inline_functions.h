
#ifndef INLINE_FUNCTIONS
#define INLINE_FUNCTIONS

inline double square(double x){
  return x*x;
}

inline arma::mat crossprod(const arma::mat& x){
    return x.t() * x;
}

inline arma::mat crossprod(const arma::mat& x, const arma::mat& y){
    return x.t() * y;
}


inline arma::mat tcrossprod(const arma::mat& x){
    return x * x.t();
}

inline arma::mat tcrossprod(const arma::mat& x, const arma::mat& y){
    return x * y.t();
}

#endif
