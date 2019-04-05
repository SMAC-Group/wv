
#include <RcppArmadillo.h>
#include "spatial_filter.h"

//' Compute the Spatial Wavelet Coefficients
//' @param X      is a matrix with row, col orientation
//' @param J1,J2  is the levels of decomposition along the rows, columns
//' @export
//' @return A \code{list} of \code{vectors} containing the wavelet coefficients.
//' @details 
//' By default this function will return the wavelet coefficient in
//' addition to the wavelet
// [[Rcpp::export]]
arma::field<arma::vec> sp_modwt_cpp(const arma::mat& X, int J1, int J2) {
  // number of rows
  int n = X.n_rows;
  // number of columns
  int m = X.n_cols;
  
  // Levels
  int nb_level = J1*J2;
  
  arma::field<arma::vec> wv_coeffs(nb_level);
  
  int i = 0;
  
  // Create a cache for sp_hfilter coefficients
  int largestJ = std::max(J1, J2);
  
  arma::field<arma::vec> hfilters_cache(largestJ);
  
  for(int j = 1; j <= largestJ; ++j){
    hfilters_cache(j-1) = sp_hfilter(j);
  }
  
  // Compute the different WVs
  for(int j1 = 1; j1 <= J1; ++j1){
    
    arma::vec hfil1 = hfilters_cache(j1-1);
    
    for(int j2 = 1; j2 <= J2; ++j2){
      
      arma::vec hfil2 = hfilters_cache(j2-1);
      
      int mm1 = std::pow(2.0, double(j1) );
      int mm2 = std::pow(2.0, double(j2) );
      int pp1 = n - mm1 + 1;
      int pp2 = m - mm2 + 1;
      
      int max = std::max(n, m); 
      
      arma::mat xh(max, max), xhh(max, max); 
      
      // NA filling Highly inefficient. 
      xh.fill(arma::datum::nan); xhh.fill(arma::datum::nan);
      
      for(int spt = 0; spt < n; ++spt) {
        for(int tpt = 0; tpt < pp1; ++tpt) {
          
          arma::vec xts = X.submat( tpt, spt, tpt + mm1 - 1, spt );
          // first_row, first_col, last_row, last_col
          
          // WV coefficients in one direction
          xh( (tpt+ mm1/2), spt ) = sum(xts % hfil1) ; 
          // WARNING: INTEGER DIVISION WAS USED HERE!
        }
      }
      
      for(int tpt = 0; tpt < pp1; ++tpt) {
        for(int spt = 0; spt < pp2; ++spt) {
          
          arma::rowvec xts = xh.submat(tpt + mm1 / 2, spt, tpt + mm1 / 2, spt + mm2 - 1 );
          // first_row, first_col, last_row, last_col
          
          // XHH are the wavelet coefficients at that level
          xhh( tpt + mm1/2, spt + mm2/2) = arma::as_scalar(xts * hfil2); // already summed.
          // WARNING: INTEGER DIVISION WAS USED HERE!
        }
      }
      
      arma::vec wv_coeff = xhh.elem( find_finite(xhh) ); //vectorise(xhh);
      
      wv_coeffs(i) = wv_coeff;
      
      // Update the level
      i = i + 1;
    }
    
  }
  
  return wv_coeffs;
}
