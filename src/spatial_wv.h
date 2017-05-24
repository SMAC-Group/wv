#ifndef SPATIAL_WV_H
#define SPATIAL_WV_H


arma::field<arma::vec> sp_wv_coeffs(const arma::mat& X, 
                                    int J1, int J2);

arma::field<arma::mat> sp_wvar_cpp(const arma::field<arma::vec>& wv_coeffs,
                                   int n, int m,
                                   int J1, int J2,
                                   bool iso = true,
                                   bool robust = true,
                                   double eff = 0.6);

arma::mat sp_wv_cov(const arma::field<arma::vec >& wv_coeffs,
                    const arma::mat& wv_rob,
                    int J1, int J2, double crob_bw, bool iso = true);

arma::field<arma::mat> spat_wavar(const arma::mat& X, int J1, int J2,
                                  bool iso = true,
                                  bool robust = true, double eff = 0.6);

#endif