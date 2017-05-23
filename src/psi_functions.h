#ifndef PSI_FUNC_H
#define PSI_FUNC_H
arma::vec psi_tuk(const arma::vec& x, double sig2_bw, double crob_bw);
arma::vec der_psi_tuk(const arma::vec& x, double sig2_bw, double crob_bw);
double der_psi_tuk(double x, double sig2_bw, double crob_bw);
#endif
