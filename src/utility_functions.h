#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H
arma::mat make_wv_iso(arma::mat wv, int min_dim);
arma::mat subset_matrix_parallel(const arma::mat& x, const arma::uvec y);
arma::umat index_mat(int m, int n, std::string ftype = "col");
arma::mat subset_matrix_parallel(const arma::mat& x, const arma::uvec y);
arma::vec lower_tri_elem(const arma::mat& X);
#endif