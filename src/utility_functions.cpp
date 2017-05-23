#include <RcppArmadillo.h>
#include "utility_functions.h"

//' Create the ISO matrix
//' @param min_dim An \code{integer} indicating the minimum value.
//' @param wv      A \code{matrix} containing the wavelet variance.
//' @export
//' @details
//' Takes the average of the lower triangular view
//' before returning. 
//' @examples
//' a = matrix(1:9, nrow = 3, byrow = TRUE)
//' make_wv_iso(a, 3 - 1)
// [[Rcpp::export]]
arma::mat make_wv_iso(arma::mat wv, int min_dim) {
    
    arma::mat WV_temp = wv.submat(          0,         0,  min_dim, min_dim );
    // first_row, first_col, last_row, last_col
    
    // Average of off diagonal elements
    arma::mat WV_temp1 = (trimatl(WV_temp) + trimatl(WV_temp.t()))/2;
    
    // Put it back into matrix
    WV_temp = trimatl(WV_temp1) - diagmat(WV_temp) + trimatu(WV_temp) ;
    // Avoid double counting the diagonal and include the old upper triangle
    
    wv.submat(          0,         0,  min_dim, min_dim ) = WV_temp;
    // first_row, first_col, last_row, last_col
    
    return wv;
}

//' Perform parallel contingous subset
//' @param x \code{matrix} to use to subset
//' @param y \code{vector} of ids
//' @export
//' @examples
//' a = matrix(1:100, nrow = 10, byrow = TRUE)
//' 
//' idx = c(1, 3, 10)
//' subset_matrix_parallel(a, idx - 1)
//' ## same as
//' # a[c(1, 3, 10), c(1, 3, 10)]
// [[Rcpp::export]]
arma::mat subset_matrix_parallel(const arma::mat& x, const arma::uvec y){
    
    int n = y.n_elem;
    
    arma::mat out(n, n);
    
    for(int i = 0; i < n; ++i){
        int id = y(i);
        for(int j = 0; j < n; ++j){
            out(i,j) = x(id, y(j));   
        }
    }
    
    return out;
}


//' Create an index matrix
//' @param m,n   An \code{integer} indicating matrix row, column dimension
//' @param ftype A \code{string} of either \code{"col"} or \code{"row"}.
//' @export
//' @details
//' Creates an index matrix filled by either row or column.
//' @examples
//' # Index matrix by column
//' index_mat(3,4, "col")
//' # index matrix by row
//' index_mat(3,4, "row")
// [[Rcpp::export]]
arma::umat index_mat(int m, int n, std::string ftype){
    
    arma::umat ind = arma::umat(m, n);
    
    int count = 0;
    
    // fill by column by default
    if(ftype == "col"){
        for(int j = 0; j < n; ++j){
            for(int i = 0; i < m; ++i){
                ++count;
                ind(i, j) = count;
            }
        }
    } else { // fill by row
        for(int i = 0; i < m; ++i){
            for(int j = 0; j < n; ++j){
                ++count;
                ind(i, j) = count;
            }
        }
    }
    
    return ind;
}

//' Extract Lower Triangular Elements
//' @param X \code{matrix}
//' @export
//' @examples
//' x = matrix(1:16, 4, 4, byrow = TRUE)
//' lower_tri_elem(x)
// [[Rcpp::export]]
arma::vec lower_tri_elem(const arma::mat& X){
    return X.elem( find( trimatl(X) ) );
}
