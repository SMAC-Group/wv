/* Copyright (C) 2017 James Balamuta, Justin Lee, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of wv R Methods Package
 *
 * The `wv` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `wv` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

#include <RcppArmadillo.h>
#include "spatial_filter.h"
#include "robust_components.h"
#include "psi_functions.h"
#include "inline_functions.h"
#include "utility_functions.h"

//' Compute the Spatial Wavelet Coefficients
//' @param X      is a matrix with row, col orientation
//' @param J1,J2  is the levels of decomposition along the rows, columns
//' @export
//' @return A \code{list} of \code{vectors} containing the wavelet coefficients.
//' @details 
//' By default this function will return the wavelet coefficient in
//' addition to the wavelet
// [[Rcpp::export]]
arma::field<arma::vec> sp_wv_coeffs(const arma::mat& X, 
                                    int J1, int J2) {
    // number of rows
    int n = X.n_rows;
    // number of columns
    int m = X.n_cols;

    // Levels
    int nb_level = J1*J2;
    
    arma::field<arma::vec> wv_coeffs(nb_level);
    
    int i = 0;
    
    // Create a cache for hfilter coefficients
    int largestJ = std::max(J1, J2);
    
    arma::field<arma::vec> hfilters_cache(largestJ);

    for(int j = 1; j <= largestJ; ++j){
        hfilters_cache(j-1) = hfilter(j);
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
            
            arma::mat xh(n, n), xhh(n, n);
            
            // The authors behind this function make a large assumption
            // that the matrix contains NA values that can be removed
            // if an element is not filled. Highly inefficient. 
            xh.fill(arma::datum::nan); xhh.fill(arma::datum::nan);

            for(int spt = 0; spt < n; ++spt) {
                for(int tpt = 0; tpt < pp1; ++tpt) {
                    
                    arma::vec xts = X.submat(       tpt,       spt, tpt + mm1 - 1, spt );
                    // first_row, first_col, last_row, last_col
                    
                    // WV coefficients in one direction
                    xh( (tpt+ mm1/2), spt ) = sum(xts % hfil1) ; 
                    // WARNING: INTEGER DIVISION WAS USED HERE!
                }
            }
            
            for(int tpt = 0; tpt < pp1; ++tpt) {
                for(int spt = 0; spt < pp2; ++spt) {
                    
                    arma::rowvec xts = xh.submat(tpt + mm1 / 2,       spt, tpt + mm1 / 2, spt + mm2 - 1 );
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


//' Spatial Wavelet Variance
//' 
//' Computes the Spatial Wavelet Variance
//' @param wv_coeffs is \code{field<vec>} containing the wavelet coefficients
//' @param n,m is \code{int} containing the number of observations in rows, cols
//' @param iso    is whether the matrix is isometric
//' @inheritParams sp_wv_coeffs
//' @param robust \code{bool} is an indicator as to whether a classic or robust
//'               estimation should occur.
//' @param eff    \code{double} is the level of efficiency
// [[Rcpp::export]]
arma::field<arma::mat> sp_wvar_cpp(const arma::field<arma::vec>& wv_coeffs,
                                   int n, int m,
                                   int J1, int J2,
                                   bool iso = true,
                                   bool robust = true,
                                   double eff = 0.6){
    
    int i = 0;
    
    // This will be the last value of the optimizer.
    double crob_bw = 0;
    
    // WV Matrix 
    arma::mat wv(J1, J2);
    arma::mat wv_rob(J1, J2);
    
    // Compute the different WVs
    for(int j1 = 1; j1 <= J1; ++j1){
        for(int j2 = 1; j2 <= J2; ++j2){
            
            arma::vec wv_coeff = wv_coeffs(i);
            
            // Update the level
            i = i + 1;
            
            int mm1 = std::pow(2.0, double(j1) );
            int mm2 = std::pow(2.0, double(j2) );
            int pp1 = n - mm1 + 1;
            int pp2 = m - mm2 + 1;
            
            
            if (robust) { // Robust
                arma::vec temp = sig_rob_bw(wv_coeff, eff);
                
                wv(j1 - 1, j2 - 1) = temp(0);
                crob_bw = temp(1);
            } else {             // Standard
                wv(j1 - 1, j2 - 1) = as_scalar(sum(pow(wv_coeff, 2.0)))/double(pp1*pp2);
            }
        }
    }
    
    // Obtain the tuning constant
    if (!robust){
        arma::vec wv_coeff = wv_coeffs(i - 1);
        
        // Fix efficiency at 1.0
        arma::vec temp = sig_rob_bw(wv_coeff, 1.0);
        
        crob_bw = temp(1);
    }
    
    // Add ISO modification here
    if (iso == true) {   
        
        // redefine number of elements
        int nb_level = J1*(J2 + 1.0)/2.0;
        
        // Not a square matrix (Must have more rows than columns)
        if (wv.n_rows < wv.n_cols) { 
            wv = wv.t();
        }
        
        int min_dim = std::min(J1, J2);    
        
        wv = make_wv_iso(wv, min_dim - 1);
    }
    
    
    arma::field<arma::mat> values(2);
    arma::vec tuning(1);
    tuning(0) = crob_bw;
    
    values(0) = wv;
    values(1) = tuning;
    
    return values;
}


// Compute the Spatial Wavelet Variance Covariance Matrix
// @param wv_coeffs \code{field} of \code{vec} containing coefficients
// @param wv_rob    \code{matrix} of robust wv
// @param J1,J2     \code{int} of level of decomposition
// @param crob_bw   \code{double} value of tuning constant
// @param iso       \code{bool} isometry
arma::mat sp_wv_cov(const arma::field<arma::vec >& wv_coeffs,
                    const arma::mat& wv_rob,
                    int J1, int J2, double crob_bw, bool iso = true) {
    
    arma::mat V_rob;
    
    int i = 0, k = 0;
    
    // Only consider lower.tri part (symmetric)
    if (iso == true) {
        
        // Compute WV covariance matrix
        
        int nJJ = J1*J2; 
        
        int nstar = wv_coeffs(nJJ - 1).n_elem;
        
        arma::mat wv_coeffs_cut(nJJ, nstar);
        
        for(int j = 0; j < nJJ; ++j){
            arma::vec temp = wv_coeffs(j);
            wv_coeffs_cut.row(j) = temp.rows(0, nstar - 1).t(); 
        }
        
        // Robust
        arma::mat S_coeffs(wv_coeffs_cut.n_rows, wv_coeffs_cut.n_cols);
        
        for (int j1 = 0; j1 < J1; ++j1) {
            for (int j2 = 0; j2 < J2; ++j2) {
                // Force to row form
                S_coeffs.row(i) = psi_tuk(wv_coeffs_cut.row(i).t(), wv_rob(j1, j2), crob_bw).t();
                i = i + 1;
            }
        }
        
        arma::umat index = index_mat(J1, J2, "row");
        //matrix(1:(J1*J2),J1,J2,byrow=T);
        
        // Two S's that are the difference.
        arma::umat temp_coeffs = index;
        
        // bool all_zero = all( Y.elem(find(trimatl(Y))) == 0 );
        
        // Remove 1 to account for C++ index.
        arma::uvec index_row = temp_coeffs.elem( find( trimatl(temp_coeffs) ) ) - 1; 
                               // S.coeffs[index[lower.tri(index,diag = T)],]
        
        // Keep an eye on this condition
        arma::uvec index_col = temp_coeffs.elem( find( trimatu(temp_coeffs.t() ) ) ) - 1; 
                               // S.coeffs[t(index)[lower.tri(index,diag = T)],]
        // Deviates from subsetting routines because the transpose does not respect
        // the lower tri view. That is, in order to extract the correct elements
        // one needs to use the upper tri view.
        
        arma::mat S_low = S_coeffs.rows( index_row ); 
        arma::mat S_up = S_coeffs.rows( index_col ); 
        arma::mat S_lower = tcrossprod(S_low.col(0));
        arma::mat S_upper = tcrossprod(S_up.col(0));
        
        for(i = 1; i < nstar; ++i){
            S_lower = S_lower + tcrossprod(S_low.col(i));
            S_upper = S_upper + tcrossprod(S_up.col(i));
        }
        
        // Take the average
        arma::mat S = (S_lower + S_upper)/(2.0*nstar);
        
        // Same w.r.t M
        arma::mat M = arma::zeros<arma::mat>(nJJ, nJJ);
        for (k = 0; k < nstar; ++k) {
            
            arma::vec Mtemp = arma::zeros<arma::vec>(nJJ);
            i = 0;
            for(int j1 = 0; j1 < J1; ++j1){
                for(int j2 = 0; j2 < J2; ++j2){
                    
                    double deriv = der_psi_tuk(wv_coeffs_cut(i, k),
                                               wv_rob(j1, j2),
                                               crob_bw);
                    
                    Mtemp(i) = (arma::is_finite(deriv) ? deriv : 0);
                    
                    // Increment after the fact
                    i = i + 1;
                }
            }
            M = M + diagmat(Mtemp);
        }
        
        arma::mat M_lower = subset_matrix_parallel(M, index_col);
        arma::mat M_upper = subset_matrix_parallel(M, index_row);
        
        // Same w.r.t M regarding average
        M = ( M_lower + M_upper ) / 2.0;
        // (M[index[lower.tri(index,diag = T)], index[lower.tri(index,diag = T)]] +
        //    M[t(index)[lower.tri(index,diag = T)], t(index)[lower.tri(index,diag = T)]])/2;
        
        M = pinv(-M / nstar, 0.0, "std");
        
        V_rob = (M*S*M.t())/nstar;
        
    } else {
        
        int nb_level = J1*J2;
        
        int nstar = wv_coeffs(nb_level - 1).n_elem;
        
        arma::mat wv_coeffs_cut(nb_level, nstar);
        
        for(int j = 0; j < nb_level; ++j){
            arma::vec temp = wv_coeffs(j);
            wv_coeffs_cut.row(j) = temp.rows(0, nstar - 1).t(); 
        }
        
        
        //  Robust
        arma::mat S_coeffs( wv_coeffs_cut.n_rows, wv_coeffs_cut.n_cols);
        i = 0;
        for(int j1 = 0; j1 < J1; ++j1) {
            for(int j2 = 0; j2 < J2; ++j2){
                S_coeffs.row(i) = psi_tuk(wv_coeffs_cut.row(i).t(), wv_rob(j1,j2), crob_bw).t();
                i = i + 1;
                // Increment after the fact
            }
        }
        arma::mat S = tcrossprod(S_coeffs.col(0));
        for(i = 1; i < nstar; ++i){
            S = S + tcrossprod(S_coeffs.col(i));
        }
        
        S = S/nstar;
        
        arma::mat M = arma::zeros<arma::mat>(nb_level, nb_level);
        
        for(k = 0; k < nstar; ++k){
            
            arma::vec Mtemp = arma::zeros<arma::vec>(nb_level);
            i = 0;
            
            for(int j1 = 0; j1 < J1; ++j1){
                for(int j2 = 0; j2 < J2; ++j2){
                    double temp = der_psi_tuk(wv_coeffs_cut(i, k), wv_rob(j1,j2), crob_bw);
                    Mtemp(i) =  arma::is_finite(temp) ? temp : 0;
                    i = i + 1; // Increment after the fact
                }
            }
            
            M = M + diagmat(Mtemp);
        }
        
        M = pinv(-M / nstar, 0.0, "std");
        V_rob = (M*S*M.t())/nstar;
        
    }
    
    return V_rob;
}


//' Compute the Spatial Wavelet Variance
//' @inheritParams sp_wv_coeffs
//' @inheritParams sp_wvar_cpp
//' @return A \code{list} with the following values:
//' \itemize{
//' \item{Classical WV}
//' \item{Robust WV}
//' \item{Covariance Matrix}
//' }
//' @export
//' @seealso \code{\link{sp_wv_coeffs}}, \code{\link{spgmwm_exp}}
//' @examples
//' eps = 0 # no contamination
//' eps = 0.05 # with contamination
//' sig.eps = 9
//' m = n = 10
//' nlag = seq(1,n,1)
//' grid = expand.grid(nlag,nlag)
//' J1 = floor(log2(n)) - 1
//' J2 = floor(log2(m)) - 1
//' 
//' # Parameters
//' sig2 = 1
//' phi = 2
//' theta = c(phi,sig2)
//' wv.theo = exp_theo(theta,J1,J2)
//' p = length(theta)
//' alpha = 0.05
//' J = J1*(J2+1)/2
//' 
//' # Covariance matrix for simulation of process
//' distM = as.matrix(dist(grid))
//' Sigma = sig2*exp(-distM/phi)
//' cholSigma = t(chol(Sigma))
//' 
//' # Simulate the matrix
//' set.seed(234)
//' sim = cholSigma%*%rnorm(m*n)
//' 
//' # If contamination is set, contaminate
//' if(eps!=0){
//'     index = sample(1:(m*n),round((m*n)*eps))
//'     sim[index] = sim[index] + rnorm(length(index),0,sqrt(sig.eps))
//' }
//' 
//' Ymle = sim
//' 
//' Ygmwm = matrix(sim,n)
//' 
//' wv = spat_wavar(Ygmwm, J1, J2, eff = 0.6)
// [[Rcpp::export]]
arma::field<arma::mat> spat_wavar(const arma::mat& X, int J1, int J2,
                                  bool iso = true,
                                  bool robust = true, double eff = 0.6) {
    
    arma::field<arma::vec > coeff_results = sp_wv_coeffs(X, J1, J2);
    
    
    arma::field<arma::mat > wv_results = sp_wvar_cpp(coeff_results,
                                                     X.n_rows, X.n_cols,
                                                     J1, J2,
                                                     iso,
                                                     robust, eff);
    
    // Tuning constant
    double crob_bw = arma::as_scalar(wv_results(1));
    
    arma::mat V_rob = sp_wv_cov(coeff_results,
                                wv_results(0),
                                J1, J2, crob_bw, iso);
    
    arma::field<arma::mat> out(3);

    if (iso == true) {      // WV
        out(0) = lower_tri_elem(wv_results(0));
    } else {
        out(0) = vectorise(wv_results(0));
    }

    out(1) = wv_results(1); // crob_bw
    out(2) = V_rob;         // Covariance Matrix
    
    return out;
    
}