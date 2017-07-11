/* Copyright (C) 2017 James Balamuta, Justin Lee, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of `wv` R Package
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

//' Haar filter for a spatial case
//' @param jscale An \code{int} of the Number of Scales
//' @export
// [[Rcpp::export]]
arma::vec sp_hfilter(int jscale){
    
    // mother wavelet (1/sqrt(2))*c(1,1)/sqrt(2)
    arma::vec g(2); 
    g.fill(0.5); 
    // father wavelet (1/sqrt(2),-1/sqrt(2))/sqrt(2)
    arma::vec h(2);
    h.fill(0.5);
    
    // number of filters
    int L = 2;
    
    if (jscale == 1) { return h; }
        
    // Create a zero vector
    int nfill_z1 = pow(2.0, jscale - 1.0) - 1.0;
    
    // Number of elements
    arma::vec hup(L + nfill_z1);
        
    hup(0) = h(0);
    hup.rows(1, nfill_z1) = arma::zeros<arma::vec>( nfill_z1 );
    hup(hup.n_elem-1) = h(1);

    // We must process the mother wavelet for this function
        
    // Always initialized to g.
    arma::vec gup = g;
    
    for(int j = 1; j < (jscale - 1); ++j){
        
        // Amount to fill
        int nfill_z2 = pow(2.0, j) - 1.0 + L;
        
        arma::vec temp(nfill_z2);
        temp(0) = g(0);
        temp.rows(1, nfill_z2 - L) = arma::zeros<arma::vec>(nfill_z2 - L);
        temp(temp.n_elem - 1.0) = g(1);
        
        // Avoid the unsigned comparison
        int gelems = gup.n_elem;
        int telems = gelems + nfill_z2 - 1;
        
        arma::vec sala = arma::zeros<arma::vec>( telems );
        
        for (int k = 0; k < telems; ++k) {
            double dummy = 0; 
            
            for (int u = 0; u < gelems; ++u) {
                if ( (k - u) >= 0 && (k - u) < nfill_z2) {
                    dummy = dummy + gup(u)*temp(k - u);
                }
            }
            sala(k) = dummy;
        }
        
        gup = sala;
    }

    arma::vec sala = arma::zeros<arma::vec>(hup.n_elem + gup.n_elem - 1 );

    // Avoids unsigned comparison
    int selems = sala.n_elem;
    int gelems = gup.n_elem;
    int helems = hup.n_elem;    
    
    for (int k = 0; k < selems; ++k) {
        double dummy = 0;
        for (int u = 0; u < helems; ++u) {
            if( (k - u) >= 0 && (k - u) < gelems){
                dummy = dummy + hup(u)*gup(k - u);
            }
        }
        sala(k) = dummy;
    }

    // sala is the new hup
    return sala;
}