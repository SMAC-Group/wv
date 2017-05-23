/* Copyright (C) 2017 James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of `spgmwm` R Package
 *
 * The `spgmwm` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `spgmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

#include <RcppArmadillo.h>
#include "psi_functions.h"

//' Psi Tukey Scoring function
//' @param x data vector
//' @param sig2_bw variance of brickwalled wv coefficients
//' @param crob_bw ?
//' @return A \code{vec} containing the scores.
//' @export
// [[Rcpp::export]]
arma::vec psi_tuk(const arma::vec& x, double sig2_bw, double crob_bw){
    arma::vec r = x/sqrt(sig2_bw);
    arma::vec w = (abs(r) < crob_bw) % pow(1 - pow(r/crob_bw, 2.0), 2.0);
    
    double crob_bw_2d = crob_bw*crob_bw;
    double crob_bw_4d = crob_bw_2d*crob_bw_2d;
    double crob_bw_6d = crob_bw_4d*crob_bw_2d;
    double crob_bw_8d = crob_bw_4d*crob_bw_4d;
    
    arma::vec out = pow(r, 2.0) % pow(w, 2.0) - (1.0/crob_bw_8d)*(
        1890*R::pnorm(crob_bw, 0, 1, true, false)-2.0*crob_bw*(
                945.0 + 315.0*crob_bw_2d + 63.0*crob_bw_4d + 9.0*crob_bw_6d + crob_bw_8d) * 
                R::dnorm(crob_bw, 0.0, 1.0, false) - 945.0) - (4.0/crob_bw_6d)*
                (210.0*R::pnorm(crob_bw, 0, 1, true, false) - 2.0*crob_bw*
                (105.0 + 35.0*crob_bw_2d+7*crob_bw_4d+crob_bw_6d)*
                R::dnorm(crob_bw, 0.0, 1.0, false)-105.0)+(6.0/crob_bw_4d)*
                (30.0*R::pnorm(crob_bw, 0, 1, true, false)-2.0*crob_bw*
                (15.0+5.0*crob_bw_2d+crob_bw_4d)*R::dnorm(crob_bw, 0.0, 1.0, false)-15.0) 
        - (4.0/crob_bw_2d)*(6.0*R::pnorm(crob_bw, 0, 1, true, false) 
                                - 2.0*crob_bw*(3.0+crob_bw_2d)*R::dnorm(crob_bw, 0.0, 1.0, false)-3.0)
                                + 2.0 *R::pnorm(crob_bw, 0, 1, true, false) 
                                - 2.0 *crob_bw*R::dnorm(crob_bw, 0.0, 1.0, false)
                                - 1;
    return out;
}

//' Derivative of Psi Tukey Scoring function
//' @inheritParams psi_tuk
//' @export
// [[Rcpp::export]]
arma::vec der_psi_tuk(const arma::vec& x, double sig2_bw, double crob_bw){
    arma::vec r = x/sqrt(sig2_bw);
    double crob_bw_2d = crob_bw*crob_bw;

    arma::vec r2 = pow(pow(r, 2.0)/crob_bw_2d - 1.0, 3.0) % (pow(r, 2.0) / sig2_bw) % (1.0 - 5.0*(pow(r, 2.0)/crob_bw_2d));
    r2 = (abs(r) < crob_bw) % r2;
    return r2;
}

double der_psi_tuk(double x, double sig2_bw, double crob_bw){
    double r = x/sqrt(sig2_bw);
    double crob_bw_2d = crob_bw*crob_bw;
    
    double r2 = pow(pow(r, 2.0)/crob_bw_2d - 1.0, 3.0) * (pow(r, 2.0) / sig2_bw) * (1.0 - 5.0*(pow(r, 2.0)/crob_bw_2d));
    r2 = (std::abs(r) < crob_bw) * r2;
    return r2;
}
