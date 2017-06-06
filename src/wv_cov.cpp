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

#include "wv_cov.h"
#include "dwt.h"

// [[Rcpp::export]]
arma::mat compute_cov_cpp(arma::field<arma::vec> coef1, arma::field<arma::vec> coef2, arma::vec variance, arma::vec lower, arma::vec upper){

  unsigned int J = coef1.n_elem; 
  
  arma::mat V(J, 4); //1st column, covariance
                     //2nd & 3rd column, confidence interval 
  
  V.col(1) = variance; 
  
  for(unsigned int i = 0; i < J; i++){
    V(i, 0) = mean(coef1(i) % coef2(i));
    V(i, 2) = V(i, 0) + lower(i); // not sure how CI is calculated 
    V(i, 3) = V(i, 0) + upper(i);
  }
  
  return V; 
}