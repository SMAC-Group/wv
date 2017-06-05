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

#ifndef WV_COV
#define WV_COV

#include <RcppArmadillo.h>

arma::mat compute_cov_cpp(arma::field<arma::vec> coef1, arma::field<arma::vec> coef2, arma::vec variance, arma::vec lower, arma::vec upper);

#endif
