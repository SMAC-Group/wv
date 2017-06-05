# Copyright (C) 2017 Justin Lee, James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of wv R Methods Package
#
# The `wv` R package is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The `wv` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Cross-covariance 
#' @description
#' Calculates the Cross-covariance between two wavelet transfomations (dwt or modwt)
#' @param x   A \code{vector} with dimensions N x 1.
#' @param y   A \code{vector} with dimensions N x 1.
#' @param decomp    A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition.
#' @param filter    A \code{string} that specifies what wavelet filter to use.
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @return
#' @details 
#' If \code{nlevels} is not specified, it is set to \eqn{\left\lfloor {{{\log }_2}\left( {length\left( x \right)} \right)} \right\rfloor}{floor(log2(length(x)))}
#' @author Justin Lee
crosswvar = function(x, y, decomp = "modwt", filter = "haar", nlevels = NULL){
  if(is.null(x) || is.null(y)){
    stop("`x` or `y` must contain a value.")
  }else if((is.data.frame(x) || is.matrix(x)) || is.data.frame(y) || is.matrix(y)){
    if(ncol(x) > 1) stop("There must be only one column of data supplied.")
  }
  
  if(length(x)!=length(y)){
    stop("`x` and `y` must be of same length.")
  }
  
  if(decomp == "modwt" && is.null(nlevels)){
    nlevels = floor(log2(length(x))-1)
  }else if(decomp == "dwt" && is.null(nlevels)){
    nlevels = floor(log2(length(x)))
  }
  
  obj =  .Call('wv_compute_cov_cpp', PACKAGE = 'wv', x = x, y = y, decomp = decomp, filter = filter, nlevels = nlevels)  
}