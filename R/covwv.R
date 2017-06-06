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

#' @title Pair Cross Covariance 
#' @description
#' Calculates the Cross-covariance between two wavelet transfomations (dwt or modwt)
#' @param x         A \code{vector} with dimensions N x 1.
#' @param y         A \code{vector} with dimensions N x 1.
#' @param decomp    A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition.
#' @param filter    A \code{string} that specifies what wavelet filter to use.
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @return Returns cross-covariance, variance of each wavelet cross-covariance and its 95% CI.
#' @importFrom coda spectrum0
#' @importFrom stats qnorm
#' @export
#' @details 
#' If \code{nlevels} is not specified, it is set to \eqn{\left\lfloor {{{\log }_2}\left( {length\left( x \right)} \right)} \right\rfloor}{floor(log2(length(x)))}
#' @author Justin Lee
crosswvar_pair = function(x, y, decomp = "modwt", filter = "haar", nlevels = NULL){
  if(is.null(x) || is.null(y)){
    stop("`x` or `y` must contain a value.")
  }else if((is.data.frame(x) || is.matrix(x)) || is.data.frame(y) || is.matrix(y)){
    if(ncol(x) > 1) stop("There must be only one column of data supplied.")
  }
  
  if(length(x)!=length(y)){
    stop("`x` and `y` must be of same length.")
  }
  
  if(decomp == "modwt"){
    if(is.null(nlevels)){
      nlevels = floor(log2(length(x))-1)
    }
    f = modwt
  }else if(decomp == "dwt"){
    if(is.null(nlevels)){
      nlevels = floor(log2(length(x)))
    }
    f = dwt
  }
  
  coef1 = f(x = x, nlevels = nlevels, filter = filter)
  coef2 = f(x = y, nlevels = nlevels, filter = filter)
  
  # Slightly inefficient!
  # Better method may be to implement spectrum0 in Rcpp 
  product = mapply("*", coef1, coef2, SIMPLIFY = FALSE) 
  const = 1/sapply(product, length)
  variance = const * unlist(sapply(product, spectrum0, max.freq = 0, order = 2, max.length = 130))
  lower = sqrt(unname(variance)) * qnorm(0.025)
  upper = sqrt(unname(variance)) * qnorm(0.975)
  
  obj =  .Call('wv_compute_cov_cpp', PACKAGE = 'wv', coef1 = coef1, coef2 = coef2, variance = variance, lower = lower, upper = upper) 
  
  # colnames(ret) = c("Cross-Covariance", "Variance", "Lower Bound", "Upper Bound")
  # ret = list(obj)
  
  return(list(obj))
}

#' @export
#' @author Justin Lee
crosswvar = function(x, decomp = "modwt", filter = "haar", nlevels = NULL){
  if(is.null(x)) stop("`x` must contain a value.")
  
  if(decomp == "modwt"){
    if(is.null(nlevels)){
      nlevels = floor(log2(nrow(x))-1)
    }
  }else if(decomp == "dwt"){
    if(is.null(nlevels)){
      nlevels = floor(log2(nrow(x)))
    }
  }
  
  mat = matrix(list(), ncol(x), ncol(x))
  
  for(i in seq_len(ncol(x))){
    j = i
    for(j in i:ncol(x)){
      mat[i,j] = crosswvar_pair(x[,i], x[,j], decomp = decomp, filter = filter, nlevels = nlevels)
    }
  }
  
  return(mat)
}