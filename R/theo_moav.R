# Copyright (C) 2017 James Balamuta, Justin Lee, Stephane Guerrier, Roberto Molinari
#
# This file is part of wv R Methods Package
#
# The `wv` R package is free software: you can redistribute it and/or modify it
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

#' @title Non-stationary Maximal-overlapping Allan Variance
#' @description 
#' Calculation of the theoritical Maximal-overlapping Allan Variance for constant-mean non-stationary time series data
#' @export
#' @param n        A \code{integer} indicating the length of each vector of consecutive observations considered for the average. 
#' @param covmat   A \code{matrix} indicating the T-by-T covariance matrix of the time series with length T.
#' @return A \code{field<numeric>} that is the theoritical Maximal-overlapping Allan Variance for constant-mean non-stationary time series data
#' @details
#' Calculate Maximal-overlapping Allan Variance based on the defination on "Xu, Haotian, et al. "A Study of the Allan Variance for Constant-Mean Non-Stationary Processes." arXiv preprint arXiv:1702.07795 (2017)".
#' n is an integer larger than 1 and smaller than \eqn{floor\left(log_2 \left(dim\left(covmat\right)[1]\right)\right)-1}{floor(log2(dim(T)[1]))-1}
#' @author Haotian Xu
#' @examples
#' set.seed(999)
#' AR1 phi=0.3
#' x = gen.gts(AR1(phi = 0.3, sigma2 = 1), N = 1000, freq = 1)
#' avar(x, type = "to")
#' 
#' a = matrix(rep(0, 1000^2), nrow = 1000)
#' for (i in 1:1000){
#'   a[,i] = seq(from = 1 - i, length.out = 1000)
#' }
#' a.diag = diag(a)
#' a[upper.tri(a,diag=T)] = 0
#' a = a + t(a) + diag(a.diag)
#' covmat = 0.3^a
#' sapply(1:8, function(y){MOAV(2^y, covmat)})
MOAV = function(n, covmat){
  # k: starting position of filter
  # n: length of each block
  # covmat: T by T covariance matrix
  # T: number of obversations
  mat = matrix(rep(0, n^2), nrow = n)
  T = dim(covmat)[1]
  for (k in 1:(T-2*n+1)){
    v1 = seq.int(from = k, to = k+n-1)
    v2 = seq.int(from = k+n, to = k+2*n-1)
    
    sigma_mat1 = covmat[v1, v1]
    sigma_mat2 = covmat[v2, v2]
    gamma_mat = covmat[v1, v2]
    
    mat = mat + sigma_mat1 + sigma_mat2 - 2*gamma_mat
  }
  sum(mat)/(n^2*2*(T-2*n+1))
}