# Copyright (C) 2017  Justin Lee, James Balamuta, Stephane Guerrier, Roberto Molinari
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

#' Calculate the Allan Variance
#' 
#' Computes the Allan Variance
#' @param x     A \code{vec} containing the time series under observation.
#' @param type  A \code{string} containing either \code{"mo"} for Maximal Overlap or \code{"to"} for Tau Overlap
#' @return av   A \code{list} that contains:
#' \itemize{
#'  \item{"clusters"}{The size of the cluster}
#'  \item{"allan"}{The Allan variance}
#'  \item{"errors"}{The error associated with the variance estimation.}
#' }
#' @details 
#' The decomposition and the amount of time it takes to perform it depends on whether you are using
#' the Tau Overlap or the Maximal Overlap.
#' 
#' @section Maximal Overlap Allan Variance:
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n | n < floor(log2(N))}}
#' Then, \eqn{M = N - 2n} samples exist. 
#' The Maximal-overlap estimator is given by:
#' \deqn{\frac{1}{{2\left( {N - 2k + 1} \right)}}\sum\limits_{t = 2k}^N {{{\left[ {{{\bar Y}_t}\left( k \right) - {{\bar Y}_{t - k}}\left( k \right)} \right]}^2}} }{See PDF Manual}
#' 
#' where \deqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }{See PDF Manual}.
#' 
#' @section Tau-Overlap Allan Variance:
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n | n < floor(log2(N))}}
#' Then, a sampling of \eqn{m = \left\lfloor {\frac{{N - 1}}{n}} \right\rfloor  - 1} samples exist. 
#' The tau-overlap estimator is given by:
#' 
#' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }{See PDF Manual}.
#' @export
#' @author James Balamuta, Stephane Gurrier
#' @references Long-Memory Processes, the Allan Variance and Wavelets, D. B. Percival and P. Guttorp
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#' 
#' # Simulate time series
#' N = 100000
#' ts = gen_gts(N, WN(sigma2 = 2) + RW(gamma2 = 1))
#' 
#' # Maximal overlap
#' av_mat_mo = avar(ts, type = "mo")
#' 
#' # Tau overlap
#' av_mat_tau = avar(ts, type = "to")
avar = function(x, type = "mo") {
  x = as.vector(x)
  
  if(type == "mo"){
    av = .Call('wv_avar_mo_cpp', PACKAGE = 'wv', x)
  }else{
    av = .Call('wv_avar_to_cpp', PACKAGE = 'wv', x)
  }
  
  av = list(clusters = av[,1], allan=av[,2], errors=av[,3])
  av$adev = sqrt(av$allan)
  av$lci = av$adev - av$errors*av$adev
  av$uci = av$adev + av$errors*av$adev
  av$type = type
  class(av) = c("avar","list")
  av
}

#' Prints Allan Variance
#' 
#' Displays the allan variance information
#' @method print avar
#' @export
#' @param x   A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @author James Balamuta
#' @return console output
#' @keywords internal
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#' 
#' # Generate time series
#' x = gen_gts(100, WN(sigma2 = 1))
#' 
#' # Compute Allan
#' out = avar(x)
#' 
#' # Print results
#' print( out )
print.avar = function(x, ...) {
  cat("\n Clusters: \n")
  print(x$clusters, digits=5)
  cat("\n Allan Variances: \n")
  print(x$allan, digits=5)
  cat("\n Errors: \n")
  print(x$errors, digits=5)
}

#' Summary Allan Variance
#' 
#' Displays the summary table of allan variance
#' @method summary avar
#' @export
#' @param object A \code{avar} object.
#' @param ...    Additional arguments affecting the summary produced.
#' @author James Balamuta
#' @return Summary table
#' @keywords internal
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#' 
#' # Generate time series
#' x = gen_gts(100, WN(sigma2 = 1))
#' 
#' # Compute Allan
#' out = avar(x)
#' 
#' # Summary
#' summary( out )
summary.avar = function(object, ...) {
  out_matrix = matrix(0, nrow = length(object$clusters), ncol = 6)
  colnames(out_matrix) = c("Time", "AVAR", "ADEV", "Lower CI", "Upper CI", "Error")
  out_matrix[,"Time"] = object$clusters
  out_matrix[,"AVAR"] = object$allan
  out_matrix[,"ADEV"] = object$adev
  out_matrix[,"Lower CI"] = object$lci
  out_matrix[,"Upper CI"] = object$uci
  out_matrix[,"Error"] = object$errors
  
  class(out_matrix) = c("summary.avar","matrix")
  out_matrix
}
