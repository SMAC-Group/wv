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

#' Maximum Overlap Discrete Wavelet Transform
#'  
#' Calculation of the coefficients for the discrete wavelet transformation
#' @param x        A \code{vector} with dimensions N x 1. 
#' @param nlevels  A \code{integer} indicating the \eqn{J} levels of decomposition.
#' @param filter   A \code{string} indicating the filter name
#' @return y       A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' The default \eqn{J} is determined by \eqn{floor\left(log_2 \left(length\left(x\right)\right)\right)}{floor(log2(length(x)))}
#' @author JJB, Justin 
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' a = modwt(x)
#' @export
modwt = function(x, nlevels = floor(log2(length(x))), filter = "haar") {
  out = modwt_cpp(x = x, filter_name = filter, nlevels)
  mostattributes(out) = list(J=nlevels, filter = filter, class=c("modwt","list"))
  out
}

#' @title Print Maximum Overlap Discrete Wavelet Transform
#' @description
#' Prints the results of the modwt list
#' @method print modwt
#' @export
#' @param x A \code{modwt} object
#' @param ... further arguments passed to or from other methods.
#' @return Prints the modwt decomposition
#' @author JJB
#' @keywords internal
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' print(modwt(x))
print.modwt = function(x, ...){
  NextMethod("print")
}

#' @title Summary Maximum Overlap Discrete Wavelet Transform
#' @description Unlists MODWT object and places it in matrix form
#' @method summary modwt
#' @export
#' @keywords internal
#' @param object A \code{modwt} object
#' @param ... additional arguments affecting the summary produced.
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' summary(modwt(x))
summary.modwt=function(object, ...) {
  cat("Results of MODWT using",attr(object,"filter"),"filter with",attr(object, "J")-1,"levels:\n")
  y = as.list(object)
  j = length(y)
  for( i in 1:j ) {
    cat("Level",i,"Wavelet Coefficients")
    y[[i]]
    cat("\n")
  }
}
