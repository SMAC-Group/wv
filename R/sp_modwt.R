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

#' @title Maximum Overlap Discrete Wavelet Transform
#' @description 
#' Calculates the coefficients for the discrete wavelet transformation
#' @param x        A \code{matrix} with dimensions N x M. 
#' @param J1       A \code{integer} indicating the \eqn{J1} row levels of decomposition.
#' @param J2       A \code{integer} indicating the \eqn{J2} column levels of decomposition
#' @return A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @author Justin Lee
#' @export
sp_modwt = function(x, J1 = floor(log2(dim(X)[1]-1)), J2 = floor(log2(dim(X)[2]-1))) {
  ret = sp_modwt_cpp(x = x, J1 = J1, J2 = J2)
  # mostattributes(ret) = list(J=nlevels, filter = filter, class=c("modwt","list"))
  ret
}