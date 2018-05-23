# Copyright (C) 2017 James Balamuta, Justin Lee, Stephane Guerrier, Roberto Molinari
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


#' @title Calculate Theoretical Allan Variance for Stationary First-Order Autoregressive 
#' (AR1) Process
#' @description 
#' This function allows us to calculate the theoretical allan variance for stationary 
#' first-order autoregressive (AR1) process.
#' @export
#' @usage av_ar1(n, phi, sigma2)
#' @param n An \code{integer} value for the size of the cluster.
#' @param phi A \code{double} value for the autocorrection parameter \eqn{\phi}{phi}.
#' @param sigma2 A \code{double} value for the variance parameter \eqn{\sigma ^2}{sigma^2}.
#' @return A \code{double} indicating the theoretical allan variance for AR1 process.
#' @note This function is based on the calculation of the theoretical allan variance 
#' for stationary AR1 process raised in "Allan Variance of Time Series Models for 
#' Measurement Data" by Nien Fan Zhang. (For more details: 
#' \url{http://iopscience.iop.org/article/10.1088/0026-1394/45/5/009/meta}.) This calculation
#' is fundamental and necessary for the study in "A Study of the Allan Variance for Constant-Mean 
#' Non-Stationary Processes" by Xu et al. (IEEE Signal Processing Letters, 2017), preprint available:
#' \url{https://arxiv.org/abs/1702.07795}.
#' @author Yuming Zhang
#' @examples
#' av1 = av_ar1(n = 5, phi = 0.9, sigma2 = 1)
#' av2 = av_ar1(n = 8, phi = 0.5, sigma2 = 2)
av_ar1 = function(n, phi, sigma2){
  numerator = n-3*phi-n*phi^2+4*phi^(n+1)-phi^(2*n+1)
  denominator = n^2 * (1-phi)^2 * (1-phi^2)
  result = numerator / denominator * sigma2
  return(result)
}




#' @title Calculate Theoretical Allan Variance for Stationary White Noise Process
#' @description 
#' This function allows us to calculate the theoretical allan variance for stationary 
#' white noise process.
#' @export
#' @usage av_wn(sigma2, n)
#' @param sigma2 A \code{double} value for the variance parameter \eqn{\sigma ^2}{sigma^2}.
#' @param n An \code{integer} value for the size of the cluster.
#' @return A \code{double} indicating the theoretical allan variance for the white noise
#' process.
#' @note This function is based on the calculation of the theoretical allan variance 
#' for stationary white noise process raised in "Allan Variance of Time Series Models for 
#' Measurement Data" by Nien Fan Zhang. (For more details: 
#' \url{http://iopscience.iop.org/article/10.1088/0026-1394/45/5/009/meta}.) This calculation
#' is fundamental and necessary for the study in "A Study of the Allan Variance for Constant-Mean 
#' Non-Stationary Processes" by Xu et al. (IEEE Signal Processing Letters, 2017), preprint available:
#' \url{https://arxiv.org/abs/1702.07795}.
#' @author Yuming Zhang
#' @examples
#' av1 = av_wn(sigma2 = 1, n = 5)
#' av2 = av_wn(sigma2 = 2, n = 8)
av_wn = function(sigma2, n){
  result = sigma2/n
  return(result)
}