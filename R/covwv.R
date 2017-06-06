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

#' @title Cross Covariance of a TS Pair 
#' @description
#' Calculates the Cross-covariance between two wavelet transfomations (dwt or modwt)
#' @param x         A \code{vector} with dimensions N x 1.
#' @param y         A \code{vector} with dimensions N x 1.
#' @param decomp    A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition.
#' @param filter    A \code{string} that specifies what wavelet filter to use.
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @return Returns a \code{list} of a \code{matrix} containing cross-covariance, variance of each wavelet cross-covariance and its 95% CI.
#' @importFrom coda spectrum0
#' @importFrom stats qnorm
#' @details 
#' If \code{nlevels} is not specified, it is set to \eqn{\left\lfloor {{{\log }_2}\left( {length\left( x \right)} \right)} \right\rfloor}{floor(log2(length(x)))}
#' @export
#' @author Justin Lee
wccv_pair = function(x, y, decomp = "modwt", filter = "haar", nlevels = NULL){
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
  
  colnames(obj) = c("Cross-Covariance", "Variance", "Lower Bound", "Upper Bound")
  ret = list(obj)
  
  return(ret)
}

#' @title Cross Covariance of Matrix 
#' @description
#' Calculates the Cross-covariance between multiple wavelet transfomations (dwt or modwt)
#' @param x         A \code{vector} with dimensions N x M.
#' @param decomp    A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition.
#' @param filter    A \code{string} that specifies what wavelet filter to use.
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @return Returns a \code{matrix} of \code{lists} of all the possible pair cross-covariance, variance of each wavelet cross-covariance and its 95% CI.
#' @details 
#' If \code{nlevels} is not specified, it is set to \eqn{\left\lfloor {{{\log }_2}\left( {length\left( x \right)} \right)} \right\rfloor}{floor(log2(length(x)))}
#' @export
#' @author Justin Lee
wccv = function(x, decomp = "modwt", filter = "haar", nlevels = NULL){
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
      mat[i,j] = wccv_pair(x[,i], x[,j], decomp = decomp, filter = filter, nlevels = nlevels)
    }
  }
  
  mostattributes(mat) = list(filter = filter, class=c("wccv","matrix","list"))
  
  return(mat)
}

#' @title Plot Cross Covariance 
#' @description
#' Plots results of the wccv list in which additional parameters can be specified
#' @method plot wccv
#' @author Haotian Xu, Stephane Guerrier, and Justin Lee
#' @keywords internal
#' @export
plot.wccv = function(wccv.obj, theo.wccv = NULL){
  log.positive = sapply(wccv.obj$crosscovariance, function(x){ifelse(x < 0, NA, log(x))})
  log.negative = sapply(wccv.obj$crosscovariance, function(x){ifelse(x > 0, NA, log(-x))})
  log.low.positive = sapply(wccv.obj$ci_low, function(x){ifelse(x < 0, NA, log(x))})
  log.low.negative = sapply(wccv.obj$ci_low, function(x){ifelse(x > 0, NA, log(-x))})
  log.high.positive = sapply(wccv.obj$ci_high, function(x){ifelse(x < 0, NA, log(x))})
  log.high.negative = sapply(wccv.obj$ci_high, function(x){ifelse(x > 0, NA, log(-x))})
  
  y.positive.min = floor(min(log.positive, na.rm = T)) - 1
  y.positive.max = ceiling(max(log.positive, na.rm = T)) + 1
  
  y.negative.min = floor(min(-log.negative, na.rm = T)) - 1
  y.negative.max = ceiling(max(-log.negative, na.rm = T)) + 1
  
  y.min = min(y.positive.min, -y.negative.max)
  y.max = max(y.positive.max, -y.negative.min)
  
  layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
  par(mar = c(0, 3.1, 4.1, 2.1))
  plot(x = a$scales, y = log.positive, log = "x", type = 'p', xaxt = 'n', ylim = c(y.min, y.max), yaxt="n",
       main = 'Sample Wavelet Cross-Covariance', ylab = "")
  ticks <- seq(y.min+2, y.max, by=5)
  labels.positive <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(2, at=ticks, labels=labels.positive)
  lines(x = a$scales, y = log.low.positive)
  lines(x = a$scales, y = log.high.positive)
  
  if (is.null(theo.wccv) == F){
    log.theo.positive = sapply(theo.wccv, function(x){ifelse(x < 0, NA, log(x))})
    log.theo.negative = sapply(theo.wccv, function(x){ifelse(x > 0, NA, log(-x))})
    lines(x = a$scales, y = log.theo.positive, lty = 3)
  }
  
  par(mar = c(4.1, 3.1, 0, 2.1))
  plot(x = a$scales, y = -log.negative, log = "x", type = 'p', xaxt = 'n', ylim = c(-y.max, -y.min), yaxt="n", ylab = "", xlab = "Scales")
  ticks.rev <- -ticks
  labels.negative <- sapply(ticks.rev, function(i) as.expression(bquote(-10^ .(-i))))
  ticks.x <- seq(0, floor(log10(N))+1, by=1)
  labels.x <- sapply(ticks.x, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=10^ticks.x, labels=labels.x)
  axis(2, at=ticks.rev, labels=labels.negative)
  lines(x = a$scales, y = -log.low.negative)
  lines(x = a$scales, y = -log.high.negative)
  if (is.null(theo.wccv) == F){
    lines(x = a$scales, y = -log.theo.negative, lty = 3)
  }
}