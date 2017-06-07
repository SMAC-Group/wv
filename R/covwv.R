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
  
  mostattributes(ret) = list(filter = filter, J = nlevels, N = length(x), class=c("wccv_pair","list","matrix"))
  
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
  
  # mostattributes(mat) = list(filter = filter, class=c("wccv","matrix","list"))
  
  return(mat)
}

#' @title Plot Cross Covariance Pair
#' @description
#' Plots results of the a wccv_pair list in which additional parameters can be specified
#' @author Justin Lee, Haotian Xu, and Stephane Guerrier
#' @method plot wvar
#' @keywords internal
#' @export
plot.wccv_pair = function(x, theo.wccv = NULL, main = NULL, col_wccv = NULL, col_ci = NULL, ...){
  J = attr(x,"J")
  N = attr(x, "N")
  scales = scales_cpp(J)
  
  x = x[[1]] # simplify 

  # Include all CI values 
  combCI = c(x[,3], x[,4])
  abscombCI = abs(combCI)
  
  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", tau, sep =""))
    }else{
      xlab = bquote(paste("Scale ", "", tau, " [", .(units), "]", sep = ""))
    }
  }
  
  if (is.null(ylab)){
    if(is.null(units)){
      ylab = expression(paste("Wavelet Cross Covariance ", "", (nu^2), "", sep = ""))
    }else{
      ylab = bquote(paste("Wavelet Cross Covariance ", "", (nu^2), " [", .(units)^2, "]", sep = ""))
    }
  }
  
  
  # Line and CI colors
  if(is.null(col_wv)){
    col_wccv = "darkblue"
  }
  
  if(is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }
  
  if(is.null(main)){
    main = "Sample Wavelet Cross-Covariance"
  }
  
  plot(NA, log = "x", xlim = c(min(scales), max(scales)), ylim = c(min(x[,3]), max(x[,4])), xaxt = "n", yaxt = "n",
       main = main, xlab = xlab, ylab = ylab)
  polygon(c(scales, rev(scales)), c(x[,3], rev(x[,4])),
          border = NA, col = col_ci)
  lines(x = scales, y = x[,1], type = "l", col = col_wccv, pch = 16, cex = 1.25)
  lines(x = scales, y = x[,1], type = "p", col = col_wccv, pch = 16)
  
  # set ticks and labels 
  ticks = seq(min(combCI), max(combCI), length = 5)
  tick.min = floor(min(log2(abscombCI))) # why not 10? 
  tick.max = floor(max(log2(abscombCI)))
  upper.ticks = seq(tick.min, tick.max, length = 2)
  lower.ticks = rev(upper.ticks)
  upper.labels = sapply(upper.ticks, function(i) as.expression(bquote(10^ .(i))))
  lower.labels = sapply(lower.ticks, function(i) as.expression(bquote(-10^ .(i))))
  labels = c(lower.labels, 0, upper.labels)
  axis(2, at=ticks, labels=labels)

  ticks.x <- seq(1, 15)
  labels.x <- sapply(ticks.x, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=2^ticks.x, labels=labels.x) # why 2 instead of 10?
    
  # not sure what this is for 
  if (is.null(theo.wccv) == F){
    log.theo.positive = sapply(theo.wccv, function(x){ifelse(x < 0, NA, log(x))})
    log.theo.negative = sapply(theo.wccv, function(x){ifelse(x > 0, NA, log(-x))})
    lines(x = scales, y = log.theo.positive, lty = 3)
  }
  
  # not sure what this is for 
  if (is.null(theo.wccv) == F){
    lines(x = scales, y = -log.theo.negative, lty = 3)
  }
  
  
}