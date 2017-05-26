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

#' Get ggplot2-like colors
#' @importFrom grDevices hcl
#' @param n number of colors.
#' @param alpha transparency.
#' @return list of colors.
ggplot_like_colors <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

#' @title Discrete Wavelet Transform
#' @description 
#' Calculation of the coefficients for the discrete wavelet transformation
#' @importFrom grDevices hcl
#' @export
#' @param x        A \code{vector} with dimensions N x 1. 
#' @param nlevels  A \code{integer} indicating the \eqn{J} levels of decomposition.
#' @param filter   A \code{string} indicating the filter name
#' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' The default \eqn{J} is determined by \eqn{floor\left(log_2 \left(length\left(x\right)\right)\right)}{floor(log2(length(x)))}
#' @author James Balamuta and Justin Lee
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' dwt(x)
dwt = function(x, nlevels = floor(log2(length(x))), filter = "haar") {
  
  if(is.vector(x) && length(x) %% 2^nlevels != 0){
    warning("The data has been truncated so that it is divisible by `nlevels` (e.g. 2^*)")
    x = x[1:2^nlevels]
  }else if(is.matrix(x) || is.data.frame(x)){
    if(ncol(x) != 1){
      stop("Only one column is allowed to be decomposed at a time.")
    }
    
    if(nrow(x) %% nlevels !=0){
      warning("The data has been truncated so that it is divisible by `nlevels` (e.g. 2^*)")
      idx = 1:2^nlevels
      x[idx,1] = x[idx,1]
    }
  }
  
  ret = dwt_cpp(x = x, filter_name = filter, nlevels)  # call to C++ version of dwt
  mostattributes(ret) = list(J=nrow(ret), filter = filter, class=c("dwt","list"))
  ret
}

#' @title Print Discrete Wavelet Transform
#' @description
#' Prints the results of the modwt list
#' @method print dwt
#' @export
#' @param x A \code{dwt} object
#' @param ... further arguments passed to or from other methods.
#' @return Prints the dwt decomposition
#' @author James Balamuta
#' @keywords internal
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' dwt(x)
print.dwt=function(x, ...){
  NextMethod("print")
}

#' @title Summary Discrete Wavelet Transform
#' @description 
#' Unlists DWT object and places it in matrix form
#' @method summary dwt
#' @import graphics 
#' @export
#' @keywords internal
#' @param object A \code{dwt} object
#' @param ... additional arguments affecting the summary produced.
#' @return Prints the dwt matrix decomposition
#' @author James Balamuta
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' summary(dwt(x))
summary.dwt=function(object, ...) {
  cat("Results of DWT using",attr(object,"filter"),"filter with",attr(object, "J")-1,"levels\n")
  y = as.list(object)
  j = length(y)
  for( i in 1:j ) {
    cat("Level",i,"Wavelet Coefficients")
    y[[i]]
    cat("\n")
  }
}

#' @title Plot Discrete Wavelet Transform
#' @description
#' Plots results of the dwt list in which additional parameters can be specified
#' @method plot dwt
#' @export
#' @param x A \code{dwt} object.
#' @param index A \code{vector} containing the indices to scales to be included in 
#' the graph. By default \code{index = 1:(min(c(J,4)))}, where \code{J} denotes the 
#' number of scales in \code{y}.
#' @param couleur A \code{vector} of colors of the same size as \code{index} used
#' for the different scales depicted in the graph. If \code{couleur} contains a single 
#' value the the same color will be used for all scales.
#' @param ... additional arguments affecting the plot produced.
#' @author Justin Lee and Stephane Guerrier
#' @keywords internal
#' @examples 
#' # Simulate a Gaussian white noise
#' n = 10^3
#' Xt = rnorm(n)
#' 
#' # dwt
#' Yt = dwt(Xt)
#' 
#' # Graph examples
#' plot(Yt)
#' plot(Yt, index = c(1,4,5,6,8,2))
#' plot(Yt, index = c(1,4,5,6), couleur = "blue")
#' plot(Yt, index = c(1,4,5,6), couleur = rep(c("blue","yellow"),2))
plot.dwt = function(x, index = NULL, couleur = NULL, ...){
  J <- attr(x,"J")
  
  if (is.null(index)){
    index <- 1:(min(c(4,J)))
  }else{
    if (max(index) > J || min(index) < 1){
      stop("Incorrect index specified")
    }
  }
  
  nb_plot <- length(index)
  
  if (is.null(couleur)){
    couleur <- ggplot_like_colors(nb_plot)
  }else{
    if (length(couleur) == 1 || length(couleur) != nb_plot){
      couleur <- rep(couleur[1],nb_plot)
    }
  }
  
  par(mfrow = c(nb_plot,1), mar = c(0,3,0,0), oma = c(5,2,1,1))
  x_range <- length(x[[1]])
  for (i in seq_len(nb_plot)){
    current_time_series <- x[[index[i]]]
    
    plot(NA, xlim = c(1,x_range), ylim = range(current_time_series), 
         bty = "n", axes = FALSE)
    box(col = "lightgrey")
    grid()
    axis(2)
    mtext(paste("Scale ",index[i], sep = ""), side = 2, line = 3, cex = 0.8)
    lines(unlist(current_time_series), col = couleur[i])
    
    # Add bottom axis
    if (i == nb_plot){
      axis(1)
      mtext("Time", side = 1, line = 3, cex = 0.8) 
    }
  }
}
