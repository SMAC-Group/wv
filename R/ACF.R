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


#' @title Auto-Covariance and Correlation Functions
#' @description The ACF function computes the estimated
#' autocovariance or autocorrelation for both univariate and multivariate cases.
#' @author Yunxiang Zhang
#' @param x      A \code{matrix} with dimensions \eqn{N \times S}{N x S} or N observations and S processes
#' @param lagmax A \code{integer} indicating the max lag.
#' @param cor    A \code{bool} indicating whether the correlation 
#' (\code{TRUE}) or covariance (\code{FALSE}) should be computed.
#' @param demean A \code{bool} indicating whether the data should be detrended
#'  (\code{TRUE}) or not (\code{FALSE})
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @details 
#' \code{lagmax} default is \eqn{10*log10(N/m)} where \eqn{N} is the number of
#' observations and \eqn{m} is the number of series being compared. If 
#' \code{lagmax} supplied is greater than the number of observations, then one
#' less than the total will be taken.
#' @examples 
#' # Get Autocorrelation
#' m = ACF(datasets::AirPassengers)
#' 
#' # Get Autocovariance and do not remove trend from signal
#' m = ACF(datasets::AirPassengers, cor = FALSE, demean = FALSE)





ACF = function(x, lagmax = 0, cor = TRUE, demean = TRUE){
  
  # Change the data to matrix form
  if(is.ts(x) || is.atomic(x)){
    x2 = data.matrix(x)        
  }
  
  #Get the acf value of the data
  acfe = acf(x, lagmax = lagmax , cor = cor, demean = demean, plot = FALSE)
  acfe = acfe$acf
  
  
  # Get the data name 
  varName = deparse(substitute(RW))
  
  
  # Adjust the name for data 
  dimnames(acfe)  = list(seq_len(nrow(acfe))-1, "ACF", varName)
  
  acfe = structure(acfe, n = nrow(x2), class = c("ACF","array"))
  acfe
  
}




#' @title Auto-Covariance and Correlation Functions
#' @description The acf function computes the estimated
#' autocovariance or autocorrelation for both univariate and multivariate cases.
#' @author Yunxiang Zhang
#' @param x,object  An \code{"ACF"} object from \code{\link{ACF}}.
#' @param show.ci   A \code{bool} indicating whether to show confidence region
#' @param ci        A \code{double} containing the 1-alpha level. Default is 0.95
#' @param ...       Additional parameters
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @rdname plot.ACF
#' @export
#' @examples 
#' # Calculate the Autocorrelation
#' m = ACF(datasets::AirPassengers)
#' 
#' # Plot with 95% CI
#' plot(m) 
#' 
#' # Plot with 90% CI
#' plot(m, ci = 0.90) 
#' 
#' # Plot without 95% CI
#' plot(m, show.ci = FALSE)


plot.ACF = function(x, show.ci = TRUE, ci = 0.95, ...){
  
  # Quiet the warnings...
  Lag = xmin = xmax = ymin = ymax = NULL 
  
  # Wide to long array transform
  x2 = as.data.frame.table(object, responseName = "ACF")
  
  colnames(x2) = c("Lag", "Signal X", "Signal Y", "ACF")
  
  # Remove character cast
  x2$Lag = as.numeric(x2$Lag)
  
  
  # Range
  x_range = range(x2$Lag)
  y_range = range(0:1)
  
  # Axie labels
  xlab = "Lags"
  ylab = "ACF"
  x_ticks = seq(x_range[1], x_range[2], by = 1)
  y_ticks = seq(y_range[1], y_range[2], by = 0.05)
  
  
  # Main plot
  plot(NA, xlim = x_range, ylim = y_range, main = paste0("The ACF plot of ",as.character((x2$`Signal Y`)[1])), xlab = xlab, ylab = ylab)
  
  
  # Add Grid
  abline(v = x_ticks, lty = 1, col = "grey95")
  abline(h = y_ticks, lty = 1, col = "grey95")
  
  
  # Plot CI 
  if(show.ci){
    
    clim0 = qnorm( (1 + ci)/2 ) / sqrt(attr(object,'n'))
    
    rect(xleft = -2, ybottom = -clim0, xright = 2*x_range[2], ytop = clim0, col = rgb(0, 1, 1, 0.3), lwd = 0)
    
  }
  
  abline(h = 0, lty = 1, lwd = 2)
  
  
  # Plot ACF
  segments(x0 = x_ticks, y0 = rep(0, x_range[2]), x1 = x_ticks, y1 = x2$ACF, lty = 1, lwd = 1)
  
  
}












