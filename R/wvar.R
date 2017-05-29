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

#' @title Wavelet Variance
#' 
#' @description
#' Calculates the (MO)DWT wavelet variance
#' @param x         A \code{vector} with dimensions N x 1. 
#' @param decomp    A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition.
#' @param filter    A \code{string} that specifies what wavelet filter to use. 
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @param robust    A \code{boolean} that triggers the use of the robust estimate.
#' @param eff       A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param alpha     A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
#' @param freq      A \code{numeric} that provides the rate of samples.
#' @param from.unit A \code{string} indicating the unit which the data is converted from.
#' @param to.unit   A \code{string} indicating the unit which the data is converted to.
#' @param ... Further arguments passed to or from other methods.
#' @return A \code{list} with the structure:
#' \describe{
#'   \item{"variance"}{Wavelet Variance}
#'   \item{"ci_low"}{Lower CI}
#'   \item{"ci_high"}{Upper CI}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"alpha"}{p value used for CI}
#'   \item{"unit"}{String representation of the unit}
#' }
#' @details 
#' If \code{nlevels} is not specified, it is set to \eqn{\left\lfloor {{{\log }_2}\left( {length\left( x \right)} \right)} \right\rfloor}{floor(log2(length(x)))}
#' @author James Balamuta and Justin Lee
#' @rdname wvar
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' # Default
#' wvar(x)
#' 
#' # Robust
#' wvar(x, robust = TRUE, eff=0.3)
#' 
#' # Classical
#' wvar(x, robust = FALSE, eff=0.3)
#' 
#' # 90% Confidence Interval 
#' wvar(x, alpha = 0.10)
#' @export
wvar = function(x, ...) {
  UseMethod("wvar")
}

#' @rdname wvar
#' @export
wvar.lts = function(x, decomp = "modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  warning('`lts` object is detected. This function can only operate on the combined process.')
  freq = attr(x, 'freq')
  unit = attr(x, 'unit')
  x = x[,ncol(x)]
  
  wvar.default(x, decomp, filter, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.gts = function(x, decomp="modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  freq = attr(x, 'freq')
  unit = attr(x, 'unit')
  x = x[,1]
  
  wvar.default(x, decomp, filter, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.ts = function(x, decomp="modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  freq = attr(x, 'tsp')[3]
  unit = NULL
  
  wvar.default(x, decomp, filter, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.default = function(x, decomp = "modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, freq = 1, from.unit = NULL, to.unit = NULL, ...){
  if(is.null(x)){
    stop("`x` must contain a value")
  }else if((is.data.frame(x) || is.matrix(x))){
    if(ncol(x) > 1) stop("There must be only one column of data supplied.")
  }
  
  if(is.null(nlevels)){
    nlevels = floor(log2(length(x)))
  }

  # check freq
  if(!is(freq,"numeric") || length(freq) != 1){ stop("'freq' must be one numeric number.") }
  if(freq <= 0) { stop("'freq' must be larger than 0.") }
  
  # check unit
  all.units = c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')
  if( (!is.null(from.unit) && !from.unit %in% all.units) || (!is.null(to.unit) && !to.unit %in% all.units) ){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
  }
  
  if(robust) {
    if(eff > 0.99) {
      stop("The efficiency specified is too close to the classical case. Use `robust = FALSE`")
    }
  }
  
  obj =  .Call('wv_modwt_wvar_cpp', PACKAGE = 'wv',
               signal=x, nlevels=nlevels, robust=robust, eff=eff, alpha=alpha, 
               ci_type="eta3", strWavelet=filter, decomp = decomp)
  
  # nlevels may be changed during modwt
  nlevels = nrow(obj)
  
  scales = .Call('wv_scales_cpp', PACKAGE = 'wv', nlevels)/freq
  
  # NO unit conversion
  if( is.null(from.unit) && is.null(to.unit)==F ){
    warning("'from.unit' is NULL. Unit conversion was not done.")
  }
  
  # unit conversion
  if (!is.null(from.unit)){
    if (!is.null(to.unit)){
      convert.obj = unitConversion(scales, from.unit = from.unit, to.unit = to.unit)
      
      if (convert.obj$converted) {
        # YES unit conversion
        scales = convert.obj$x
        message(paste0('Unit of object is converted from ', from.unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  if(!is.null(from.unit) && !is.null(to.unit)){ 
    unit = to.unit
  }else{
    unit = from.unit}
  
  create_wvar(obj, decomp, filter, robust, eff, alpha, scales, unit)
}

#' @title Create a \code{wvar} object
#' 
#' @description
#' Structures elements into a \code{wvar} object
#' @param obj    A \code{matrix} with dimensions N x 3, that contains the wavelet variance, low ci, hi ci.
#' @param decomp A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition
#' @param filter A \code{string} that specifies the type of wavelet filter used in the decomposition
#' @param robust A \code{boolean} that triggers the use of the robust estimate.
#' @param eff    A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param alpha  A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
#' @param scales A \code{vec} that contains the amount of decomposition done at each level.
#' @param unit   A \code{string} that contains the unit expression of the frequency.
#' @return A \code{list} with the structure:
#' \describe{
#'   \item{"variance"}{Wavelet Variance}
#'   \item{"ci_low"}{Lower CI}
#'   \item{"ci_high"}{Upper CI}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"alpha"}{p value used for CI}
#'   \item{"unit"}{String representation of the unit}
#' }
#' @keywords internal
create_wvar = function(obj, decomp, filter, robust, eff, alpha, scales, unit){
  structure(list(variance = obj[,1],
                       ci_low = obj[,2], 
                       ci_high = obj[,3], 
                       robust = robust, 
                       eff = eff,
                       alpha = alpha,
                       scales = scales,
                       decomp = decomp,
                       unit = unit,
                       filter = filter), class = "wvar")
}

#' @title Print Wavelet Variances
#' 
#' @description
#' Displays the summary table of wavelet variance.
#' @author James Balamuta 
#' @method print wvar
#' @export
#' @keywords internal
#' @param x A \code{wvar} object.
#' @param ... further arguments passed to or from other methods.
#' @return Summary table
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = wvar(x)
#' print( out )
print.wvar = function(x, ...){
  mat = matrix(unlist(x[1:3]),ncol=3,byrow=F)
  colnames(mat) = c("Variance", "Low CI", "High CI")
  rownames(mat) = x$scales
  print(mat)
}

#' @title Summary of Wavelet Variances
#' 
#' @description 
#' Displays the summary table of wavelet variance accounting for CI values and supplied efficiency.
#' @method summary wvar
#' @export
#' @keywords internal
#' @param object A \code{wvar} object.
#' @return Summary table and other properties of the object.
#' @param ... additional arguments affecting the summary produced.
#' @author James Balamuta
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' ret = wvar(x)
#' summary(ret)
summary.wvar = function(object, ...){
  name = if(object$robust){
    "robust" 
  }else{
    "classical"
  }
  cat("Results of the wavelet variance calculation using the ",name, " method.\n",sep="")
  if(object$robust){
    cat("Robust was created using efficiency=",object$eff,"\n",sep="")
  }
  
  cat("The confidence interval was generated using (1-",object$alpha,")*100 \n",sep="")
  
  print(object)
}

#' @title Plot Wavelet Variances
#' 
#' @description 
#' Displays a plot of wavelet variance accounting for CI values and supplied efficiency.
#' @method plot wvar
#' @export
#' @keywords internal
#' @param \code{x} A \code{wvar} object.
#' @param \code{xlab} A title for the x axis.
#' @param \code{ylab} A title for the y axis.
#' @param \code{main} An overall title for the plot.
#' @param \code{col_wv} Color of the wavelet variance line.
#' @param \code{col_ci} Color of the confidence interval shade.
#' @param \code{nb_ticks_x} Number of ticks for the x-axis.
#' @param \code{nb_ticks_y} Number of ticks for the y-axis.
#' @param \code{legend_position} Position of the legend.
#' @param \code{...} Additional arguments affecting the plot.
#' @return Plot of wavelet variance and confidence interval for each scale.
#' @author Stephane Gurrier, Nathanael Claussen, and Justin Lee
plot.wvar = function(x, units = NULL, xlab = NULL, ylab = NULL, main = NULL, 
                     col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                     legend_position = NULL, ...){
  # Control margins (bottom, left, top, right)
  par(oma = c(1,0,1,0), mar=c(4,5,2,1))
  
  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", "(", tau, ")", sep =""))
    }else{
      xlab = bquote(paste("Scale ", "(", tau, ") [", .(units), "]", sep = ""))
    }
  }
  
  if (is.null(ylab)){
    if(is.null(units)){
      ylab = expression(paste("Wavelet Variance ", "(", hat(nu), ")", sep = ""))
    }else{
      ylab = bquote(paste("Wavelet Variance ", "(", hat(nu), ") [", .(units)^2, "]", sep = ""))
    }
  }
  
  # Main Title
  if (is.null(main)){
    main = "Haar Wavelet Variance Representation"
  }
  
  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = "darkblue"
  }
  
  if (is.null(col_ci)){
    col_ci <- hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }
  
  # Range
  x_range = range(x$scales)
  x_low = floor(log2(x_range[1]))
  x_high = ceiling(log2(x_range[2]))
  
  y_range = range(c(x$ci_low, x$ci_high))
  y_low = floor(log2(y_range[1]))
  y_high = ceiling(log2(y_range[2]))
  
  # Axes
  if (is.null(nb_ticks_x)){
    nb_ticks_x = 6
  }
  
  if (is.null(nb_ticks_y)){
    nb_ticks_y = 5
  }
  
  x_ticks = seq(x_low, x_high, by = 1)
  if (length(x_ticks) > nb_ticks_x){
    x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
  }
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))
  
  y_ticks <- seq(y_low, y_high, by = 1)
  if (length(y_ticks) > nb_ticks_y){
    y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
  }
  y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
  
  # Legend position
  if (is.null(legend_position)){
    if (which.min(abs(c(y_low, y_high) - log2(x$variance[1]))) == 1){
      legend_position = "topleft"
    }else{
      legend_position = "bottomleft"
    }
  }   
  
  # Main plot                     
  plot(NA, xlim = x_range, ylim = y_range, xlab = xlab, ylab = ylab, 
       log = "xy", xaxt = 'n', yaxt = 'n', bty = "n")
  abline(v = 2^x_ticks, lty = 1, col = "grey95")
  abline(h = 2^y_ticks, lty = 1, col = "grey95")
  axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  axis(2, at = 2^y_ticks, labels = y_labels, padj = -0.2)  
  polygon(c(x$scales, rev(x$scales)), c(x$ci_low, rev(x$ci_high)),
          border = NA, col = col_ci)
  
  
  polygon(c(x$scales, rev(x$scales)), c(x$ci_low, rev(x$ci_high)),
          border = NA, col = col_ci)
  
  CI_conf = 1 - x$alpha
  legend(legend_position,
         legend=c(expression(paste("Empirical WV ",hat(nu))), bquote(paste("CI(",hat(nu),", ",.(CI_conf),")"))),
         pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
  
  lines(x$scales, x$variance, type = "l", col = col_wv, pch = 16)
  lines(x$scales, x$variance, type = "p", col = col_wv, pch = 16, cex = 1.25)
  
  par(xpd = NA)
  win_dim = par("usr")
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = 10^c(win_dim[4], win_dim[4],
               win_dim[4] + 0.09*(win_dim[4] - win_dim[3]), 
               win_dim[4] + 0.09*(win_dim[4] - win_dim[3]))
  lines(x_vec[1:2], rep(10^win_dim[3],2), col = 1)
  lines(rep(x_vec[1], 2), 10^win_dim[3:4], col = 1)
  lines(rep(x_vec[2], 2), 10^win_dim[3:4], col = 1)
  polygon(x_vec, y_vec, col = "grey95")
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = 10^(win_dim[4] + 0.09/2*(win_dim[4] - win_dim[3])), main)
  par(xpd = FALSE)
}
