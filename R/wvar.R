# Copyright (C) 2014 - 2017  James Balamuta, Stephane Guerrier, Roberto Molinari
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

#' Wavelet Variance
#' 
#' Calculates the (MODWT) wavelet variance
#' @param x         A \code{vector} with dimensions N x 1, or a \code{lts} object, or a \code{gts} object, or a \code{imu} object. 
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
#' @author JJB
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

#' Create a \code{wvar} object
#' 
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

#' Print Wavelet Variances
#' 
#' Displays the summary table of wavelet variance.
#' @method print wvar
#' @export
#' @keywords internal
#' @param x A \code{wvar} object.
#' @param ... further arguments passed to or from other methods.
#' @author JJB
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

#' Summary of Wavelet Variances
#' 
#' Displays the summary table of wavelet variance in addition to CI values and
#' supplied efficiency.
#' @method summary wvar
#' @export
#' @keywords internal
#' @param object A \code{wvar} object.
#' @return Summary table and other properties of the object.
#' @param ... additional arguments affecting the summary produced.
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = wvar(x)
#' summary( out )
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
