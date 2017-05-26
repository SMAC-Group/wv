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
#' # Robust
#' wvar(x, robust = TRUE, eff=0.3)
#' # 90% confidence interval
#' wvar(x, alpha = 0.10)
#' 
#' # IMU Object
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyros = 1:3, accels = 4:6, freq = 100)
#' df = wvar.imu(test)
#' }
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

#' Compare Wavelet Variances
#' 
#' Compare the estimates given by the classical and robust methods of 
#' calculating the wavelet variance.
#' @param ... Any number of \code{wvar} or \code{wvar.imu} objects.
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE)
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param auto.label.wvar A \code{boolean} that indicates whether legend label should indicate the objects are robust or classical.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of confidence interval.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of point.
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of point. 
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param axis.x.label A \code{string} that indicates the label on x axis
#' @param axis.y.label A \code{string} that indicates the label on y axis
#' @param units A two-element vector indicating the units of gyroscope and accelerometer sensor. Set it to \code{NULL} if units are not needed. 
#' @param facet.label.size An \code{integer} that indicates the size of facet label
#' @param facet.label.background A \code{string} that indicates the background color of the facet label
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend
#' @param legend.title.size An \code{integer} that indicates the size of title on legend
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param nrow An \code{integer} that indicates number of rows
#' @author JJB, Wenchao
#' @note
#' If you meet the error "polygon edge not found", RStudio is complaining that you don't have enough 
#' space to plot the graph. You can adjust the graphics window, or open an external window. The 
#' function \code{\link{external_graphs}} can be used and it works for all operating systems.
#' 
#' When \code{wvar.imu} objects are supplied, some parameters, e.g. \code{split} and \code{nrow}, are 
#' invalid, since the graph is plot seperately and put in 2 rows by default.
#' 
#' @examples
#' \dontrun{
#' ## Case1: Supplied objects are \code{wvar}:
#' data1 = gen_gts(1000, AR1(phi = .32, sigma2=.01))
#' data2 = gen_gts(2000, ARMA(ar=c(.8,.1), ma=c(.3), sigma2=1))
#' data3 = gen_gts(4000, AR1(phi = .32, sigma2=1))
#' 
#' wv1 = wvar(data1, robust = TRUE)
#' wv2 = wvar(data2)
#' wv3 = wvar(data3)
#' 
#' compare_wvar(wv1, wv2)
#' compare_wvar(wv1, wv2, CI = FALSE)
#' compare_wvar(wv1, wv2, split = FALSE)
#' compare_wvar(wv1, wv2, wv3, split = FALSE)
#' 
#' # Change default setting
#' color = c('green','red','blue')
#' label = c('1','2','3')
#' compare_wvar(wv1, wv2, wv3, 
#'              line.color = color, CI.color = color,
#'               legend.label = label)
#' compare_wvar(wv1, wv2, wv3, 
#'              line.color = color, CI.color = color, 
#'              legend.label = label, split = FALSE)
#' 
#' ## Case2: Supplied objects are \code{wvar.imu}:
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test1 = imu(imu6, gyros = 1:3, accels = 4:6, axis = c('X', 'Y', 'Z'), freq = 100)
#' wv1 = wvar(test1)
#' 
#' test2 = imu(imu6, gyros = 1, accels = 3:4, axis = c('X','X','Y'), freq = 100)
#' wv2 = wvar(test2, robust = T)
#' 
#' compare_wvar(wv1, wv2)
#' compare_wvar(wv1, wv2, auto.label.wvar = F, legend.label = c('data1', 'data2'))
#' }
compare_wvar = function(..., background = 'white', split = TRUE, CI = TRUE, auto.label.wvar = T, transparence = 0.1, line.color = NULL, 
                        CI.color = NULL, line.type = NULL,  point.size = NULL, point.shape = NULL,
                        title = "Haar Wavelet Variance Representation", title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = expression(paste("Scale ", tau)),
                        axis.y.label = expression(paste("Wavelet Variance ", nu)),
                        units = c(bquote(rad^2/s^2), bquote(m^2/s^4)),
                        facet.label.size = 13, facet.label.background = "#003C7D33",
                        legend.label = NULL,
                        legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                        legend.text.size = 13, nrow = 1 ){
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  obj_list = list(...)
  object.names = as.character(substitute(...()))
  numObj = length(obj_list)
  
  if (numObj == 0){
    stop('At least one wvar object should be given')
    
  }else if (numObj == 1){
    ## just plot
    plot(...)
    
  }else{
    
    # Case1: Supplied objects are 'wvar.imu'
    is.wvar.imu = sapply(obj_list, FUN = is, class2 = 'wvar.imu')
    
    # Case2: Supplied objects are 'wvar'
    is.wvar = sapply(obj_list, FUN = is, class2 = 'wvar')
    
    if(!all(is.wvar == T) && !all(is.wvar.imu == T)  ){
      stop("Supplied objects must be either 'wvar' or 'wvar.imu' objects.")
    }
    
    # check parameter
    if(all(is.wvar.imu == T)){
      # Case1: Supplied objects are 'wvar.imu'
      if(CI){
        params = c('line.color', 'line.type', 'CI.color', 'point.size', 'point.shape', 'legend.label')
        requireLength = c(3*numObj, 3*numObj, numObj, 3*numObj, 3*numObj, numObj)
        default = list(NULL, NULL,  NULL, NULL, NULL, NULL)
        nullIsFine = c(rep(T,6))
        
      }else{
        params = c('line.color', 'line.type', 'point.size', 'point.shape', 'legend.label')
        requireLength = c(numObj, numObj, numObj, numObj, numObj)
        default = list(NULL, NULL, NULL, NULL, NULL)
        nullIsFine = c(rep(T,5))
      }
    }else{
      # Case2: Supplied objects are 'wvar'
      # line.type is set here
      params = c('line.type', 'line.color', 'CI.color', 'point.size', 'point.shape', 'legend.label')
      requireLength = c(2, numObj, numObj, numObj, numObj, numObj)
      default = list(c('solid','dotted'), NULL,  NULL, rep(5, numObj), rep(20, numObj), NULL)
      nullIsFine = c(rep(T,6))
    }

    # legend.label
    if(is.null(legend.label)){
      legend.label = object.names
    }
    
    if(auto.label.wvar){
      
      rob = sapply( obj_list, FUN = function(x){
        if(is(x, 'wvar.imu')){
          control = x$dataobj[[1]]$robust
        }else{control = x$robust }
         
        if(control){'(Robust)'}else{'(Classical)'} }
      )
      
      legend.label = paste(legend.label, rob)
    }
    # make sure legend.label does not have duplicates
    legend.label = addSpaceIfDuplicate(legend.label)
    
    
    if(all(is.wvar.imu == T)){
      
      return( compare_wvar.imu(obj.list = obj_list,
                       background = background, CI = CI, auto.label.wvar = auto.label.wvar, transparence = transparence, 
                       line.color = line.color, CI.color = CI.color, line.type = line.type, point.size = point.size, point.shape = point.shape,
                       title = title, title.size = title.size, 
                       axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                       axis.x.label = axis.x.label,
                       axis.y.label = axis.y.label,
                       units = units,
                       facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                       legend.label = legend.label,
                       legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                       legend.text.size = legend.text.size ) )
      
    }
    
    total.len = 0
    each.len = numeric(numObj)
    for (i in 1:numObj){
      each.len[i] = length(obj_list[[i]]$variance)
      total.len = total.len + each.len[i]
    }
    #Initialize empty data frame with right number of rows
    obj = data.frame(WV = numeric(total.len),
                     scales = numeric(total.len),
                     low = numeric(total.len),
                     high = numeric(total.len),
                     dataset = 'XYZ', stringsAsFactors=FALSE)
    
    #put data into data frame
    t = 1
    for (i in 1:numObj){
      d = each.len[i]
      obj[t:(t+d-1),] = data.frame(WV = obj_list[[i]]$variance,
                                   scales = obj_list[[i]]$scales,
                                   low = obj_list[[i]]$ci_low,
                                   high = obj_list[[i]]$ci_high,
                                   dataset = legend.label[i], stringsAsFactors=FALSE)
      t = t +d
    }
    
    if (numObj == 2 ){
      if(is.null(line.color)){
        line.color = c("#003C7D","#F47F24")
      }
      if(is.null(CI.color)){
        CI.color = c("#003C7D","#F47F24")
      }
    }
    
    autoplot.wvarComp(obj, split = split, CI = CI, background = background, transparence = transparence, line.color =line.color, 
                        CI.color = CI.color, line.type = line.type,  point.size = point.size, point.shape = point.shape,
                        title = title, title.size= title.size, 
                        axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                        axis.x.label = axis.x.label,
                        axis.y.label = axis.y.label,
                        facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                        legend.label = legend.label,
                        legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                        legend.text.size = legend.text.size,
                        nrow = nrow)
  }
  
}
