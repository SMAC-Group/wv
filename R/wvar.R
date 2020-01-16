#' @title Convert Unit of Time Series Data
#' @description Manipulate the units of time to different ones
#' @keywords internal
#' @param x          A \code{vector} containing the values on x-axis.
#' @param from.unit  A \code{string} indicating the unit which the data is converted from.
#' @param to.unit    A \code{string} indicating the unit which the data is converted to.
#' @details
#' The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' Make sure \code{from.unit} and \code{to.unit} are not \code{NULL} before it is passed to this function.
#' @return A \code{list} with the following structure:
#' \itemize{
#'  \item "x": Data
#'  \item "converted": A \code{boolean} indicating whether conversion is made
#' }
#' @export
#' @examples
#' x = seq(60, 3600, 60)
#' unitConversion(x, 'sec', 'min')
#' y = 1:10
#' unitConversion(y, 'hour', 'sec')
unitConversion = function(x, from.unit, to.unit){
  
  #ns, ms, second, min, hour, day, month, year
  unit = c(ns = 1, ms = 2,se = 3, mi = 4, ho = 5, da = 6, mo = 7, ye = 8)
  
  #assume 1 month = 30 days
  ratio = c(1E6, 1E3, 60, 60, 24, 30, 12)
  from.unit.1 = substr(from.unit, 1, 2)
  to.unit.1 = substr(to.unit, 1, 2)
  
  #check unit:
  no.convert = F
  if(from.unit.1 == to.unit.1){no.convert = T}
  if(is.na(unit[from.unit.1]) ) {
    message = paste('No such unit: ', from.unit, '. Supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year". Conversion is terminated.', sep = '')
    warning(message); no.convert = T}
  if(is.na(unit[to.unit.1]) ) {
    message = paste('No such unit: ', to.unit, '. Supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year". Conversion is terminated.', sep = '')
    warning(message); no.convert = T}
  
  if(!no.convert){
    #print out warning when day is convert to month, or month is converted to day.
    conversionRange = unit[from.unit.1] : unit[to.unit.1]
    if(6 %in% conversionRange && 7 %in% conversionRange){
      warning('Unit conversion might be wrong because this function simply assumes 1 month = 30 days.')
    }
  }
  
  if(!no.convert){
    if(unit[from.unit.1] > unit[to.unit.1]){
      temp = ratio[unit[to.unit.1]: (unit[from.unit.1]-1)]
      multiplier = prod(temp)
      x = x*multiplier
    }else{
      temp = ratio[unit[from.unit.1]: (unit[to.unit.1]-1) ]
      multiplier = prod(temp)
      x = x/multiplier
    }
  }
  obj = list(x = x, converted = !no.convert)  
  return(obj)
}

#' @title Wavelet Variance
#' @description Calculates the (MO)DWT wavelet variance
#' @param x         A \code{vector} with dimensions N x 1.
#' @param decomp    A \code{string} that indicates whether to use a "dwt" or "modwt" decomposition.
#' @param filter    A \code{string} that specifies which wavelet filter to use.
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @param alpha     A \code{double} that specifies the significance level which in turn specifies the \eqn{1-\alpha} confidence level.
#' @param robust    A \code{boolean} that triggers the use of the robust estimate.
#' @param eff       A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param freq      A \code{numeric} that provides the rate of samples.
#' @param from.unit A \code{string} indicating the unit from which the data is converted.
#' @param to.unit   A \code{string} indicating the unit to which the data is converted.
#' @param ...       Further arguments passed to or from other methods.
#' @return A \code{list} with the structure:
#' \itemize{
#'   \item "variance": Wavelet Variance
#'   \item "ci_low": Lower CI
#'   \item "ci_high": Upper CI
#'   \item "robust": Robust active 
#'   \item "eff": Efficiency level for Robust calculation
#'   \item "alpha": p value used for CI
#'   \item "unit": String representation of the unit
#' }
#' @details 
#' The default value of \code{nlevels} will be set to \eqn{\left\lfloor {{{\log }_2}\left( {length\left( x \right)} \right)} \right\rfloor}{floor(log2(length(x)))}, unless otherwise specified.
#' @author James Balamuta, Justin Lee and Stephane Guerrier
#' @rdname wvar
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' 
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
  if (sum(class(x) %in% "gts") == 1){
    x = as.numeric(x)
  }
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
wvar.imu = function(x, decomp="modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  # Retrive sensor name
  if (!is.null(attr(x, "stype"))){
    sensor_name = attr(x, "stype")
  }else{
    warning("Unknown sensor name. IMU object is missing some information.")
    sensor_name = NULL
  }
  
  # Retrive freq
  if (!is.null(attr(x, "freq"))){
    freq = attr(x, "freq")
  }else{
    warning("Unknown frequency. IMU object is missing some information. Freq is set to 1 by default.")
    freq = 1
  }
  
  # Retrive sample size
  if (!is.null(attr(x, "dim"))){
    n = attr(x, "dim")[1]
  }else{
    warning("Unknown sample size. IMU object is missing some information.")
    n = NULL
  }
  
  # Retrive col names
  if (!is.null(attr(x, "dimnames")[[2]])){
    col_names = attr(x, "dimnames")[[2]]
  }else{
    stop("Unknown colunms names. IMU object is missing some information.")
    col_names = NULL
  }
  
  # Retrive sensor
  if (!is.null(attr(x, "sensor"))){
    sensor = attr(x, "sensor")
  }else{
    warning("Unknown sensor. IMU object is missing some information.")
    sensor = NULL
  }
  
  # Retrive axis
  if (!is.null(attr(x, "axis"))){
    ax = attr(x, "axis")
  }else{
    warning("Unknown axes. IMU object is missing some information.")
    ax = NULL
  }
  
  # Compute wvar
  m = length(col_names)
  wvariance = list()
  for (i in 1:m){
    wvariance[[i]] = wvar.default(x[,i], decomp, filter, nlevels, alpha, robust, eff, freq = freq, to.unit = to.unit)
  }
  names(wvariance) = col_names
  out = list(sensor = sensor_name, freq = freq, n = n, type = sensor, axis = ax, wvar = wvariance)
  class(out) = "imu_wvar"
  invisible(out)
}


#' @rdname wvar
#' @export
#' @importFrom methods is
wvar.default = function(x, decomp = "modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, freq = 1, from.unit = NULL, to.unit = NULL, ...){
  if(is.null(x)){
    stop("`x` must contain a value")
  }else if((is.data.frame(x) || is.matrix(x))){
    if(ncol(x) > 1) stop("There must be only one column of data supplied.")
  }
  
  if(decomp == "modwt" && is.null(nlevels)){
    nlevels = floor(log2(length(x))-1)
  }else if(decomp == "dwt" && is.null(nlevels)){
    nlevels = floor(log2(length(x)))
  }
  
  # Check Freq
  if(!is(freq,"numeric") || length(freq) != 1){ stop("'freq' must be one numeric number.") }
  if(freq <= 0) { stop("'freq' must be larger than 0.") }
  
  # Check Unit
  all.units = c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')
  if( (!is.null(from.unit) && !from.unit %in% all.units) || (!is.null(to.unit) && !to.unit %in% all.units) ){
    stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
  }
  
  if(robust) {
    if(eff > 0.99) {
      stop("The efficiency specified is too close to the classical case. Use `robust = FALSE`")
    }
  }
  
  obj = modwt_wvar_cpp(signal=x, nlevels=nlevels, robust=robust, eff=eff, alpha=alpha, 
              ci_type="eta3", strWavelet=filter, decomp = decomp)
  
  # nlevels may be changed during modwt
  nlevels = nrow(obj)
  
  scales = scales_cpp(nlevels)/freq
  
  # NO Unit Conversion
  if( is.null(from.unit) && is.null(to.unit)==F ){
    warning("'from.unit' is NULL. Unit conversion was not done.")
  }
  
  # Unit Conversion
  if (!is.null(from.unit)){
    if (!is.null(to.unit)){
      convert.obj = unitConversion(scales, from.unit = from.unit, to.unit = to.unit)
      
      if (convert.obj$converted) {
        # YES Unit Conversion
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
#' @description Structures elements into a \code{wvar} object
#' @param obj    A \code{matrix} with dimensions N x 3 that contains Wavelet Variance, Lower CI, and Upper CI.
#' @param decomp A \code{string} that indicates whether to use a "dwt" or "modwt" decomposition.
#' @param filter A \code{string} that specifies the type of wavelet filter used in the decomposition.
#' @param robust A \code{boolean} that triggers the use of the robust estimate.
#' @param eff    A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param alpha  A \code{double} that specifies the significance level which in turn specifies the \eqn{1-\alpha} confidence level.
#' @param scales A \code{vec} that contains the amount of decomposition performed at each level.
#' @param unit   A \code{string} that indicates the unit expression of the frequency.
#' @return A \code{list} with the structure:
#' \itemize{
#'   \item "variance": Wavelet variance
#'   \item "ci_low": Lower CI
#'   \item "ci_high": Upper CI
#'   \item "robust": Robust active
#'   \item "eff": Efficiency level for robust calculation
#'   \item "alpha": p value used for CI
#'   \item "unit": String representation of the unit
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
#' @description Displays the summary table of wavelet variance.
#' @author James Balamuta 
#' @method print wvar
#' @export
#' @keywords internal
#' @param x     A \code{wvar} object.
#' @param ...   Further arguments passed to or from other methods.
#' @return      Summary table
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
#' @description Displays the summary table of wavelet variance accounting for CI values and supplied efficiency.
#' @method summary wvar
#' @export
#' @keywords internal
#' @param object  A \code{wvar} object.
#' @param ...     Additional arguments affecting the summary produced.
#' @return        Summary table and other properties of the object.
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

#' @title Plot Wavelet Variance
#' @description Displays a plot of wavelet variance accounting for CI values and supplied efficiency.
#' @method plot wvar
#' @keywords internal
#' @param x                A \code{wvar} object.
#' @param units            A \code{string} that specifies the units of time plotted on the x axis.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_wv           A \code{string} that specifies the color of the wavelet variance line.
#' @param col_ci           A \code{string} that specifies the color of the confidence interval polygon.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position  A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param ci_wv            A \code{boolean} that determines whether a confidence interval polygon will be drawn.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of wavelet variance and confidence interval for each scale.
#' @author Stephane Guerrier, Nathanael Claussen, and Justin Lee
#' @export
#' @examples 
#' set.seed(999)
#' n = 10^4
#' Xt = rnorm(n)
#' wv = wvar(Xt)
#' plot(wv)
#' plot(wv, main = "Simulated white noise", xlab = "Scales")
#' plot(wv, units = "sec", legend_position = "topright")
#' plot(wv, col_wv = "darkred", col_ci = "pink")
plot.wvar = function(x, units = NULL, xlab = NULL, ylab = NULL, main = NULL, 
                     col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                     legend_position = NULL, ci_wv = NULL, point_cex = NULL, 
                     point_pch = NULL, ...){
  
  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", tau, sep =""))
    }else{
      xlab = bquote(paste("Scale ", tau, " [", .(units), "]", sep = " "))
    }
  }
  
  if (is.null(ylab)){
    ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
  }else{
    ylab = ylab
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
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }
  
  # Range
  x_range = range(x$scales)
  x_low = floor(log2(x_range[1]))
  x_high = ceiling(log2(x_range[2]))
  
  y_range = range(c(x$ci_low, x$ci_high))
  y_low = floor(log10(y_range[1]))
  y_high = ceiling(log10(y_range[2]))
  
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
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))
  
  y_ticks <- seq(y_low, y_high, by = 1)
  if (length(y_ticks) > nb_ticks_y){
    y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
  }
  y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
  
  # Legend Position
  if (is.null(legend_position)){
    if (which.min(abs(c(y_low, y_high) - log2(x$variance[1]))) == 1){
      legend_position = "topleft"
    }else{
      legend_position = "bottomleft"
    }
  }   
  
  # Main Plot                     
  plot(NA, xlim = x_range, ylim = y_range, xlab = xlab, ylab = ylab, 
       log = "xy", xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE)
  win_dim = par("usr")
  
  par(new = TRUE)
  plot(NA, xlim = x_range, ylim = 10^c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
       xlab = xlab, ylab = ylab, log = "xy", xaxt = 'n', yaxt = 'n', bty = "n")
  win_dim = par("usr")
  
  # Add Grid
  abline(v = 2^x_ticks, lty = 1, col = "grey95")
  abline(h = 10^y_ticks, lty = 1, col = "grey95")
  
  # Add Title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = 10^c(win_dim[4], win_dim[4],
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]), 
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = 10^(win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)
  
  # Add Axes and Box
  lines(x_vec[1:2], rep(10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)
  #y_ticks = y_ticks[(2^y_ticks) < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
  y_labels = y_labels[1:length(y_ticks)]
  box()
  axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2)  
  
  # CI for WV
  if (ci_wv == TRUE || is.null(ci_wv)){
    polygon(c(x$scales, rev(x$scales)), c(x$ci_low, rev(x$ci_high)),
            border = NA, col = col_ci)
  }
  
  # Add legend
  CI_conf = 1 - x$alpha
  
  if (x$robust == TRUE){
    wv_title_part1 = "Empirical Robust WV "
  }else{
    wv_title_part1 = "Empirical WV "
  }
  
  if (!is.na(legend_position)){
    if (legend_position == "topleft"){
      legend_position = 10^c(1.1*win_dim[1], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
      legend(x = legend_position[1], y = legend_position[2],
             legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))),
                        as.expression(bquote(paste("CI(",hat(nu)^2,", ",.(CI_conf),")")))),
             pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
    }else{
      if (legend_position == "topright"){
        legend_position = 10^c(0.7*win_dim[2], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
        legend(x = legend_position[1], y = legend_position[2],
               legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))), 
                          as.expression(bquote(paste("CI(",hat(nu)^2,", ",.(CI_conf),")")))),
               pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
      }else{
        legend(legend_position,
               legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))), 
                          as.expression(bquote(paste("CI(",hat(nu)^2,", ",.(CI_conf),")")))),
               pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
      }
    }
  }
  
  # Add WV
  lines(x$scales, x$variance, type = "l", col = col_wv, pch = 16)
  
  if (is.null(point_pch)){
    point_pch = 16
  }
  
  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$scales, x$variance, type = "p", col = col_wv, pch = point_pch, cex = point_cex)
}


#' @title Plot Wavelet Variance based on IMU Data
#' @description Displays a plot of wavelet variance accounting for CI values and supplied efficiency.
#' @method plot imu_wvar
#' @keywords internal
#' @param x                A \code{wvar} object.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_wv           A \code{string} that specifies the color of the wavelet variance line.
#' @param col_ci           A \code{string} that specifies the color of the confidence interval polygon.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param ci_wv            A \code{boolean} that determines whether a confidence interval polygon will be drawn.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of wavelet variance and confidence interval for each scale.
#' @author Stephane Guerrier and Yuming Zhang
#' @export
#' @examples 
#' data("kvh1750_wv")
#' plot(kvh1750_wv)
plot.imu_wvar = function(x, xlab = NULL, ylab = NULL, main = NULL,
                         col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                         ci_wv = NULL, point_cex = NULL, point_pch = NULL, ...){
  type = unique(x$type)
  if ("Gyroscope" %in% type){
    gyro_index = which(x$type == "Gyroscope")
  }else{
    gyro_index = NULL
  }
  if ("Accelerometer" %in% type){
    accel_index = which(x$type == "Accelerometer")
  }else{
    accel_index = NULL
  }
  
  ncol = length(unique(x$axis))
  nrow = length(type)
  m = length(x$wvar)
  J = length(x$wvar[[1]]$variance)
  
  # remove negative CI values
  index_to_remove = c()
  for (i in 1:m) {
    if(length(which(x$wvar[[i]]$lci<0)) > 0){
      index_to_remove = c(index_to_remove, which(x$wvar[[i]]$lci<0))
    }
  }
  if (!is.null(index_to_remove)){
    index_to_remove = unique(index_to_remove)
    index_to_keep = which(seq(1:J) != index_to_remove)
  }else{
    index_to_keep = 1:J
  }
  J = length(index_to_keep)
  scales = x$wvar[[1]]$scales[index_to_keep]
  ci_up = ci_lw = av = matrix(NA, J, m)
  for (i in 1:m){
    ci_up[,i] = x$wvar[[i]]$ci_high[index_to_keep]
    ci_lw[,i] = x$wvar[[i]]$ci_low[index_to_keep]
    av[,i] = x$wvar[[i]]$variance[index_to_keep]
  }
  
  # Axes
  if (is.null(nb_ticks_x)){
    nb_ticks_x = 6
  }
  if (is.null(nb_ticks_y)){
    nb_ticks_y = 5
  }
  
  # Range
  x_range = range(scales)
  x_low = floor(log10(x_range[1]))
  x_high = ceiling(log10(x_range[2]))
  x_ticks = seq(x_low, x_high, by = 1)
  if (length(x_ticks) > nb_ticks_x){
    x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
  }
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))
  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = "darkblue"
  }
  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }
  if (is.null(point_pch)){
    point_pch = 16
  }
  if (is.null(point_cex)){
    point_cex = 1.25
  }
  # Main Title
  if (is.null(main)){
    main = paste("Wavelet Variance Representation - ", x$sensor, " @ ", x$freq, " Hz", sep="")
  }
  # Labels
  if (is.null(xlab)){
    xlab = bquote(paste("Averaging time ", tau, " [sec]", sep = " "))
  }
  if (is.null(ylab)){
    ylab = expression(paste("Wavelet Variance ", nu, sep = ""))
  }
  
  # Main plot
  par(omi=rep(1, 4), mar=c(0,0,0,0), mfrow=c(nrow,ncol))
  
  # Gyro
  if (!is.null(gyro_index)){
    y_range = c(min(ci_lw[,gyro_index]), max(ci_up[,gyro_index]))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))
    y_ticks <- seq(y_low, y_high, by = 1)
    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
    for (i in seq_along(gyro_index)){
      plot(NA, xlim = range(scales), ylim = y_range, xaxt="n", yaxt="n", log = "xy", bty = "n")
      box(col = "grey")
      mtext(paste("Axis - ", x$axis[gyro_index][i], sep = ""), 3, line = 0.5)
      if (i == 1){
        axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2, cex = 1.25)
      }
      if (i == 1){
        mtext("Gyroscope", 2, line = 4.5)
        mtext(ylab, 2, line = 2.5)
      }
      abline(h = 10^y_ticks, col = "grey85")
      abline(v = 10^x_ticks, col = "grey85")
      # CI for AD
      if(ci_wv == TRUE || is.null(ci_wv)){
        polygon(c(scales, rev(scales)), c(ci_lw[,gyro_index[i]], rev(ci_up[,gyro_index[i]])),
                border = NA, col = col_ci)
      }
      # Add AD
      lines(scales, (av[,gyro_index[i]]), type = "l", col = col_wv, pch = 16)
      lines(scales, (av[,gyro_index[i]]), type = "p", col = col_wv, pch = point_pch, cex = point_cex)
      if (is.null(accel_index)){
        axis(1, at = 10^x_ticks, labels = x_labels, padj = -0.2, cex = 1.25)
      }
    }
  }
  # Accel
  if (!is.null(accel_index)){
    y_range = c(min(ci_lw[,accel_index]), max(ci_up[,accel_index]))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))
    y_ticks <- seq(y_low, y_high, by = 1)
    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
    for (i in seq_along(accel_index)){
      plot(NA, xlim = range(scales), ylim = y_range, xaxt="n", yaxt="n", log = "xy", bty = "n")
      box(col = "grey")
      if (i == 1){
        axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2, cex = 1.25)
      }
      if (i == 1){
        mtext("Accelerometer", 2, line = 4.5)
        mtext(ylab, 2, line = 2.5)
      }
      if (length(accel_index) == 3 && i == 2){
        mtext(xlab, 1, line = 3.5)
      }
      if (is.null(gyro_index)){
        mtext(paste("Axis - ", x$axis[gyro_index][i], sep = ""), 3, line = 0.5)
      }
      abline(h = 10^y_ticks, col = "grey85")
      abline(v = 10^x_ticks, col = "grey85")
      # CI for AD
      if(ci_wv == TRUE || is.null(ci_wv)){
        polygon(c(scales, rev(scales)), c(ci_lw[,accel_index[i]], rev(ci_up[,accel_index[i]])),
                border = NA, col = col_ci)
      }
      # Add AD
      lines(scales, (av[,accel_index[i]]), type = "l", col = col_wv, pch = 16)
      lines(scales, (av[,accel_index[i]]), type = "p", col = col_wv, pch = point_pch, cex = point_cex)
      axis(1, at = 10^x_ticks, labels = x_labels, padj = -0.2, cex = 1.25)
    }
  }
  # Add main title
  mtext(main, side = 3, line = 3, outer = TRUE)
  par(mfrow = c(1,1))
}



#' @title Comparison between classical and robust Wavelet Variances
#' @description Displays a plot of the wavelet variances (classical and robust) for a given time series accounting for CI values.
#' @param x               A time series objects.
#' @param eff             An \code{integer} that specifies the efficiency of the robust estimator.
#' @param units           A \code{string} that specifies the units of time plotted on the x axis.
#' @param xlab            A \code{string} that gives a title for the x axis.
#' @param ylab            A \code{string} that gives a title for the y axis.
#' @param main            A \code{string} that gives an overall title for the plot.
#' @param col_wv          A \code{string} that specifies the color of the wavelet variance line.
#' @param col_ci          A \code{string} that specifies the color of the confidence interval shade.
#' @param nb_ticks_x      An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y      An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param ...             Additional arguments affecting the plot.
#' @return Plot of wavelet variance and confidence interval for each scale.
#' @author Stephane Guerrier, Nathanael Claussen, and Justin Lee
#' @examples 
#' set.seed(999)
#' n = 10^4
#' Xt = rnorm(n)
#' wv = wvar(Xt)
#' 
#' plot(wv)
#' plot(wv, main = "Simulated white noise", xlab = "Scales")
#' plot(wv, units = "sec", legend_position = "topright")
#' plot(wv, col_wv = "darkred", col_ci = "pink")
#' @export
robust_eda = function(x, eff = 0.6, units = NULL, xlab = NULL, ylab = NULL, main = NULL, 
                      col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                      legend_position = NULL, ...){
  wv_cl  = wvar(x)
  wv_rob = wvar(x, robust = TRUE, eff = eff)
  
  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", tau, sep =""))
    }else{
      xlab = bquote(paste("Scale ", "", tau , " [", .(units), "]", sep = ""))
    }
  }
  
  if (is.null(ylab)){
    ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
  }else{
    ylab = ylab
  }
  
  # Main Title
  if (is.null(main)){
    main = "Classical vs Robust WV"
  }
  
  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = c("darkblue", "darkorange2")
  }
  
  if (is.null(col_ci)){
    col_ci = c(hcl(h = 210, l = 65, c = 100, alpha = 0.2), hcl(h = 60, l = 65, c = 100, alpha = 0.2))
  }
  
  # Range
  x_range = range(wv_cl$scales)
  x_low = floor(log2(x_range[1]))
  x_high = ceiling(log2(x_range[2]))
  
  y_range = range(c(wv_cl$ci_low, wv_cl$ci_high, wv_rob$ci_low, wv_rob$ci_high))
  y_low = floor(log10(y_range[1]))
  y_high = ceiling(log10(y_range[2]))
  
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
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))
  
  y_ticks <- seq(y_low, y_high, by = 1)
  if (length(y_ticks) > nb_ticks_y){
    y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
  }
  y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
  
  # Legend position
  if (is.null(legend_position)){
    if (which.min(abs(c(y_low, y_high) - log2(wv_rob$variance[1]))) == 1){
      legend_position = "topleft"
    }else{
      legend_position = "bottomleft"
    }
  }   
  
  # Main plot                     
  plot(NA, xlim = x_range, ylim = y_range, xlab = xlab, ylab = ylab, 
       log = "xy", xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE)
  win_dim = par("usr")
  
  par(new = TRUE)
  plot(NA, xlim = x_range, ylim = 10^c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
       xlab = xlab, ylab = ylab, log = "xy", xaxt = 'n', yaxt = 'n', bty = "n")
  win_dim = par("usr")
  
  # Add grid
  abline(v = 2^x_ticks, lty = 1, col = "grey95")
  abline(h = 10^y_ticks, lty = 1, col = "grey95")
  
  # Add title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = 10^c(win_dim[4], win_dim[4],
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]), 
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = 10^(win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)
  
  # Add Axes and Box
  lines(x_vec[1:2], rep(10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)
  y_ticks = y_ticks[(2^y_ticks) < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
  y_labels = y_labels[1:length(y_ticks)]
  box()
  axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2)  
  
  # CI for WV
  polygon(c(wv_cl$scales, rev(wv_cl$scales)), c(wv_cl$ci_low, rev(wv_cl$ci_high)),
          border = NA, col = col_ci[1])
  
  polygon(c(wv_rob$scales, rev(wv_rob$scales)), c(wv_rob$ci_low, rev(wv_rob$ci_high)),
          border = NA, col = col_ci[2])
  
  # Legend Position
  if (!is.na(legend_position)){
    if (legend_position == "topleft"){
      legend_position = 10^c(1.1*win_dim[1], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
      legend(x = legend_position[1], y = legend_position[2],
             legend = c("Classical WV", "Classical CI", "Robust WV", "Robust CI"),
             pch = c(16, 15, 16, 15), lty = c(1, NA, 1, NA), col = c(col_wv[1], col_ci[1], col_wv[2], col_ci[2]), 
             cex = 1, pt.cex = c(1.25, 3, 1.25, 3), bty = "n")
    }else{
      if (legend_position == "topright"){
        legend_position = 10^c(0.7*win_dim[2], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
        legend(x = legend_position[1], y = legend_position[2],
               legend = c("Classical WV", "Classical CI", "Robust WV", "Robust CI"),
               pch = c(16, 15, 16, 15), lty = c(1, NA, 1, NA), col = c(col_wv[1], col_ci[1], col_wv[2], col_ci[2]), 
               cex = 1, pt.cex = c(1.25, 3, 1.25, 3), bty = "n")
      }else{
        legend(legend_position,
               legend = c("Classical WV", "Classical CI", "Robust WV", "Robust CI"),
               pch = c(16, 15, 16, 15), lty = c(1, NA, 1, NA), col = c(col_wv[1], col_ci[1], col_wv[2], col_ci[2]), 
               cex = 1, pt.cex = c(1.25, 3, 1.25, 3), bty = "n")
      }
    }
  }
  
  lines(wv_cl$scales, wv_cl$variance, type = "l", col = col_wv[1], pch = 16)
  lines(wv_cl$scales, wv_cl$variance, type = "p", col = col_wv[1], pch = 16, cex = 1.25)
  
  lines(wv_cl$scales, wv_rob$variance, type = "l", col = col_wv[2], pch = 16)
  lines(wv_cl$scales, wv_rob$variance, type = "p", col = col_wv[2], pch = 16, cex = 1.25)
  
}



#' @title Multi-Plot Comparison Between Multiple Wavelet Variances
#' @description 
#' This is a helper function for the \code{compare_var()} function. 
#' This method accepts the same set of arguments as \code{compare_wvar} and returns a comparision 
#' of multiple wavelet variances of different time series accounting for CI values as a set of different plots.
#' 
#' @param graph_details       List of inputs
#' 
#' @author Stephane Guerrier, Justin Lee, and Nathanael Claussen
#' @export
#'
compare_wvar_split = function(graph_details){
  
  old_pars = par(mfrow = c(graph_details$obj_len, graph_details$obj_len),
             mar = c(0.5,0.5,0.5,1.5), oma = c(4,4,4,4))
  on.exit(par(old_pars))
  
  for (i in 1:graph_details$obj_len){
    for (j in 1:graph_details$obj_len){
      
      # Main plot                     
      plot(NA, xlim = graph_details$x_range, ylim = graph_details$y_range, log = "xy", xaxt = 'n', 
           yaxt = 'n', bty = "n", ann = FALSE)
      win_dim = par("usr")
      kill_y_tick = graph_details$y_at < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
      
      # Add grid
      abline(v = graph_details$x_at, lty = 1, col = "grey95")
      abline(h = graph_details$y_at, lty = 1, col = "grey95")
      
      # Add axes and box
      box(col = "grey")
      
      # Corner left piece
      if (j == 1){
        axis(2, at = graph_details$y_at[kill_y_tick], 
             labels = graph_details$y_labels[kill_y_tick], padj = -0.2, cex.axis = 1/log(graph_details$obj_len))
      }
      
      # Corner bottom
      if (i == graph_details$obj_len){
        
        axis(1, at = graph_details$x_at, labels = graph_details$x_labels, 
             padj = 0.1, cex.axis = 1/log(graph_details$obj_len))           # figure out how to size these things for smaller plots
      }
      # Diag graph
      if (i == j){
        scales   = graph_details$obj_list[[i]]$scales
        ci_low   = graph_details$obj_list[[i]]$ci_low
        ci_high  = graph_details$obj_list[[i]]$ci_high
        variance = graph_details$obj_list[[i]]$variance
        if(graph_details$ci_wv[i] == TRUE){  
          polygon(c(scales, rev(scales)), c(ci_low, rev(ci_high)),
                  border = NA, col = graph_details$col_ci[i])
        }
        lines(scales, variance, type = "l", col = graph_details$col_wv[i], pch = 16)
        lines(scales, variance, type = "p", col = graph_details$col_wv[i], 
              pch = 17, cex = graph_details$point_cex[i]/1.25)
        
        win_dim = par("usr")
        x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
        y_vec = 10^c(win_dim[4], win_dim[4],
                     win_dim[4] - 0.09*(win_dim[4] - win_dim[3]), 
                     win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
        
        
        box()
        
#        if (graph_details$add_legend){
#          if (i == j){
#            legend(graph_details$legend_position, graph_details$names, bty = "n",
#                   lwd = 1, pt.cex = graph_details$point_cex, pch = graph_details$point_pch,
#                   col = graph_details$col_wv, cex=0.7)
#          }
#        }
        
      }
      
      if (i != j){
        scales   = graph_details$obj_list[[i]]$scales
        ci_low   = graph_details$obj_list[[i]]$ci_low
        ci_high  = graph_details$obj_list[[i]]$ci_high
        variance = graph_details$obj_list[[i]]$variance
        
        win_dim = par("usr")
        x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
        y_vec = 10^c(win_dim[4], win_dim[4],
                     win_dim[4] - 0.09*(win_dim[4] - win_dim[3]), 
                     win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
        
        if (is.null(graph_details$main[i,j])){
          main = paste("WV:", graph_details$names[i], 
                       "vs", graph_details$names[j])
        }
        box()
        
        if (i < j && graph_details$ci_wv[i] == TRUE){
          polygon(c(scales, rev(scales)), c(ci_low, rev(ci_high)),
                  border = NA, col = graph_details$col_ci[i])
        }
        
        lines(scales, variance, type = "l", col = graph_details$col_wv[i], pch = 16)
        lines(scales, variance, type = "p", col = graph_details$col_wv[i], 
              pch = 17, cex = graph_details$point_cex[i]/1.25)
        
        scales   = graph_details$obj_list[[j]]$scales
        ci_low   = graph_details$obj_list[[j]]$ci_low
        ci_high  = graph_details$obj_list[[j]]$ci_high
        variance = graph_details$obj_list[[j]]$variance
        
        if (i < j && graph_details$ci_wv[i] == TRUE){ # don't show confidence intervals 
          polygon(c(scales, rev(scales)), c(ci_low, rev(ci_high)),
                  border = NA, col = graph_details$col_ci[j])
        }
        lines(scales, variance, type = "l", col = graph_details$col_wv[j], pch = 16)
        lines(scales, variance, type = "p", col = graph_details$col_wv[j], 
              pch = 17, cex = graph_details$point_cex[j]/1.25)
      }
      
      # Add Details
      # @todo: expand win_dim and position $names
      if(j==4){
        
        x_vec = 13^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
        
        #mtext(graph_details$names[i], side = 4, line = 0.1, cex = 0.8)
        
        par(xpd = TRUE) #Draw outside plot area
        text(x = x_vec[2], y = 0.03, graph_details$names[i], srt = 270, cex = 1.3, col = graph_details$col_wv[i])
        par(xpd = FALSE)
      }
      
      if(i==1){
        mtext(graph_details$names[j], side = 3, line = 0.2, cex = 0.8, col = graph_details$col_wv[j])
      }
    }
  }
  
  mtext(graph_details$ylab, side = 2, line = 2.80, cex = graph_details$cex_labels, outer = T)
  mtext(graph_details$xlab, side = 1, line = 2.75, cex = graph_details$cex_labels, outer = T)
  
  if (is.null(graph_details$main)){
    main = "Haar Wavelet Variance Representation"
  }else{
    mtext(main, side = 3, line = 1)
  }
}



#' @title Combined Plot Comparison Between Multiple Wavelet Variances
#' @description 
#' This is a helper function for the \code{compare_var()} function. 
#' This method accepts the same set of arguments as \code{compare_wvar} and returns a single plot
#' that compares multiple wavelet variances of different time series accounting for CI values.
#' 
#' @param graph_details       List of inputs
#' 
#' @author Stephane Guerrier, Justin Lee, and Nathanael Claussen
#' @export
#'
compare_wvar_no_split = function(graph_details){
  # Main plot                     
  plot(NA, xlim = graph_details$x_range, ylim = graph_details$y_range, log = "xy", xaxt = 'n', 
       yaxt = 'n', bty = "n", ann = FALSE)
  win_dim = par("usr")
  
  # Main Plot
  
  par(new = TRUE)
  plot(NA, xlim = graph_details$x_range, ylim = 10^c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
       log = "xy", xaxt = 'n', yaxt = 'n', bty = "n", 
       xlab = graph_details$xlab, ylab = graph_details$ylab,
       cex.lab = graph_details$cex_labels)
  win_dim = par("usr")
  
  # Add Grid
  abline(v = graph_details$x_at, lty = 1, col = "grey95")
  abline(h = graph_details$y_at, lty = 1, col = "grey95")
  
  # Add Title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = 10^c(win_dim[4], win_dim[4],
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]), 
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = 10^(win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), 
       graph_details$main)
  
  # Add Axes and Box
  lines(x_vec[1:2], rep(10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)
  y_ticks = graph_details$y_ticks[(10^graph_details$y_ticks) < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
  kill_y_tick = graph_details$y_at < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  
  box()
  axis(1, at = graph_details$x_at, labels = graph_details$x_labels, padj = 0.3)
  axis(2, at = graph_details$y_at[kill_y_tick], 
       labels = graph_details$y_labels[kill_y_tick], padj = -0.2)
  
  for (i in 1:graph_details$obj_len){
    scales   = graph_details$obj_list[[i]]$scales
    ci_low   = graph_details$obj_list[[i]]$ci_low
    ci_high  = graph_details$obj_list[[i]]$ci_high
    variance = graph_details$obj_list[[i]]$variance
    
    if (graph_details$ci_wv[i]){
      polygon(c(scales, rev(scales)), c(ci_low, rev(ci_high)),
              border = NA, col = graph_details$col_ci[i])
    }
    
    lines(scales, variance, type = "l", col = graph_details$col_wv[i], pch = 16)
    lines(scales, variance, type = "p", col = graph_details$col_wv[i], 
          pch = graph_details$point_pch[i], cex = graph_details$point_cex[i])
  }
  
  if (graph_details$add_legend){
    legend(graph_details$legend_position, graph_details$names, bty = "n",
           lwd = 1, pt.cex = graph_details$point_cex, pch = graph_details$point_pch,
           col = graph_details$col_wv)
  }
  
}



#' @title Comparison Between Multiple Wavelet Variances
#' @description 
#' Displays plots of multiple wavelet variances of different time series accounting for CI values.
#' 
#' @param ...             One or more time series objects.
#' @param split           A \code{boolean} that, if TRUE, arranges the plots into a matrix-like format.
#' @param add_legend      A \code{boolean} that, if TRUE, adds a legend to the plot.
#' @param units           A \code{string} that specifies the units of time plotted on the x axes. Note: This argument will not be used if xlab is specified.
#' @param xlab            A \code{string} that gives a title for the x axes.
#' @param ylab            A \code{string} that gives a title for the y axes.
#' @param main            A \code{string} that gives an overall title for the plot.
#' @param col_wv          A \code{string} that specifies the color of the wavelet variance lines. 
#' @param col_ci          A \code{string} that specifies the color of the confidence interval shade.
#' @param nb_ticks_x      An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y      An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param ci_wv           A \code{boolean} that determines whether confidence interval polygons will be drawn.
#' @param point_cex       A \code{double} that specifies the size of each symbol to be plotted.
#' @param point_pch       A \code{double} that specifies the symbol type to be plotted.
#' @param names           A \code{string} that specifies the name of the WVAR objects. 
#' @param cex_labels      A \code{double} that specifies the magnification of the labels (x and y).
#' @author Stephane Guerrier and Justin Lee
#' @export
#' @examples
#' set.seed(999)
#' n = 10^4
#' Xt = arima.sim(n = n, list(ar = 0.10))
#' Yt = arima.sim(n = n, list(ar = 0.35))
#' Zt = arima.sim(n = n, list(ar = 0.70))
#' Wt = arima.sim(n = n, list(ar = 0.95))
#' 
#' wv_Xt = wvar(Xt)
#' wv_Yt = wvar(Yt)
#' wv_Zt = wvar(Zt)
#' wv_Wt = wvar(Wt)
#' 
#' compare_wvar(wv_Xt, wv_Yt, wv_Zt, wv_Wt)
compare_wvar = function(... , split = FALSE, add_legend = TRUE, units = NULL, xlab = NULL, 
                        ylab = NULL, main = NULL, col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, 
                        nb_ticks_y = NULL, legend_position = NULL, ci_wv = NULL, point_cex = NULL, 
                        point_pch = NULL, names = NULL, cex_labels = 0.8){
  
  obj_list = list(...)
  obj_name = as.character(substitute(...()))
  obj_len  = length(obj_list)
  
  # Check if passed objects are of class wvar
  is_wvar = sapply(obj_list, FUN = is, class2 = 'wvar')
  if(!all(is_wvar == T)){
    stop("Supplied objects must be 'wvar' objects.")
  }
  
  # Check length of time series argument
  if (obj_len == 0){
    stop('No object given!')
  }else if (obj_len == 1){
    # -> plot.wvar
    plot.wvar(..., nb_ticks_x = nb_ticks_x, nb_ticks_y = nb_ticks_y)
  }else{
    
    if (is.null(xlab)){
      if (is.null(units)){
        xlab = expression(paste("Scale ", tau, sep =""))
      }else{
        xlab = bquote(paste("Scale ", "(", .(units), ")", sep = " "))
      }
    }else{
      xlab = xlab
    }
    
    if (is.null(ylab)){
      ylab = bquote(paste("Wavelet Variance ", nu^2, sep = " "))
    }else{
      ylab = ylab
    }
    
    if (is.null(ci_wv)){
      ci_wv = rep(TRUE, obj_len)
    }else{
      ci_wv = rep(ci_wv, obj_len)
    }
    
    # Main Title
    if (split == FALSE){
      if (is.null(main)){
        main = "Haar Wavelet Variance Representation"
      }
    }else{
      if (!is.null(main) && (dim(main)[1] != obj_len || dim(main)[2] != obj_len)){
        main = NULL
      }
    }
    
    hues = seq(15, 375, length = obj_len + 1)
    # Line Colors
    if (is.null(col_wv)){
      col_wv = hcl(h = hues, l = 65, c = 200, alpha = 1)[seq_len(obj_len)]
    }else{
      if (length(col_wv) != obj_len){
        col_wv = hcl(h = hues, l = 65, c = 200, alpha = 1)[seq_len(obj_len)]
      }
    }
    # CI Colors
    if (is.null(col_ci)){
      col_ci = hcl(h = hues, l = 80, c = 100, alpha = 0.2)[seq_len(obj_len)]
    }else{
      if (length(col_ci) != obj_len){
        col_ci = hcl(h = hues, l = 80, c = 100, alpha = 0.2)[seq_len(obj_len)]
      }
    }
    
    # X and Y Limits
    x_range = y_range = rep(NULL, 2)
    for (i in 1:obj_len){
      x_range = range(c(x_range, obj_list[[i]]$scales))
      y_range = range(c(y_range, obj_list[[i]]$ci_low, obj_list[[i]]$ci_high))
    }
    
    x_low = floor(log10(x_range[1]))
    x_high = ceiling(log10(x_range[2]))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))
    
    # Axes Labels and Ticks
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
    x_at = 10^x_ticks
    x_actual_length = sum((x_at < x_range[2])*(x_at > x_range[1]))
    if (x_actual_length < (3 + as.numeric(split == FALSE))){
      x_low = floor(log2(x_range[1]))
      x_high = ceiling(log2(x_range[2]))
      x_ticks = seq(x_low, x_high, by = 1)
      if (length(x_ticks) > 8){
        x_ticks = seq(x_low, x_high, by = 2)
      }
      x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))
      x_at = 2^x_ticks
    }
    
    y_ticks <- seq(y_low, y_high, by = 1)
    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels = sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
    y_at = 10^y_ticks
    
    # Legend position
    if (is.null(legend_position)){
      inter = rep(NA, obj_len)
      for (i in 1:obj_len){
        inter[i] = obj_list[[i]]$variance[1]
      }
      mean_wv_1 = mean(inter)
      if (which.min(abs(c(y_low, y_high) - log2(mean_wv_1))) == 1){
        legend_position = "topleft"
      }else{
        legend_position = "bottomleft"
      }
    }
    
    # Type of Points
    if (is.null(point_pch)){
      inter = rep(15:18, obj_len)
      point_pch = inter[1:obj_len]
    }else{
      if (length(point_pch) != obj_len){
        inter = rep(15:18, obj_len)
        point_pch = inter[1:obj_len]
      }
    }
    
    #Size of Points
    if (is.null(point_cex)){
      inter = rep(c(1.25,1.25,1.25,1.25), obj_len)
      point_cex = inter[1:obj_len]
    }else{
      if (length(point_pch) != obj_len){
        inter = rep(c(1.25,1.25,1.25,1.25), obj_len)
        point_cex = inter[1:obj_len]
      }
    }
    
    # Names of WVAR Objects
    if (is.null(names)){
      names = obj_name
    }else{
      if (length(names) != obj_len){
        names = obj_name
      }
    }
    
    # Arguments passed into compare_wvar_split or compare_wvar_no_split
    graph_details = list(obj_list = obj_list, obj_len = obj_len, names = names, xlab = xlab, 
                         ylab = ylab, col_wv = col_wv, add_legend = add_legend,
                         col_ci = col_ci, main = main, legend_position = legend_position,
                         ci_wv = ci_wv, point_cex = point_cex, point_pch = point_pch,
                         x_range = x_range, y_range = y_range, x_ticks = x_ticks, 
                         x_labels = x_labels, y_labels = y_labels, x_at = x_at, y_at = y_at,
                         y_ticks = y_ticks, nb_ticks_x = nb_ticks_x, nb_ticks_y = nb_ticks_y,
                         cex_labels = cex_labels)
    
    if (split == FALSE){
      # -> compare_wvar_no_split
      compare_wvar_no_split(graph_details)
    }else{
      # -> compare_wvar_split
      compare_wvar_split(graph_details)
    }
  }
}