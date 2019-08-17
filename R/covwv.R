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
  if (sum(class(x) %in% "gts") == 1){
    x = as.numeric(x)
  }
  
  if (sum(class(y) %in% "gts") == 1){
    y = as.numeric(y)
  }
  
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

  obj =  compute_cov_cpp(coef1 = coef1, coef2 = coef2, variance = variance, lower = lower, upper = upper)

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
    if (sum(class(x) %in% "gts") == 1){
      x = as.numeric(x)
    }
  
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
#' @method plot wccv_pair
#' @keywords internal
#' @export
#' @examples
#' n = 10^5
#' Xt = cumsum(rnorm(n, 0, 0.01))
#' Wt = Xt + rnorm(n)
#' Yt = Xt + rnorm(n)
#' wcov = wccv_pair(Wt, Yt)
#' plot(wcov)
plot.wccv_pair = function(x, theo.wccv = NULL, main = NULL, xlab = NULL, ylab = NULL,
                          units = NULL, col_wccv = NULL, col_ci = NULL,
                          nb_ticks_x = NULL, nb_ticks_y = NULL, ...){
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
  if(is.null(col_wccv)){
    col_wccv = "darkblue"
  }

  if(is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }

  if(is.null(main)){
    main = "Sample Wavelet Cross-Covariance"
  }


  # Axes
  if (is.null(nb_ticks_x)){
    nb_ticks_x = 4
  }

  if (is.null(nb_ticks_y)){
    nb_ticks_y = 7
  }

  # set ticks and labels
  tick_y_max  = ceiling(max(log10(abscombCI)))
  tick_y_min  = floor(min(log10(abscombCI)))
  tick_y_step =  2*(tick_y_max - tick_y_min)/(nb_ticks_y - 1)

  if (tick_y_step < 0.75){
    tick_y_step = 0.5
  }else{
    tick_y_step = round(tick_y_step)
  }

  y_at_lower = y_at_upper = y_at = seq(tick_y_min, tick_y_max, by = tick_y_step)

  upper_labels = sapply(y_at_upper, function(i) as.expression(bquote(10^ .(i))))
  lower_labels = sapply(y_at_lower, function(i) as.expression(bquote(-10^ .(i))))

  m = length(y_at_lower)
  ticks_y = c(-(m:1), 0, 1:m)
  labels = c(rev(lower_labels), 0, upper_labels)


  x_high = ceiling(log10(scales[J]))
  x_low = floor(log10(scales[1]))
  x_ticks = seq(x_low, x_high, by = 1)
  if (length(x_ticks) > nb_ticks_x){
    x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
  }
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))
  x_at = 10^x_ticks
  x_actual_length = sum((x_at < x_high)*(x_at > x_low))

  if (x_actual_length < 4){
    x_low = floor(log2(scales[1]))
    x_high = ceiling(log2(scales[J]))
    x_ticks = seq(x_low, x_high, by = 1)
    if (length(x_ticks) > 8){
      x_ticks = seq(x_low, x_high, by = 2)
    }
    x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))
    x_at = 2^x_ticks
  }

  plot(NA, log = "x", xlim = c(scales[1], scales[J]),
       ylim = c(min(ticks_y), max(1.09*ticks_y)), xaxt = "n", yaxt = "n",
       main = main, xlab = xlab, ylab = ylab, ann = FALSE, bty = "n")

  # Main plot
  win_dim = par("usr")

  par(new = TRUE)
  plot(NA, log = "x", xlim = c(scales[1], scales[J]),
       ylim = c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])), xaxt = "n", yaxt = "n",
       main = main, xlab = xlab, ylab = ylab, ann = FALSE, bty = "n")
  win_dim = par("usr")



  # Add grid
  abline(v = x_at, lty = 1, col = "grey95")
  abline(h = ticks_y, lty = 1, col = "grey95")
  abline(h = 0)

  # Add title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = c(win_dim[4], win_dim[4],
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = (win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)
  lines(x_vec[1:2], rep((win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)

  box()
  axis(2, at = ticks_y, labels = labels)
  axis(1, at = x_at, labels = x_labels)

  # Add CI
  y_low_ci  = wccv_get_y(x[,3], tick_y_min, tick_y_step)
  y_high_ci = wccv_get_y(x[,4], tick_y_min, tick_y_step)
  polygon(c(scales, rev(scales)), c(y_low_ci, rev(y_high_ci)),
          border = NA, col = col_ci)

  # Add wccv
  y_wccv = wccv_get_y(x[,1], tick_y_min, tick_y_step)
  lines(x = scales, y = y_wccv, type = "l", col = col_wccv, pch = 16, cex = 1.25)
  lines(x = scales, y = y_wccv, type = "p", col = col_wccv, pch = 16)

  # not sure what this is for
  if (is.null(theo.wccv) == F){
    # log.theo.positive = sapply(theo.wccv, function(x){ifelse(x < 0, NA, log(x))})
    # log.theo.negative = sapply(theo.wccv, function(x){ifelse(x > 0, NA, log(-x))})
    # lines(x = scales, y = log.theo.positive, lty = 3)
    lines(scales, wccv_get_y(theo.wccv, tick_y_min, tick_y_step), col="orange", lty = 3, lwd = 2)
  }

  # # not sure what this is for
  # if (is.null(theo.wccv) == F){
  #   lines(x = scales, y = -log.theo.negative, lty = 3)
  # }


}

#' @title Mapping to log10 scale
#' @description
#' Map x to the value in log10 scale
#' @export
#' @param x        A \code{vector} with dimensions J x 1.
#' @param tick_y_min A \code{negtive integer} the minimum power of 10, which corresponds to the smallest scale on y-axis.
#' @param tick_y_step An \code{integer} indicating the increment of the sequence.
#' @return A \code{field<vec>} that contains values in log10 scale.
#' @details
#' \code{tick_y_min} is usually chosen as \eqn{floor(min(log10(abs(x))))}
#' @author James Balamuta and Justin Lee
#' @examples
#' x = 2^(-1:-9)
#' y.min = floor(min(log10(abs(x))))
#' y.step = 2
#' wccv_get_y(x, y.min, y.step)
wccv_get_y = function(x, tick_y_min, tick_y_step){
  if (sum(class(x) %in% "gts") == 1){
    x = as.numeric(x)
  }
  
  n = length(x)
  for (i in 1:n){
    if (x[i] > 0){
      x[i] = (log10(x[i]) - tick_y_min)/tick_y_step + 1
    }else{
      x[i] = -(log10(abs(x[i])) - tick_y_min)/tick_y_step - 1
    }
  }
  x
}
