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
#' @export
#' @examples 
#' # Get Autocorrelation
#' m = ACF(datasets::AirPassengers)
#' 
#' # Get Autocovariance and do not remove trend from signal
#' m = ACF(datasets::AirPassengers, cor = FALSE, demean = FALSE)
#' @importFrom stats is.ts
#' @importFrom stats acf
ACF = function(x, lagmax = 0, cor = TRUE, demean = TRUE){
  
  if (sum(class(x) %in% "gts") == 1){
    x = as.numeric(x)
  }
  
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
  
  if (is.null(attr(x, "data_name"))){
    acfe = structure(acfe, n = nrow(x2), class = c("auto_corr", "array"))
  }else{
    acfe = structure(acfe, n = nrow(x2), main = attr(x, "data_name"), class = c("ACF", "array"))
  }
  
  acfe
  
}




#' @title Auto-Covariance and Correlation Functions
#' @description The acf function computes the estimated
#' autocovariance or autocorrelation for both univariate and multivariate cases.
#' @author Yunxiang Zhang
#' @param x         An \code{"ACF"} object from \code{\link{ACF}}.
#' @param show.ci   A \code{bool} indicating whether to show confidence region
#' @param ci        A \code{double} containing the 1-alpha level. Default is 0.95
#' @param ...       Additional parameters
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @rdname plot.auto_corr
#' @keywords internal
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
#' @importFrom grDevices rgb
plot.auto_corr = function(x, show.ci = TRUE, alpha = 0.05, main = NULL, ...){
  if (sum(class(x) %in% "gts") == 1){
    x = as.numeric(x)
  }
  
  # TO ADD AS INPUTS
  xlab = "Lags"
  ylab = "ACF"
  col_ci = rgb(0, 0.6, 1, 0.2)
  alpha = 0.05
  
  
  # Quiet the warnings...
  Lag = xmin = xmax = ymin = ymax = NULL 
  
  # Wide to long array transform
  x2 = as.data.frame.table(x, responseName = "ACF")
  
  colnames(x2) = c("Lag", "Signal X", "Signal Y", "ACF")
  
  # Remove character cast
  x2$Lag = as.numeric(x2$Lag)
  
  # Range
  x_range = range(x2$Lag)
  
  if (show.ci == TRUE){
    n = attr(x,"n")
    mult = qnorm(1-alpha/2)
    y_range = range(c(x2$ACF, 1/sqrt(n)*mult*c(-1,1)))
  }else{
    y_range = range(0:1)
  }
  
  
  x_ticks = seq(x_range[1], x_range[2], by = 1)
  y_ticks = seq(y_range[1], y_range[2], by = 0.05)
  
  old_pars = par(mar = c(5.1, 5.1, 1, 2.1))
  on.exit(par(old_pars))
  
  
  # Title 
  if (is.null(main)){
    if (is.null(attr(x,"data_name"))){
      main = paste0("ACF of ",as.character((x2$`Signal Y`)[1]))
    }else{
      main = paste0("ACF of ", attr(x,"data_name"))
    }
  }
  else {
    main = main
  }
  
  
  # Main plot
  plot(NA, xlim = c(1, max(x2$Lag)), ylim = y_range, 
       xlab = xlab, ylab = ylab, xaxt = 'n', 
       yaxt = 'n', bty = "n", ann = FALSE)
  win_dim = par("usr")
  
  par(new = TRUE)
  plot(NA, xlim = c(0, max(x2$Lag)), ylim = c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
       xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n', bty = "n")
  win_dim = par("usr")
  
  # Add grid
  grid(NULL, NULL, lty = 1, col = "grey95")
  
  
  # Add title
  
  x_vec = c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = c(win_dim[4], win_dim[4],
            win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
            win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = mean(c(win_dim[1], win_dim[2])),
       y = (win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), 
       main)
  
  # Add axes and box
  lines(x_vec[1:2], rep((win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)
  box()
  axis(1, padj = 0.3)
  y_axis = axis(2, labels = FALSE, tick = FALSE)
  y_axis = y_axis[y_axis < (win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
  axis(2, padj = -0.2, at = y_axis)
  
  abline(h = 0, lty = 1, lwd = 2)
  # Plot CI 
  if(show.ci){
    
    clim0 = 1/sqrt(n)*mult
    rect(xleft = -2, ybottom = -clim0, xright = 2*x_range[2], 
         ytop = clim0, col = col_ci, lwd = 0)
    
  }
  
  
  
  
  # Plot ACF
  segments(x0 = x_ticks, y0 = rep(0, x_range[2]), x1 = x_ticks, y1 = x2$ACF, lty = 1, lwd = 1)
  
  
}







