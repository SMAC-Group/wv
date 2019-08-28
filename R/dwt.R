#' @title Discrete Wavelet Transform
#' @name dwt
#' @description 
#' Calculation of the coefficients for the discrete wavelet transformation
#' @export
#' @param x        A \code{vector} with dimensions N x 1. 
#' @param nlevels  A \code{integer} indicating the \eqn{J} levels of decomposition.
#' @param filter   A \code{string} indicating the filter name
#' @return A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' The default \eqn{J} is determined by \eqn{floor\left(log_2 \left(length\left(x\right)\right)\right)}{floor(log2(length(x)))}
#' @author James Balamuta, Justin Lee and Stephane Guerrier
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' ret = dwt(x)
#' 
#' summary(ret)
#' 
#' plot(ret)
dwt = function(x, nlevels = floor(log2(length(x))), filter = "haar") {
  if (sum(class(x) %in% "gts") == 1){
    x = as.numeric(x)
  }
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
#' @name print.dwt
#' @description
#' Prints the results of the modwt list
#' @method print dwt
#' @export
#' @param x A \code{dwt} object
#' @param ... further arguments passed to or from other methods.
#' @return Prints the dwt decomposition
#' @author James Balamuta and Nathanael Claussen 
#' @keywords internal
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' print(dwt(x))
print.dwt=function(x, ...){
  NextMethod("print")
}

#' @title Summary Discrete Wavelet Transform
#' @name summary.dwt
#' @description 
#' Prints DWT object in a concise format
#' @method summary dwt
#' @importFrom utils head
#' @export
#' @keywords internal
#' @param object A \code{dwt} object
#' @param ... additional arguments affecting the summary produced.
#' @return Prints the dwt matrix decomposition
#' @author Nathanael Claussen and Justin Lee
#' @examples
#' set.seed(999)
#' x = rnorm(2^8)
#' summary(dwt(x))
summary.dwt=function(object, ...) {
  cat("\n")
  cat("Results of DWT using",attr(object,"filter"),"filter with",attr(object, "J"),"levels:\n")
  cat("Displaying only the first 6 coefficients...\n")
  y = as.list(object)
  j = length(y)
  for( i in 1:j ) {
    cat("Level",i,"Wavelet Coefficients\n", c(head(y[[i]])), "...\n")
  }
}

#' @title Plot Discrete Wavelet Transform
#' @name plot.dwt
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
  J = attr(x,"J")
  
  if (is.null(index)){
    index = 1:(min(c(4,J)))
  }else{
    if (max(index) > J || min(index) < 1){
      stop("Incorrect index specified")
    }
  }
  
  nb_plot = length(index)
  
  if (is.null(couleur)){
    hues = seq(15, 375, length = nb_plot + 1)
    couleur = hcl(h = hues, l = 65, c = 100, alpha = 1)[seq_len(nb_plot)]
  }else{
    if (length(couleur) == 1 || length(couleur) != nb_plot){
      couleur = rep(couleur[1],nb_plot)
    }
  }
  
  old_pars = par(mfrow = c(nb_plot,1), mar = c(0,3,0,0), oma = c(5,2,1,1))
  on.exit(par(old_pars))
  
  x_range = length(x[[1]])
  for (i in seq_len(nb_plot)){
    current_time_series = x[[index[i]]]
    
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
