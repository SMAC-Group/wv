#' @title Maximum Overlap Discrete Wavelet Transform
#' @description 
#' Calculates the coefficients for the discrete wavelet transformation
#' @param x        A \code{vector} with dimensions N x 1. 
#' @param nlevels  A \code{integer} indicating the \eqn{J} levels of decomposition.
#' @param filter   A \code{string} indicating the filter name
#' @return A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' The default \eqn{J} is determined by \eqn{floor\left(log_2 \left(length\left(x\right)\right)\right)}{floor(log2(length(x)))}
#' @author James Balamuta, Justin Lee and Stephane Guerrier
#' @export
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' ret = modwt(x)
#' 
#' summary(ret)
#' 
#' plot(ret)
modwt = function(x, nlevels = floor(log2(length(x)-1)), filter = "haar") {
  if (sum(class(x) %in% "gts") == 1){
    x = as.numeric(x)
  }
  ret = modwt_cpp(x = x, filter_name = filter, nlevels)
  mostattributes(ret) = list(J=nlevels, filter = filter, class=c("modwt","list"))
  ret
}

#' @title Print Maximum Overlap Discrete Wavelet Transform
#' @description
#' Prints the results of the modwt list
#' @method print modwt
#' @export
#' @param x A \code{modwt} object
#' @param ... further arguments passed to or from other methods.
#' @return Prints the modwt decomposition
#' @author James Balamuta and Nathanael Claussen 
#' @keywords internal
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' print(modwt(x))
print.modwt = function(x, ...){
  NextMethod("print")
}

#' @title Summary Maximum Overlap Discrete Wavelet Transform
#' @description 
#' Prints MODWT object in a concise format
#' @method summary modwt
#' @importFrom utils head
#' @export
#' @keywords internal
#' @param object A \code{modwt} object
#' @param ... additional arguments affecting the summary produced.
#' @return Prints the modwt matrix decomposition
#' @author Nathanael Claussen and Justin Lee
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' summary(modwt(x))
summary.modwt=function(object, ...) {
  cat("\n")
  cat("Results of MODWT using",attr(object,"filter"),"filter with",attr(object, "J"),"levels:\n")
  cat("Displaying only the first 6 coefficients...\n")
  y = as.list(object)
  j = length(y)
  for( i in 1:j ) {
    cat("Level",i,"Wavelet Coefficients\n", c(head(y[[i]])), "...\n")
  }
}

#' @title Plot Maximum Overlap Discrete Wavelet Transform
#' @description
#' Plots results of the modwt list in which additional parameters can be specified
#' @method plot modwt
#' @export
#' @param x A \code{modwt} object.
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
#' # MODWT
#' Yt = modwt(Xt)
#' 
#' # Graph examples
#' plot(Yt)
#' plot(Yt, index = c(1,4,5,6,8,2))
#' plot(Yt, index = c(1,4,5,6), couleur = "blue")
#' plot(Yt, index = c(1,4,5,6), couleur = rep(c("blue","yellow"),2))
plot.modwt = function(x, index = NULL, couleur = NULL, ...){
  J = attr(x,"J")
  
  if (is.null(index)){
    index = 1:(min(c(4,J)))
  }else{
      if (sum(index == "all") == 1){     # i.e. index == "all"
        index = 1:J
      }
    
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

