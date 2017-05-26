#' Get ggplot2-like colors
#'
#' @param n number of colors.
#' @param alpha transparency.
#' @return list of colors.
ggplot_like_colors <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}


#' Plot function for dwt objects.
#' 
#' @param \code{y} A \code{dwt} object.
#' @param \code{index} A \code{vector} containing the indices to scales to be included in 
#' the graph. By default \code{index = 1:(min(c(J,4)))}, where \code{J} denotes the 
#' number of scales in \code{y}.
#' @param \code{couleur} A \code{vector} of colors of the same size as \code{index} used
#' for the different scales depicted in the graph. If \code{couleur} contains a single 
#' value thn the same color will be used for all scales.
#' @author Justin Lee and Stéphane Guerrier
#' @examples 
#' # Simulate a Gaussian white noise
#' n = 10^3
#' Xt = rnorm(n)
#' 
#' # dwt
#' Xt.dwt = dwt(Xt)
#' 
#' # Graph examples
#' plot(Xt.dwt)
#' plot(Xt.dwt, index = c(1,4,5,6,8,2))
#' plot(Xt.dwt, index = c(1,4,5,6), couleur = "blue")
#' plot(Xt.dwt, index = c(1,4,5,6), couleur = rep(c("blue","yellow"),2))
#' @export

plot.dwt = function(object, index = NULL, couleur = NULL){
  J <- attr(object,"J")
  
  if (is.null(index)){
    index <- 1:(min(c(4,J)))
  }else{
    if (max(index) > J || min(index) < 1){
      print("dsfa")
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
  x_range <- length(object[[1]])
  for (i in seq_len(nb_plot)){
    current_time_series <- object[[index[i]]]
    
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


#' Plot function for modwt objects.
#' 
#' @param \code{y} A \code{modwt} object.
#' @param \code{index} A \code{vector} containing the indices to scales to be included in 
#' the graph. By default \code{index = 1:(min(c(J,4)))}, where \code{J} denotes the 
#' number of scales in \code{y}.
#' @param \code{couleur} A \code{vector} of colors of the same size as \code{index} used
#' for the different scales depicted in the graph. If \code{couleur} contains a single 
#' value thn the same color will be used for all scales.
#' @author Justin Lee and Stéphane Guerrier
#' @examples 
#' # Simulate a Gaussian white noise
#' n = 10^3
#' Xt = rnorm(n)
#' 
#' # MODWT
#' Xt.modwt = modwt(Xt)
#' 
#' # Graph examples
#' plot(Xt.modwt)
#' plot(Xt.modwt, index = c(1,4,5,6,8,2))
#' plot(Xt.modwt, index = c(1,4,5,6), couleur = "blue")
#' plot(Xt.modwt, index = c(1,4,5,6), couleur = rep(c("blue","yellow"),2))
#' @export

plot.modwt = function(object, index = NULL, couleur = NULL){
  J <- attr(object,"J")

  if (is.null(index)){
    index <- 1:(min(c(4,J)))
  }else{
    if (max(index) > J || min(index) < 1){
      print("dsfa")
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
  x_range <- length(object[[1]])
  for (i in seq_len(nb_plot)){
    current_time_series <- object[[index[i]]]
    
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

