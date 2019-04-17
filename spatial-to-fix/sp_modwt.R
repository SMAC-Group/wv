#' @title Maximum Overlap Discrete Wavelet Transform
#' @description 
#' Calculates the coefficients for the discrete wavelet transformation
#' @param x        A \code{matrix} with dimensions N x M. 
#' @param J1       A \code{integer} indicating the \eqn{J1} row levels of decomposition.
#' @param J2       A \code{integer} indicating the \eqn{J2} column levels of decomposition
#' @return A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @author Justin Lee
#' @export
sp_modwt = function(x, J1 = NULL, J2 = NULL) {
  if(is.null(J1)){
    J1 = floor(log2(dim(x)[1]-1))
  }
  if(is.null(J2)){
    J2 = floor(log2(dim(x)[2]-1))
  }
  ret = sp_modwt_cpp(X = x, J1 = J1, J2 = J2)
  # mostattributes(ret) = list(J=nlevels, filter = filter, class=c("modwt","list"))
  ret
}