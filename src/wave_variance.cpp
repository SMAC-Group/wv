/* Copyright (C) 2014 - 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
 *
 * This file is part of GMWM R Methods Package
 *
 * The `gmwm` R package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The `gmwm` R package is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 */

#include <RcppArmadillo.h>
#include "wave_variance.h"
// Uses filters
#include "wv_filters.h"

// Uses brick_wall
#include "dwt.h"

// Uses robust components...
#include "robust_components.h"


//' @title Generate a Wave Variance for a Univariate Time Series
//' @description Computes an estimate of the wave variance
//' @param signal_modwt_bw A \code{field<vec>} that contains the brick walled modwt or dwt decomposition
//' @param robust          A \code{boolean} to determine the type of wave estimation.
//' @param eff             A \code{double} that indicates the efficiency.
//' @return A \code{vec} that contains the wave variance.
//' @keywords internal
//' @examples
//' set.seed(1337)
//' x = rnorm(100)
//' decomp = modwt_cpp(x, filter_name = "haar", nlevels = 4, boundary = "periodic", brickwall = TRUE)
//' wave_variance(decomp)
//' 
//' wave_variance(decomp, robust = TRUE, eff = 0.6)
// [[Rcpp::export]]
arma::vec wave_variance(const arma::field<arma::vec>& signal_modwt_bw, bool robust = false, double eff = 0.6){
  
  unsigned int nb_level = signal_modwt_bw.n_elem;
  arma::vec y(nb_level);
  
  if(robust){
    // Robust wavelet variance estimation
    for(unsigned int i=0; i < nb_level; i++){
      arma::vec wav_coef = sort(signal_modwt_bw(i));
      y(i) = sig_rob_bw(wav_coef, eff);
    }
  }else{
    // Classical wavelet variance estimation
    for(unsigned int i=0; i < nb_level;i++){
      arma::vec temp = signal_modwt_bw(i);
      y(i) = dot(temp,temp)/temp.n_elem;
    }
  }
  
  return y;
}


//' @title Computes the (MODWT) wavelet variance
//' @description Calculates the (MODWT) wavelet variance
//' @param signal_modwt_bw  A \code{field<vec>} that contains the modwt decomposition after it has been brick walled.
//' @param robust           A \code{boolean} that triggers the use of the robust estimate.
//' @param eff              A \code{double} that indicates the efficiency as it relates to an MLE.
//' @param alpha            A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @param ci_type          A \code{String} indicating the confidence interval being calculated. Valid value: "eta3"
//' @return A \code{mat} with the structure:
//' \itemize{
//'   \item{"variance"}{Wavelet Variance}
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' }
//' @keywords internal
//' @details 
//' This function does the heavy lifting with the signal_modwt_bw
//' @examples
//' x = rnorm(100)
//' decomp = modwt_cpp(x, filter_name = "haar", nlevels = 4, boundary = "periodic", brickwall = TRUE)
//' wvar_cpp(decomp, robust=FALSE, eff=0.6, alpha = 0.05, ci_type="eta3")
// [[Rcpp::export]]
arma::mat wvar_cpp(const arma::field<arma::vec>& signal_modwt_bw, 
                   bool robust, double eff, double alpha, 
                   std::string ci_type) {
  double alpha_ov_2 = alpha/2.0;

  // Wavelet Variance
  arma::vec y = wave_variance(signal_modwt_bw, robust, eff);
  
  // Confidence Interval
  return ci_wave_variance(signal_modwt_bw, y, ci_type, alpha_ov_2, robust, eff);
}



//' @title Computes the (MODWT) wavelet variance
//' @description Calculates the (MODWT) wavelet variance
//' @param signal     A \code{vec} that contains the data.
//' @param robust     A \code{boolean} that triggers the use of the robust estimate.
//' @param eff        A \code{double} that indicates the efficiency as it relates to an MLE.
//' @param alpha      A \code{double} that indicates the \eqn{\left(1-p\right)\times \alpha}{(1-p)*alpha} confidence level 
//' @param ci_type    A \code{string} indicating the confidence interval being calculated. Valid value: "eta3"
//' @param strWavelet A \code{string} indicating the type of wave filter to be applied. Must be "haar"
//' @param decomp     A \code{string} indicating whether to use "modwt" or "dwt" decomp
//' @return A \code{mat} with the structure:
//' \itemize{
//'   \item{"variance"}{Wavelet Variance}
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' }
//' @keywords internal
//' @details 
//' This function powers the wvar object. It is also extendable...
//' @examples
//' x=rnorm(100)
//' modwt_wvar_cpp(x, nlevels=4, robust=FALSE, eff=0.6, alpha = 0.05,
//'                ci_type="eta3", strWavelet="haar", decomp="modwt")
// [[Rcpp::export]]
arma::mat modwt_wvar_cpp(const arma::vec& signal, unsigned int nlevels, bool robust, double eff, double alpha, 
                         std::string ci_type, std::string strWavelet, std::string decomp) {
  
  // signal_modwt
  arma::field<arma::vec> signal_modwt_bw;
  if(decomp == "modwt"){
    signal_modwt_bw = modwt_cpp(signal, strWavelet, nlevels, "periodic", true);
  }else{
    signal_modwt_bw = dwt_cpp(signal, strWavelet, nlevels, "periodic", true);
  }

  arma::mat o = wvar_cpp(signal_modwt_bw,robust, eff, alpha, ci_type);
  
  return o;
}


//' @title Computes the MO/DWT wavelet variance for multiple processes
//' @description Calculates the MO/DWT wavelet variance
//' @param signal     A \code{matrix} that contains the same number of observations per dataset
//' @param robust     A \code{boolean} that triggers the use of the robust estimate.
//' @param eff        A \code{double} that indicates the efficiency as it relates to an MLE.
//' @param alpha      A \code{double} that indicates the \eqn{\left(1-p\right)\times \alpha}{(1-p)*alpha} confidence level 
//' @param ci_type    A \code{string} indicating the confidence interval being calculated. Valid value: "eta3"
//' @param strWavelet A \code{string} indicating the type of wave filter to be applied. Must be "haar"
//' @param decomp     A \code{string} indicating whether to use "modwt" or "dwt" decomp
//' @return A \code{field<mat>} with the structure:
//' \itemize{
//'   \item{"variance"}{Wavelet Variance}
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' }
//' @keywords internal
//' @details 
//' This function processes the decomposition of multiple signals quickly
//' @examples
//' x = cbind(rnorm(100),rnorm(100))
//' batch_modwt_wvar_cpp(x, nlevels=4, robust=FALSE, eff=0.6, 
//'                      alpha = 0.05, ci_type="eta3", strWavelet="haar", 
//'                      decomp="modwt")
// [[Rcpp::export]]
arma::field<arma::mat> batch_modwt_wvar_cpp(const arma::mat& signal, unsigned int nlevels, bool robust, double eff, double alpha, 
                                          std::string ci_type, std::string strWavelet, std::string decomp) {
  
  unsigned int num = signal.n_cols;
  // signal_modwt
  arma::field<arma::mat> wvars(num);
  for(unsigned int i = 0; i < num; i++){
    wvars(i) = modwt_wvar_cpp(signal.col(i), nlevels, robust, eff, alpha, 
                              ci_type, strWavelet, decomp);
    
  }
  
  return wvars;
}

//' @title Computes the MODWT scales
//' @description Calculates the MODWT scales
//' @param nb_level  A \code{integer} that contains the level of decomposition J.
//' @return A \code{vec} that contains 2^1, ... , 2^J
//' @keywords internal
//' @details 
//' Used in wvar object.
//' @examples
//' scales_cpp(5)
// [[Rcpp::export]]
arma::vec scales_cpp(unsigned int nb_level){
  // Define scales
  arma::vec scales(nb_level);
  for(unsigned int i=0; i< nb_level;i++){
    scales(i) = pow(2,i+1);
  }
  return scales;
}