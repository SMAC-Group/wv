
#include <RcppArmadillo.h>
#include "wave_variance.h"
// Uses filters
#include "wv_filters.h"
#include "dwt.h"

// Uses robust components...
#include "robust_components.h"

//' @title Generate eta3 confidence interval
//' @description Computes the eta3 CI 
//' @param y          A \code{vec} that computes the modwt dot product of each wavelet coefficient divided by their length.
//' @param dims       A \code{String} indicating the confidence interval being calculated.
//' @param alpha_ov_2 A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
//' @return A \code{matrix} with the following columns:
//' \item{Column 1}{Wavelet Variance.}
//' \item{Column 2}{Chi-squared Lower Bounds.}
//' \item{Column 3}{Chi-squared Upper Bounds.}
//' @keywords internal
// [[Rcpp::export]]
arma::mat ci_eta3(const arma::vec& y, const arma::vec& dims, double alpha_ov_2) {
    
    unsigned int num_elem = dims.n_elem;

    arma::mat out(num_elem, 3);

    for(unsigned int i = 0; i<num_elem;i++){
      double eta3 = std::max(dims(i)/pow(2.0,double(i+1)),1.0);
      out(i,1) = eta3 * y(i)/R::qchisq(1-alpha_ov_2, eta3, 1, 0); // Lower CI
      out(i,2) = eta3 * y(i)/R::qchisq(alpha_ov_2, eta3, 1, 0); // Upper CI
    }

    out.col(0) = y;

    return out;
}

//' @title Generate eta3 robust confidence interval
//' @description Computes the eta3 robust CI
//' @param wv_robust   A \code{vec} that computes the modwt dot product of each wavelet coefficient divided by their length.
//' @param wv_ci_class A \code{mat} that contains the CI mean, CI Lower, and CI Upper
//' @param alpha_ov_2  A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level
//' @param eff         A \code{double} that indicates the efficiency.
//' @return A \code{matrix} with the following columns:
//'  \item{Column 1}{Robust Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' @details
//' Within this function we are scaling the classical 
//' @keywords internal
// [[Rcpp::export]]
arma::mat ci_eta3_robust(const arma::vec& wv_robust, const arma::mat& wv_ci_class, double alpha_ov_2, double eff) {
    unsigned int num_elem = wv_robust.n_elem;

    arma::mat out(num_elem, 3);
    
    //double q1 = R::qnorm(1-alpha_ov_2, 0.0, 1.0, true, false);
    
    double coef = sqrt(1.0/eff);

    for(unsigned int i = 0; i<num_elem;i++){
      
      double wv_ci = wv_ci_class(i,0); // store WV for i-th case
      
      
      double lci = (wv_ci - wv_ci_class(i,1))/wv_ci; // lower classical ci scale
      double uci = (wv_ci_class(i,2) - wv_ci)/wv_ci; // upper classical ci scale
      
      
      double wv_ri = wv_robust(i); // store WV for i-th case
      
      // Avoid negative result
      lci = wv_ri - lci*coef*wv_ri;
      if(lci > 0){
        out(i,1) = lci;        
      }else{
        // Take the minimum from either Robust WV / 2 or a previous minimum value.
        // This logic should no longer be needed once we remove the J-1th scale. 
        out(i,1) = std::min(wv_ri / 2.0, arma::as_scalar(arma::min(out.submat(0, 1, i-1, 1)))); 
        // Replaced DPL_EPSILON (Drop interval to 0) and made graph look odd
      }

      out(i,2) = wv_ri + uci*coef*wv_ri;
    }

    out.col(0) = wv_robust;

    return out;
}

//' @title Generate a Confidence interval for a Univariate Time Series
//' @description Computes an estimate of the multiscale variance and a chi-squared confidence interval
//' @param signal_modwt_bw A \code{field<vec>} that contains the modwt or dwt decomposition
//' @param wv              A \code{vec} that contains the wave variance.
//' @param type            A \code{String} indicating the confidence interval being calculated.
//' @param alpha_ov_2      A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level.
//' @param robust          A \code{boolean} to determine the type of wave estimation.
//' @param eff             A \code{double} that indicates the efficiency.
//' @return A \code{matrix} with the following structure:
//'  \item{Column 1}{Wavelet Variance}
//'  \item{Column 2}{Chi-squared Lower Bounds}
//'  \item{Column 3}{Chi-squared Upper Bounds}
//' @keywords internal
//' @details 
//' This function can be expanded to allow for other confidence interval calculations.
// [[Rcpp::export]]
arma::mat ci_wave_variance(const arma::field<arma::vec>& signal_modwt_bw, const arma::vec& wv,
                            std::string type = "eta3", double alpha_ov_2 = 0.025, bool robust = false, double eff = 0.6){
    
  unsigned int nb_level = wv.n_elem;
  arma::vec dims(nb_level);
  
  for(unsigned int i = 0; i < nb_level; i++){
    dims(i) = signal_modwt_bw(i).n_elem;
  }

  arma::mat out(nb_level , 3);

  
  if(type == "eta3"){
      
      if(!robust){
        out = ci_eta3(wv, dims, alpha_ov_2);  
      }else{
        // per the WV Robust change... 
        // We need to obtain the classical CI first, then modify it.
        arma::vec wv_class = wave_variance(signal_modwt_bw, false, eff); // Requires the next function....
        arma::mat wv_ci_class = ci_eta3(wv_class, dims, alpha_ov_2);  // calculate the CI
    
        // wv is the wave robust
        out = ci_eta3_robust(wv, wv_ci_class, alpha_ov_2, eff);
      }
  }
  else{
    Rcpp::stop("The wave variance type supplied is not supported. Please use: eta3");
  }

  return out;
}


//' @title Generate a Wave Variance for a Univariate Time Series
//' @description Computes an estimate of the wave variance
//' @param signal_modwt_bw A \code{field<vec>} that contains the modwt or dwt decomposition
//' @param robust          A \code{boolean} to determine the type of wave estimation.
//' @param eff             A \code{double} that indicates the efficiency.
//' @return A \code{vec} that contains the wave variance.
//' @keywords internal
// [[Rcpp::export]]
arma::vec wave_variance(const arma::field<arma::vec>& signal_modwt_bw, bool robust = false, double eff = 0.6){
  
  unsigned int nb_level = signal_modwt_bw.n_elem;
  arma::vec y(nb_level);
  
  if(robust){
    // Robust wavelet variance estimation
    for(unsigned int i=0; i < nb_level; i++){
      arma::vec wav_coef = sort(signal_modwt_bw(i));
      y(i) = sig_rob_bw(wav_coef, eff)(0);
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
//' @return A \code{mat} with the following structure:
//'   \item{"variance"}{Wavelet Variance}
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' @keywords internal
//' @details 
//' This function does the heavy lifting with the signal_modwt_bw
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
//' @return A \code{mat} with the following structure:
//'   \item{"variance"}{Wavelet Variance}
//'   \item{"low"}{Lower CI}
//'   \item{"high"}{Upper CI}
//' @keywords internal
//' @details 
//' This function powers the wvar object. It is also extendable...
// [[Rcpp::export]]
arma::mat modwt_wvar_cpp(const arma::vec& signal, unsigned int nlevels, bool robust, double eff, double alpha, 
                         std::string ci_type, std::string strWavelet, std::string decomp) {
  
  // signal_modwt
  arma::field<arma::vec> signal_modwt_bw;
  if(decomp == "modwt"){
    signal_modwt_bw = modwt_cpp(signal, strWavelet, nlevels);
  }else{
    signal_modwt_bw = dwt_cpp(signal, strWavelet, nlevels);
  }

  arma::mat o = wvar_cpp(signal_modwt_bw, robust, eff, alpha, ci_type);
  
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
//' \item{variance}{Wavelet Variance.}
//' \item{low}{Lower Confidence Interval.}
//' \item{high}{Upper Confidence Interval.}
//' @keywords internal
//' @details 
//' This function processes the decomposition of multiple signals quickly
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
// [[Rcpp::export]]
arma::vec scales_cpp(unsigned int nb_level){
  // Define scales
  arma::vec scales(nb_level);
  for(unsigned int i=0; i< nb_level;i++){
    scales(i) = pow(2.0,double(i+1));
  }
  return scales;
}