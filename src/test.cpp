#include <Rcpp.h>
using namespace Rcpp
using namespace std

// This is a simple example of exporting a C++ function to R.
// Source this function into an R session with 
// Rcpp::sourceCpp [function]


NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// R code block in C++ ...processed with sourceCpp
//

/*** R
timesTwo(42)
*/
