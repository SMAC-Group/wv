using namespace Rcpp;
using std::cout;

#include <iostream>
#include <Rcpp.h>


// This is a simple example of exporting a C++ function to R.
// Source this function into an R session with 
// Rcpp::sourceCpp [function]


// [[Rcpp::export]]
int main() {
  for (int hashNum = 1; hashNum <= 5; hashNum++) {
    cout << "#";
  }
  cout << "\n";
  return 0;
}

/*** R
main()
  */


NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// R code block in C++ ...processed with sourceCpp
//

/*** R
timesTwo(42)
*/
