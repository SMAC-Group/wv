## Test environments

* local OS X install, R 3.6.1
* local Windows install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* rhub (linux, windows, solaris)

## R CMD check results

There was 1 NOTE:

* checking installed package size ... NOTE
  
It appears that within the Windows and Linux architectures, the CHECK procedure returns only one NOTE regarding the fact that the libs subdirectory is beyond the 1MB threshold. However, this NOTE doesn't occur to the OS X. Our understanding is that this size inflation of the libs subdirectory is due to the use of the Rcpp package. Indeed, some functions of the simts package have been written in C++ using Rcpp without which various functions would lose a considerable amount of computational efficiency leading to major parts of the package becoming impractical to use.

## Downstream dependencies

There are currently no downstream dependencies for this package.