## Local check on Ubuntu 22.04

==> Rcpp::compileAttributes()

* Updated R/RcppExports.R

==> devtools::check(document = FALSE)

══ Building ════════════════════════════════════════════════════════════════════════════
Setting env vars:
• CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
• CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
• CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX14FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX17FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX20FLAGS: -Wall -pedantic -fdiagnostics-color=always
── R CMD build ─────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/home/lionel/github_repo/wv/DESCRIPTION’ ...
─  preparing ‘wv’:
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘wv_0.1.3.tar.gz’
   
══ Checking ════════════════════════════════════════════════════════════════════════════
Setting env vars:
• _R_CHECK_CRAN_INCOMING_USE_ASPELL_           : TRUE
• _R_CHECK_CRAN_INCOMING_REMOTE_               : FALSE
• _R_CHECK_CRAN_INCOMING_                      : FALSE
• _R_CHECK_FORCE_SUGGESTS_                     : FALSE
• _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
• NOT_CRAN                                     : true
── R CMD check ─────────────────────────────────────────────────────────────────────────
─  using log directory ‘/home/lionel/github_repo/wv.Rcheck’
─  using R version 4.5.1 (2025-06-13)
─  using platform: x86_64-pc-linux-gnu
─  R was compiled by
       gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
       GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
─  running under: Ubuntu 22.04.5 LTS
─  using session charset: UTF-8
─  using options ‘--no-manual --as-cran’
✔  checking for file ‘wv/DESCRIPTION’
─  checking extension type ... Package
─  this is package ‘wv’ version ‘0.1.3’
─  package encoding: UTF-8
✔  checking package namespace information
✔  checking package dependencies (811ms)
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files ...
✔  checking for hidden files and directories ...
✔  checking for portable file names ...
✔  checking for sufficient/correct file permissions
─  checking whether package ‘wv’ can be installed ... [54s/54s] OK (53.6s)
─  used C++ compiler: ‘g++ (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0’
─  checking installed package size ... INFO
     installed size is 10.6Mb
     sub-directories of 1Mb or more:
       libs   9.7Mb
✔  checking package directory ...
✔  checking for future file timestamps (55.6s)
✔  checking DESCRIPTION meta-information ...
✔  checking top-level files
✔  checking for left-over files
✔  checking index information
✔  checking package subdirectories ...
✔  checking code files for non-ASCII characters ...
✔  checking R files for syntax errors ...
✔  checking whether the package can be loaded (702ms)
✔  checking whether the package can be loaded with stated dependencies (585ms)
✔  checking whether the package can be unloaded cleanly (568ms)
✔  checking whether the namespace can be loaded with stated dependencies (594ms)
✔  checking whether the namespace can be unloaded cleanly (716ms)
✔  checking loading without being on the library search path (700ms)
✔  checking dependencies in R code (1.3s)
✔  checking S3 generic/method consistency (733ms)
✔  checking replacement functions (579ms)
✔  checking foreign function calls (672ms)
✔  checking R code for possible problems (5.6s)
✔  checking Rd files ...
✔  checking Rd metadata ...
✔  checking Rd line widths ...
✔  checking Rd cross-references ...
✔  checking for missing documentation entries (734ms)
✔  checking for code/documentation mismatches (1.8s)
✔  checking Rd \usage sections (978ms)
✔  checking Rd contents ...
✔  checking for unstated dependencies in examples ...
✔  checking contents of ‘data’ directory ...
✔  checking data for non-ASCII characters ...
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves ...
✔  checking line endings in C/C++/Fortran sources/headers ...
✔  checking line endings in Makefiles
✔  checking compilation flags in Makevars ...
✔  checking for GNU extensions in Makefiles ...
✔  checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS)
✔  checking use of PKG_*FLAGS in Makefiles
✔  checking use of SHLIB_OPENMP_*FLAGS in Makefiles ...
✔  checking pragmas in C/C++ headers and code ...
✔  checking compilation flags used
✔  checking compiled code ...
✔  checking examples (2.2s)
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
   
── R CMD check results ─────────────────────────────────────────────────── wv 0.1.3 ────
Duration: 2m 11.4s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded



## R-CMD-check on GitHub actions 

All jobs pass on 

- macOS-latest (release)
- windows-latest (release)
- ubuntu-latest (devel)
- ubuntu-latest (release)
- ubuntu-latest (oldrel-1)

see https://github.com/SMAC-Group/wv/actions/workflows/R-CMD-check.yaml


## Downstream dependencies

There are currently no downstream dependencies for this package.