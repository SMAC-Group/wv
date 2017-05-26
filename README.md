
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Linux Build Status](https://travis-ci.org/SMAC-Group/wv.svg?branch=master)](https://travis-ci.org/SMAC-Group/wv)

`wv` R Package <img src="man/figures/logo.png" align="right" />
===============================================================

This repository holds the Wavelet Variance (wv) R package. This estimation technique computes the classical and robust wavelet variance for time series and regular lattices.

Below are examples of the capabilities of the `wv` package.

Install Instructions
====================

To install the `gmwm` package, there is currently one option: GitHub (Developmental).

Recommended R Interface
-----------------------

We firmly recommend that any users of this package use the [RStudio IDE](https://www.rstudio.com/products/rstudio/download/) over the default R GUI.

**All Systems**

With the system dependency taken care of, we continue on by installing the R specific package dependencies and finally the package itself by doing the following in an R session:

``` r
# Install dependencies
install.packages(c("RcppArmadillo","ggplot2","reshape2","devtools","knitr","rmarkdown"))

# Install the package from GitHub without Vignettes/User Guides
devtools::install_github("SMAC-Group/wv")
```

(Maximum Overlap) Discrete Wavelet Transform
============================================

For the moment this only works with Haar wavelets

``` r
# Simulate a Gaussian random walk
n = 10^3
Xt = cumsum(rnorm(n))
# DWT
# TODO:
# - summary
# - plot
Xt.dwt = dwt(Xt)
#> Warning in dwt(Xt): The data has been truncated so that it is divisible by
#> `nlevels` (e.g. 2^*)

# MODWT
# TODO:
# - summary
Xt.modwt = modwt(Xt)
#plot(Xt.modwt)
```

User Guides
===========

Various guides ship with package or are available on <http://smac-group.com/> to provide insight into how to use the different methods. At the present time, the following vignettes are available:

1.  Process to Haar Wavelet Variance [(Online)](https://smac-group.com/computing/2016/05/23/process-to-haar-wavelet-variance-formulae.html)
