# check if packages are installed, and install them when needed 
inst_pkg = load_pkg = c('rbenchmark', 'microbenchmark', 'ggplot2', 'Rcpp')
inst_pkg = inst_pkg[!(inst_pkg %in% installed.packages()[,'Package'])]
if (length(inst_pkg)>0) install.packages(inst_pkg)
# load all necessary packages 
pkgs_loaded = lapply(load_pkg, require, character.only = TRUE)

sourceCpp("src/dwt.cpp")
source("src/modwt.R")

set.seed(1)
x = rnorm(10000)


# Test for equality 
a = modwt(x, nlevels = 2)
b = modwt_bw(x, nlevels = 2)
d = modwt_test(x, nlevels = 2)

all.equal(b, d)
# Result: Equal except for class name 


# Run benchmark
out = benchmark(modwt_bw(x), modwt_test(x))
out
# Results: 
# Note that modwt_test(x) runs relatively faster than modwt_bw(x) 
# Some significant relative difference. In terms of bigger datasets this may be hugely beneficial


# Run microbenchmark
out = microbenchmark(modwt_bw(x), modwt_test(x))
# Violin Plot
autoplot(out)



# Tests for "reflection" 
# d = modwt(x, boundary = "reflection", nlevels = 5)
# e = modwt_bw(x, boundary = "reflection", nlevels = 5)
