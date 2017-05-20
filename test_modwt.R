# check if packages are installed, and install them when needed 
inst_pkg = load_pkg = c('rbenchmark', 'microbenchmark', 'ggplot2', 'Rcpp')
inst_pkg = inst_pkg[!(inst_pkg %in% installed.packages()[,'Package'])]
if (length(inst_pkg)>0) install.packages(inst_pkg)
# load all necessary packages 
pkgs_loaded = lapply(load_pkg, require, character.only = TRUE)

set.seed(1)
x = rnorm(100)


# Test for equality 
a = modwt(x, nlevels = 5)
b = modwt_bw(x, nlevels = 5)

d = modwt(x, boundary = "reflection", nlevels = 5)
e = modwt_bw(x, boundary = "reflection", nlevels = 5)
