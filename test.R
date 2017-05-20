# check if packages are installed, and install them when needed 
inst_pkg = load_pkg = c('rbenchmark', 'microbenchmark', 'ggplot2', 'Rcpp')
inst_pkg = inst_pkg[!(inst_pkg %in% installed.packages()[,'Package'])]
if (length(inst_pkg)>0) install.packages(inst_pkg)
# load all necessary packages 
pkgs_loaded = lapply(load_pkg, require, character.only = TRUE)

set.seed(1)
x = rnorm(2^20)

# Run benchmark
out = benchmark(dwt(x), dwt_bw(x))
# Table Object
out
# Results: 
# Note that dwt(x) runs relatively faster than dwt_bw(x) 
# Some significant relative difference. In terms of bigger datasets this may be hugely beneficial


# Test for equality 
all.equal(dwt(x),dwt_bw(x)) 
# Result: Equal except for class name 


# Run microbenchmark
out = microbenchmark(dwt(x), dwt_bw(x))
# Table Object
summary(out)
# Violin Plot
autoplot(out)

# Results: 
# Again, looking at the microbenchmark results, dwt(x) runs signficantly faster than dwt_bw(x) 
# We can assume that dwt(x) will be more suitable for our purposes in the future 