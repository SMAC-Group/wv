# check if packages are installed, and install them when needed 
inst_pkg = load_pkg = c('rbenchmark', 'microbenchmark', 'ggplot2', 'Rcpp')
inst_pkg = inst_pkg[!(inst_pkg %in% installed.packages()[,'Package'])]
if (length(inst_pkg)>0) install.packages(inst_pkg)
# load all necessary packages 
pkgs_loaded = lapply(load_pkg, require, character.only = TRUE)

set.seed(1)
x = rnorm(2^16)

# Run benchmark
out = benchmark(dwt(x), dwt.BrickWall(x))
# Table Object
out

a = dwt(x)
b = dwt_bw(x)



# Run microbenchmark
out = microbenchmark(dwt(x), dwt.BrickWall(x))
# Table Object
summary(out)
# Violin Plot
autoplot(out)

