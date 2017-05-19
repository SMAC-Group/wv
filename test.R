# Load Library
library('rbenchmark')

set.seed(1)
x = rnorm(2^8)

# Run benchmark
out = benchmark(dwt(x), dwt.BrickWall(x))