install.packages("benchmark")
install.packages("microbenchmark")
install.packages("ggplot2")
# Load Libraries
library('rbenchmark')
library('microbenchmark')
library('ggplot2')

set.seed(1)
x = rnorm(2^8)

# Run benchmark
out = benchmark(dwt(x), dwt.BrickWall(x))
# Table Object
out

# Run microbenchmark
out = microbenchmark(dwt(x), dwt.BrickWall(x))
# Table Object
summary(out)
# Violin Plot
autoplot(out)

