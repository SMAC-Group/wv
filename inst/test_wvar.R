library(wv)
par(mfrow = c(1,2))

# Test white noise
Xt = rnorm(10000)
wv = wvar(Xt)
plot(wv, title = "White noise")
lines(wv$scales, 1/wv$scales, col = "darkorange")

# Test random walk
Yt = cumsum(Xt)
wv = wvar(Yt)
plot(wv, title = "Random walk")
lines(wv$scales, (wv$scales^2 + 2)/(12*wv$scales), col = "darkorange")
