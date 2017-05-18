set.seed(1)
x = rnorm(2^8)
me = dwt_cpp(x, filter_name = "haar", nlevels = 4, boundary = "periodic", brickwall = FALSE)
you = brick_wall(me, select_filter("haar"), "dwt")