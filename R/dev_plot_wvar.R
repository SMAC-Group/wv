xlab = NULL; ylab = NULL; main = NULL;
col_ad = NULL; col_ci = NULL; nb_ticks_x = NULL; nb_ticks_y = NULL;
ci_ad = NULL; point_pch = NULL; point_cex = NULL

x = adis_wvar

#comparison with plot.avar
x2 = adis_av
library(avar)
plot(adis_av)
adis_av$avar[[1]]
adis_wvar$wvar[[1]]

plot.imu_wvar = function(x, xlab = NULL, ylab = NULL, main = NULL,
                         col_ad = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                         ci_ad = NULL, point_pch = NULL, point_cex = NULL, ...){
  type = unique(x$type)
  
  if ("Gyroscope" %in% type){
    gyro_index = which(x$type == "Gyroscope")
  }else{
    gyro_index = NULL
  }
  
  if ("Accelerometer" %in% type){
    accel_index = which(x$type == "Accelerometer")
  }else{
    accel_index = NULL
  }
  
  ncol = length(unique(x$axis))
  nrow = length(type)
  
  m = length(x$wvar)
  J = length(x$wvar[[1]]$variance)

  # remove negative CI values
  index_to_remove = c()
  for (i in 1:m) {
    x$wvar[[1]]$ci_low
    if(length(which(x$wvar[[i]]$ci_low<0)) > 0){
      index_to_remove = c(index_to_remove, which(x$wvar[[i]]$ci_low<0))
    }
  }
  if (!is.null(index_to_remove)){
    index_to_remove = unique(index_to_remove)
    index_to_keep = which(seq(1:J) != index_to_remove)
  }else{
    index_to_keep = 1:J
  }
  
  J = length(index_to_keep)
  scales = x$wvar[[1]]$scales[index_to_keep]

  ci_up = ci_lw = wvar_mat = matrix(NA, J, m)
  
  for (i in 1:m){
    ci_up[,i] = x$wvar[[i]]$ci_high[index_to_keep]
    ci_lw[,i] = x$wvar[[i]]$ci_low[index_to_keep]
    wvar_mat[,i] = x$wvar[[i]]$variance[index_to_keep]
  }
  
  # Axes
  if (is.null(nb_ticks_x)){
    nb_ticks_x = 6
  }
  
  if (is.null(nb_ticks_y)){
    nb_ticks_y = 5
  }
  
  # Range
  x_range = range(scales)
  x_low = floor(log10(x_range[1]))
  x_high = ceiling(log10(x_range[2]))
  
  x_ticks = seq(x_low, x_high, by = 1)
  if (length(x_ticks) > nb_ticks_x){
    x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
  }
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))
  
  
  # Line and CI colors
  if (is.null(col_ad)){
    col_ad = "darkblue"
  }
  
  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }
  
  if (is.null(point_pch)){
    point_pch = 16
  }
  
  if (is.null(point_cex)){
    point_cex = 1.25
  }
  
  # Main Title
  if (is.null(main)){
    main = paste("Allan Variance Representation - ", x$sensor, " @ ", x$freq, " Hz", sep="")
  }
  
  # Labels
  if (is.null(xlab)){
    xlab = bquote(paste("Averaging time ", tau, " [sec]", sep = " "))
  }
  
  if (is.null(ylab)){
    ylab = expression(paste("Allan Deviation ", phi, sep = ""))
  }
  
  
  # Main plot
  par(omi=rep(1.0, 4), mar=c(0,0,0,0), mfrow=c(nrow,ncol))
  
  # Gyro
  if (!is.null(gyro_index)){
    y_range = c(min(ci_lw[,gyro_index]), max(ci_up[,gyro_index]))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))
    
    y_ticks <- seq(y_low, y_high, by = 1)
    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
    
    
    for (i in seq_along(gyro_index)){
      plot(NA, xlim = range(scales), ylim = y_range, xaxt="n", yaxt="n", log = "xy", bty = "n")
      box(col = "grey")
      
      mtext(paste("Axis - ", x$axis[gyro_index][i], sep = ""), 3, line = 0.5)
      
      if (i == 1){
        axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2, cex = 1.25)
      }
      
      if (i == 1){
        mtext("Gyroscope", 2, line = 4.5)
        mtext(ylab, 2, line = 2.5)
      }
      
      abline(h = 10^y_ticks, col = "grey85")
      abline(v = 10^x_ticks, col = "grey85")
      
      # CI for AD
      if(ci_ad == TRUE || is.null(ci_ad)){
        polygon(c(scales, rev(scales)), c(ci_lw[,gyro_index[i]], rev(ci_up[,gyro_index[i]])),
                border = NA, col = col_ci)
      }
      
      # Add AD
      lines(scales, sqrt(av[,gyro_index[i]]), type = "l", col = col_ad, pch = 16)
      lines(scales, sqrt(av[,gyro_index[i]]), type = "p", col = col_ad, pch = point_pch, cex = point_cex)
      
      if (is.null(accel_index)){
        axis(1, at = 10^x_ticks, labels = x_labels, padj = -0.2, cex = 1.25)
      }
    }
  }
  
  # Accel
  if (!is.null(accel_index)){
    y_range = c(min(ci_lw[,accel_index]), max(ci_up[,accel_index]))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))
    
    y_ticks <- seq(y_low, y_high, by = 1)
    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
    
    
    for (i in seq_along(accel_index)){
      plot(NA, xlim = range(scales), ylim = y_range, xaxt="n", yaxt="n", log = "xy", bty = "n")
      box(col = "grey")
      if (i == 1){
        axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2, cex = 1.25)
      }
      
      
      if (i == 1){
        mtext("Accelerometer", 2, line = 4.5)
        mtext(ylab, 2, line = 2.5)
      }
      
      if (length(accel_index) == 3 && i == 2){
        mtext(xlab, 1, line = 3.5)
      }
      
      if (is.null(gyro_index)){
        mtext(paste("Axis - ", x$axis[gyro_index][i], sep = ""), 3, line = 0.5)
      }
      
      abline(h = 10^y_ticks, col = "grey85")
      abline(v = 10^x_ticks, col = "grey85")
      
      # CI for AD
      if(ci_ad == TRUE || is.null(ci_ad)){
        polygon(c(scales, rev(scales)), c(ci_lw[,accel_index[i]], rev(ci_up[,accel_index[i]])),
                border = NA, col = col_ci)
      }
      
      # Add AD
      lines(scales, sqrt(av[,accel_index[i]]), type = "l", col = col_ad, pch = 16)
      lines(scales, sqrt(av[,accel_index[i]]), type = "p", col = col_ad, pch = point_pch, cex = point_cex)
      axis(1, at = 10^x_ticks, labels = x_labels, padj = -0.2, cex = 1.25)
    }
  }
  
  # Add main title
  mtext(main, side = 3, line = 3, outer = TRUE)
  par(mfrow = c(1,1))
}


