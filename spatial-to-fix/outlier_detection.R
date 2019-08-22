# Rough implementation only based on first level decompoostion (set "detect = TRUE" and specify value for "w" between 0 and 1)

obj.fun.find.biwc <- function(crob, eff) {
  
  mu20c=1309458150*pnorm(crob)-2*crob*(654729075+218243025*crob^2+43648605*crob^4+6235515*crob^6+692835*crob^8+62985*crob^10+4845*crob^12+323*crob^14+19*crob^16+crob^18)*dnorm(crob)-654729075
  mu18c=68918850*pnorm(crob)-2*crob*(34459425+11486475*crob^2+2297295*crob^4+328185*crob^6+36465*crob^8+3315*crob^10+255*crob^12+17*crob^14+crob^16)*dnorm(crob)-34459425
  mu16c=4054050*pnorm(crob)-2*crob*(2027025+675675*crob^2+135135*crob^4+19305*crob^6+2145*crob^8+195*crob^10+15*crob^12+crob^14)*dnorm(crob)-2027025
  mu14c=270270*pnorm(crob)-2*crob*(135135+45045*crob^2+9009*crob^4+1287*crob^6+143*crob^8+13*crob^10+crob^12)*dnorm(crob)-135135
  mu12c=20790*pnorm(crob)-2*crob*(10395+3465*crob^2+693*crob^4+99*crob^6+11*crob^8+crob^10)*dnorm(crob)-10395
  mu10c=1890*pnorm(crob)-2*crob*(945+315*crob^2+63*crob^4+9*crob^6+crob^8)*dnorm(crob)-945
  mu8c=210*pnorm(crob)-2*crob*(105+35*crob^2+7*crob^4+crob^6)*dnorm(crob)-105
  mu6c=30*pnorm(crob)-2*crob*(15+5*crob^2+crob^4)*dnorm(crob)-15
  mu4c=6*pnorm(crob)-2*crob*(3+crob^2)*dnorm(crob)-3
  mu2c=2*pnorm(crob)-2*crob*dnorm(crob)-1
  ac=(1/crob^8)*mu10c-(4/crob^6)*mu8c+(6/crob^4)*mu6c-(4/crob^2)*mu4c+mu2c
  Q=(28/crob^4)*mu8c-(8/crob^2)*mu6c+(1/crob^16)*mu20c+(8/crob^14)*mu18c+(28/crob^12)*mu16c-(56/crob^10)*mu14c+(70/crob^8)*mu12c-(56/crob^6)*mu10c+mu4c-ac^2
  M=(1/crob^8)*mu12c-(4/crob^6)*mu10c+(6/crob^4)*mu8c-(4/crob^2)*mu6c+mu4c-ac
  out = ((0.5*M^2)/Q - eff)^2
  
  return(out)
  
}

find.biwc = function(eff){

  crob.t = optimize(obj.fun.find.biwc, c(0,1000), eff = eff)$minimum
  return(crob.t)
  
}

obj.fun.sig.rob.bw = function(sig2.bw, x, crob.bw, a.of.c) {
  
  r = x/sqrt(sig2.bw)
  ivec <- (abs(r) > crob.bw)
  w = ((1 - ivec) * (1 - (r/crob.bw)^2)^2)
  out = (mean(r^2*w^2) - a.of.c)^2
  
  return(out)
  
}

sig.rob.bw = function(eff, y) {
  
  crob.bw <- find.biwc(eff)
  x <- y/sd(y)
  a.of.c <- (1/crob.bw^8)*(1890*pnorm(crob.bw)-2*crob.bw*(945+315*crob.bw^2+63*crob.bw^4+9*crob.bw^6+crob.bw^8)*dnorm(crob.bw)-945)-(4/crob.bw^6)*(210*pnorm(crob.bw)-2*crob.bw*(105+35*crob.bw^2+7*crob.bw^4+crob.bw^6)*dnorm(crob.bw)-105)+(6/crob.bw^4)*(30*pnorm(crob.bw)-2*crob.bw*(15+5*crob.bw^2+crob.bw^4)*dnorm(crob.bw)-15)-(4/crob.bw^2)*(6*pnorm(crob.bw)-2*crob.bw*(3+crob.bw^2)*dnorm(crob.bw)-3)+2*pnorm(crob.bw)-2*crob.bw*dnorm(crob.bw)-1
  sig2.hat.rob.bw <- optimize(obj.fun.sig.rob.bw, c(0,2), x = x, crob.bw = crob.bw, a.of.c = a.of.c)$minimum*var(y)
  
  return(sig2.hat.rob.bw)
  
}

detect <- function(signal, eff = 0.6, w = 0, identify = FALSE) {
  
  # Length of the time series
  n <- length(signal)
  
  # Compute number of scales considered
  nb.level <- floor(log(n, 2))
  
  # MODWT transform
  signal.modwt <- waveslim::modwt(signal, "haar", nb.level)
  signal.modwt.bw <- waveslim::brick.wall(signal.modwt, "haar")
    
  # Construct matrix with wavelet coef
  signal.modwt.rob <- matrix(unlist(signal.modwt.bw), n, (nb.level + 1))
    
  # Construct robust wavelet variance
  wv.rob.bw <- NA
  weight <- vector("list", nb.level)
  crob.bw <- find.biwc(eff)
    
  for (i in 1:nb.level) {
      
    # Non NA coef
    y <- na.omit(signal.modwt.rob[, i])
      
    if(length(y) > 1) {
        
      # Robust estimators
      wv.rob.bw[i] <- sig.rob.bw(eff, y)
          
      index <- (abs(y/sqrt(wv.rob.bw[i])) > crob.bw)
      weight[[i]] <- ((1 - index) * (1 - (y/sqrt(wv.rob.bw[i])/crob.bw)^2)^2)
          
    }
    
  }
    
  outliers <- NA
    
  if(identify == TRUE) {
      
    if(any(weight[[1]] <= w)) {
        
      ind <- which(weight[[1]] <= w) + 0.5
        
        for(i in 1:length(ind)) {
          
          pos <- which(c(abs(signal[floor(ind[i])]), abs(signal[ceiling(ind[i])])) == max(abs(signal[floor(ind[i])]), abs(signal[ceiling(ind[i])])))
          outliers[i] <- c(floor(ind[i]), ceiling(ind[i]))[pos]
          
        }
        
        outliers <- unique(outliers)
        
      }
      
  }
  
  diagnostics <- list(weights = weight, outliers = outliers)
  class(diagnostics) = "detect"
  
  return(diagnostics)
    
}

plot.detect <- function(x) {
  
  quartz()
  par(mfrow = c(length(x$weights), 1), mar = c(3, 3, 1, 1))
  
  for(i in 1:length(x$weights)) {
    
    plot(x$weights[[i]], type = "l")
    
  }
  
}