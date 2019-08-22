# Rough implementation only based on first level decompoostion (set "detect = TRUE" and specify value for "w" between 0 and 1)

wavar = function(signal, nb.level = NA, robust = FALSE, eff=0.6, detect = FALSE, w=0, compute.v = FALSE, type.cov="theoretical") {
  
  # Length of the time series
  n <<- length(signal)
  
  # Compute number of scales considered
  if(is.na(nb.level)){
    nb.level = floor(log(n,2))
  }
  
  
  # MODWT transform
  signal.modwt <- waveslim::modwt(signal, "haar", nb.level)
  signal.modwt.bw <- waveslim::brick.wall(signal.modwt, "haar")
  
  # Compute wavelet variance  
  vmod = na.omit(wave.variance(signal.modwt.bw, type = "eta3")[1:nb.level,1])
  
  if(robust==TRUE){
    
    # Construct matrix with wavelet coef
    signal.modwt.rob <- matrix(unlist(signal.modwt.bw),n,(nb.level + 1))
    
    # Construct robust wavelet variance
    wv.rob.h = NA
    wv.rob.bw = NA
    wv.rob.p = NA
    
    for (i in 1:nb.level){
      # Non NA coef
      y = na.omit(signal.modwt.rob[,i])
      
      if(length(y)>1){
        
        # Robust estimators
        wv.rob.p[i] = percival(y)
        wv.rob.h[i] = sig.rob(eff,y)
        wv.rob.bw[i] = sig.rob.bw(eff,y)
        
        if(i==1){
          index = (abs(y/sqrt(wv.rob.bw[i])) > crob.bw)
          weight = ((1 - index) * (1 - (y/sqrt(wv.rob.bw[i])/crob.bw)^2)^2)
        }
        
      }
    }
    
    outliers = NA
    
    if(detect==TRUE){
      
      if(any(weight<=w)){
        
        ind = which(weight<=w)+0.5
        
        for(i in 1:length(ind)){
          
          pos = which(c(abs(signal[floor(ind[i])]),abs(signal[ceiling(ind[i])]))==max(abs(signal[floor(ind[i])]),abs(signal[ceiling(ind[i])])))
          outliers[i] = c(floor(ind[i]),ceiling(ind[i]))[pos]
          
        }
        
        outliers = unique(outliers)
        
      }
      
    }
    
    V = NA
    V.rob = NA
    
    if(compute.v == TRUE){
      
      if(type.cov=="theoretical"){
        
        nstar = sum(!is.na(signal.modwt.bw[[nb.level]]))
        wv.coeffs = lapply(signal.modwt.bw, function(d) d[n-nstar+1:n])
        wv.coeffs = do.call("rbind",lapply(wv.coeffs, function(d) d[!is.na(d)]))
        
        # Standard
        #V = compute.V(wv.coeffs)
        
        # Robust
        S.coeffs = matrix(NA,dim(wv.coeffs)[1]-1,dim(wv.coeffs)[2])
        for(i in 1:nb.level){
          S.coeffs[i,] = psi.tuk(wv.coeffs[i,],wv.rob.bw[i])
        }
        S = tcrossprod(S.coeffs[,1])
        for(i in 2:nstar){
          S = S + tcrossprod(S.coeffs[,i])
        }
        S = S/nstar
        M = matrix(0,nb.level,nb.level)
        for(k in 1:nstar){
          Mtemp = rep(NA,nb.level)
          for(i in 1:nb.level){
            temp = der.psi.tuk(wv.coeffs[i,k],wv.rob.bw[i])
            Mtemp[i] =  ifelse(is.na(temp),0,temp)   
          }
          M = M + diag(Mtemp)
        }
        M.temp = try(solve(-M/nstar),silent=T)
        if(class(M.temp)=="try-error"){
          M = ginv(-M/nstar)
        } else {
          M = M.temp
        }
        V.rob = M%*%S%*%t(M)
        
      } else {
        
        wv.coeffs = do.call("rbind",signal.modwt.bw)
        nstar = apply(wv.coeffs,1,function(d) sum(!is.na(d)))[1:nb.level]
        wv.coeffs[is.na(wv.coeffs)] = 0 
        
        # Standard
        #V = compute.V(wv.coeffs)
        
        # Robust
        S.coeffs = matrix(NA,dim(wv.coeffs)[1]-1,dim(wv.coeffs)[2])
        for(i in 1:nb.level){
          S.coeffs[i,] = psi.tuk(wv.coeffs[i,],wv.rob.bw[i])
        }
        S = tcrossprod(S.coeffs[,1])
        for(i in 2:n){
          S = S + tcrossprod(S.coeffs[,i])
        }
        S = S/n
        M = matrix(0,nb.level,nb.level)
        for(k in 1:n){
          Mtemp = rep(NA,nb.level)
          for(i in 1:nb.level){
            temp = der.psi.tuk(wv.coeffs[i,k],wv.rob.bw[i])
            Mtemp[i] =  ifelse(is.na(temp),0,temp)   
          }
          M = M + diag(Mtemp)
        }
        M.temp = try(solve(-M/n),silent=T)
        if(class(M.temp)=="try-error"){
          M = ginv(-M/n)
        } else {
          M = M.temp
        }
        V.rob = M%*%S%*%t(M)/n
        V.rob[V.rob==diag(V.rob)] = diag(V.rob)*n/nstar
        
      }
      
    }
    
  }
  
  if(robust==TRUE){
    wav.var = list(wv = vmod, wv.rob = wv.rob.bw, wv.mp = wv.rob.p, wv.hub = wv.rob.h, J = length(vmod), out = outliers, weight = weight, V = V, V.rob = V.rob)
  } else {
    wav.var = list(wv = vmod, wv.rob = NA, wv.mp = NA, wv.hub = NA, J = length(vmod), out = NA, V=V, V.rob = V.rob)
  }
  
}

signal <- rnorm(100)
nb.level <- floor(log(length(signal), 2))
signal.modwt <- waveslim::modwt(signal, "haar", nb.level)
signal.modwt.bw <- waveslim::brick.wall(signal.modwt, "haar")
signal.modwt.rob <- matrix(unlist(signal.modwt.bw), length(signal), (nb.level + 1))

detect_huber <- function(wv_coeffs, wv_vec, c, lim) {
  
  r = x/sqrt(sig2)
  w = pmax(0,pmin(1,crob/abs(r)))
  
}