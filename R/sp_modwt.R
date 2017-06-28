## Function to compute filter for Spatial decomposition
sp_hfilter = function(jscale){
  
  g = (1/sqrt(2))*c(1,1)/sqrt(2)
  h = c(1/sqrt(2),-1/sqrt(2))/sqrt(2)
  L = 2
  
  if(jscale==1) hup=h
  
  if(jscale >1){
    zero=c(rep(0,(2^(jscale-1)-1)))
    hup=h[1]
    for(i in 1:(L-1)) hup=c(hup,zero,h[(i+1)])
    
    for( j in 0:(jscale-2)){
      if(j==0) gup=g
      
      zero=c(rep(0,(2^j-1)))
      temp=g[1]
      for(i in 1:(L-1)) temp=c(temp,zero,g[(i+1)]) 
      
      if(j >0){
        sala=rep(0,(length(gup)+length(temp)-1))
        for(k in 1:length(sala)){
          dummy=0
          for( u in 1:length(gup)){
            if((k-u+1)>0 && (k-u+1)<=length(temp)) dummy=dummy + gup[u]*temp[(k-u+1)]
          }
          sala[k]=dummy
        }
        gup=sala
      } 
    }
    
    sala=rep(0,(length(hup)+length(gup)-1))
    for(k in 1:length(sala)){
      dummy=0
      for( u in 1:length(hup)){
        if((k-u+1)> 0 && (k-u+1)<=length(gup)) dummy=dummy+hup[u]*gup[(k-u+1)]
      }
      sala[k]=dummy
    }
    hup=sala      
  }
  
  hup
}

## General function to compute wavelet coefficients for spatial cases 
sp_modwt = function(X, J1, J2, compute.v = FALSE){
  
  n = dim(X)[1]
  m = dim(X)[2]
  wv = matrix(NA,J1,J2)
  wv.rob = matrix(NA,J1,J2)
  wv.mp = matrix(NA,J1,J2)
  sig.ci = matrix(NA,J1,J2)
  mn = matrix(NA,J1,J2)
  CI = matrix(NA,min(J1,J2),3)
  CI.rob = matrix(NA,min(J1,J2),3)
  nb.level = J1*J2
  
  if(compute.v==TRUE) wv.coeffs = vector("list",nb.level)
  
  i = 0
  k = 0
  
  for(j1 in 1:J1){
    for(j2 in 1:J2){
      
      hfil1 = sp_hfilter(j1)
      hfil2 = sp_hfilter(j2)
      mm1 = 2^j1
      mm2 = 2^j2
      pp1 = n - mm1 + 1
      pp2 = m - mm2 + 1
      xh = xhh = matrix(NA,n,n)
      
      # perhaps not n but m 
      for(spt in 1:n){
        for(tpt in 1:pp1) {
          
          xts=X[tpt:(tpt+mm1-1),spt] 
          xh[ (tpt+mm1%/%2), spt] = sum(xts*hfil1)
          
        }
      }
      
      for(tpt in 1:pp1){
        for(spt in 1:pp2) {
          
          xts=xh[(tpt+mm1%/%2),spt:(spt+mm2-1)] 
          xhh[(tpt+mm1%/%2), (spt+mm2%/%2)] = sum(xts*hfil2)
          
        }
      }
      
      i = i+1
      
      wv.coeff = c(xhh)
      wv.coeff = wv.coeff[!is.na(wv.coeff)]
      
      if(compute.v==TRUE) wv.coeffs[[i]] = wv.coeff
      
      # Standard
      wv[j1,j2] = sum(wv.coeff^2,na.rm=T)/(pp1*pp2)
      
      # Robust
      wv.mp[j1,j2] = percival(wv.coeff)
      wv.rob[j1,j2] = sig.rob.bw(eff,wv.coeff)
      
    }
    
  }
}