## General function to compute wavelet coefficients for spatial cases 
spat.wavar = function(X, J1, J2, compute.v = FALSE){
  
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
      
      hfil1 = hfilter(j1)
      hfil2 = hfilter(j2)
      mm1 = 2^j1
      mm2 = 2^j2
      pp1 = n - mm1 + 1
      pp2 = m - mm2 + 1
      xh = xhh = matrix(NA,n,n)
      
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