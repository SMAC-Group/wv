# Libraries for wavelets and tapers
library(waveslim)
library(RSEIS)
library(MASS)

mat2vec = function(A){
  index = lower.tri(A, diag = TRUE)
  A[index]
}

find.first = function(v) {
  na.length = sum(is.na(v))
  return(v[na.length + 1])
}

# Functions for robust estimation
obj.fun.find.hubc = function(crob){
  mu4c=6*pnorm(crob)-2*crob*(3+crob^2)*dnorm(crob)-3
  mu2c=2*pnorm(crob)-2*crob*dnorm(crob)-1
  ac=mu2c+2*(crob^2)*(1-pnorm(crob))
  Q=2*(crob^4)*(1-pnorm(crob))+mu4c-ac^2
  M=2*(crob^2)*(1-pnorm(crob)+crob*dnorm(crob))+mu4c-ac
  out = ((0.5*M^2)/Q - eff)^2
  return(out)
}

find.hubc = function(eff){
  eff <<- eff
  crob.h = optimize(obj.fun.find.hubc, c(0,15))$minimum
  return(crob.h)
}

obj.fun.find.biwc = function(crob){
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
  eff <<- eff
  crob.t = optimize(obj.fun.find.biwc, c(0,1000))$minimum
  return(crob.t)
}

obj.fun.sig.rob = function(sig2){
  r = x/sqrt(sig2)
  w = pmax(0,pmin(1,crob/abs(r)))
  out = (mean(r^2*w^2) - a.of.c)^2
  return(out)
}

sig.rob = function(eff,y,gauss=TRUE){
  crob <<- find.hubc(eff)
  x <<- y/sd(y)
  a.of.c <<- 2*pnorm(crob) - 1 - 2*crob*dnorm(crob) + 2*crob^2*(1 - pnorm(crob))
  sig2.hat.rob = optimize(obj.fun.sig.rob, c(0,2))$minimum*var(y)
  return(sig2.hat.rob)
}

obj.fun.sig.rob.bw = function(sig2.bw){
  r = x/sqrt(sig2.bw)
  ivec <- (abs(r) > crob.bw)
  w=((1 - ivec) * (1 - (r/crob.bw)^2)^2)
  out = (mean(r^2*w^2) - a.of.c)^2
  return(out)
}

sig.rob.bw = function(eff,y){
  crob.bw <<- find.biwc(eff)
  x <<- y/sd(y)
  a.of.c <<- (1/crob.bw^8)*(1890*pnorm(crob.bw)-2*crob.bw*(945+315*crob.bw^2+63*crob.bw^4+9*crob.bw^6+crob.bw^8)*dnorm(crob.bw)-945)-(4/crob.bw^6)*(210*pnorm(crob.bw)-2*crob.bw*(105+35*crob.bw^2+7*crob.bw^4+crob.bw^6)*dnorm(crob.bw)-105)+(6/crob.bw^4)*(30*pnorm(crob.bw)-2*crob.bw*(15+5*crob.bw^2+crob.bw^4)*dnorm(crob.bw)-15)-(4/crob.bw^2)*(6*pnorm(crob.bw)-2*crob.bw*(3+crob.bw^2)*dnorm(crob.bw)-3)+2*pnorm(crob.bw)-2*crob.bw*dnorm(crob.bw)-1
  sig2.hat.rob.bw = optimize(obj.fun.sig.rob.bw, c(0,2))$minimum*var(y)
  return(sig2.hat.rob.bw)
}

percival = function(x){
  Tn=log(median(x^2))
  beta=get.slepians(npoints=length(x),nwin=5)/length(x)
  J=t(beta)%*%sign(log(x^2)-Tn)
  mu=(t(J)%*%colSums(beta))/(t(colSums(beta))%*%colSums(beta))
  Ahat=mean((J-mu*colSums(beta))^2)
  muhat=Tn-2*log(qnorm(3/4))-(Ahat/(-2*dnorm(qnorm(3/4))*qnorm(3/4))^2)/(2*length(x))
  sighat=exp(muhat)
  return(sighat)
}

psi.tuk = function(x,sig2.bw){
  r = x/sqrt(sig2.bw)
  ivec <- (abs(r) > crob.bw)
  w=((1 - ivec) * (1 - (r/crob.bw)^2)^2)
  out = r^2*w^2 - (1/crob.bw^8)*(1890*pnorm(crob.bw)-2*crob.bw*(945+315*crob.bw^2+63*crob.bw^4+9*crob.bw^6+crob.bw^8)*dnorm(crob.bw)-945)-(4/crob.bw^6)*(210*pnorm(crob.bw)-2*crob.bw*(105+35*crob.bw^2+7*crob.bw^4+crob.bw^6)*dnorm(crob.bw)-105)+(6/crob.bw^4)*(30*pnorm(crob.bw)-2*crob.bw*(15+5*crob.bw^2+crob.bw^4)*dnorm(crob.bw)-15)-(4/crob.bw^2)*(6*pnorm(crob.bw)-2*crob.bw*(3+crob.bw^2)*dnorm(crob.bw)-3)+2*pnorm(crob.bw)-2*crob.bw*dnorm(crob.bw)-1
  return(out)
}

der.psi.tuk = function(x,sig2.bw){
  r = x/sqrt(sig2.bw)
  r2 = ((r^2/crob.bw^2-1)^3)*(r^2/sig2.bw)*(1-5*(r^2/crob.bw^2))
  ivec <- (abs(r) > crob.bw)
  r2 =(1 - ivec) * r2
  return(r2)
}

My.acf = function(x){
  x=x[!is.na(x)]
  return(my.acf(x))
}

compute.cov1 = function(W1,W2){
  # Remove missing values
  W1 = na.omit(W1)
  W2 = na.omit(W2)
  
  # Compute length to check integrity
  w1.length = length(W1)
  w2.length = length(W2) 
  
  # Compute lagMax
  lagMax = min(c(w1.length,w2.length))
  lagMin = max(c(w1.length,w2.length))
  
  # Compute covariance between two wavelet levels
  cross.cov.12 = ccf(W1,W2,type = "covariance", lag.max = lagMax, plot = FALSE)$acf
  cross.cov.21 = ccf(W2,W1,type = "covariance", lag.max = lagMax, plot = FALSE)$acf
  cov.W1W2 = (sum(as.vector(cross.cov.12)^2) + sum(as.vector(cross.cov.21)^2))/lagMin/2
  return(cov.W1W2)  
}

compute.cov = function(W1,W2){
  # Remove missing values
  W1 = na.omit(W1)
  W2 = na.omit(W2)
  
  # Compute length to check integrity
  w1.length = length(W1)
  
  # Compute covariance between two wavelet levels
  cross.cov = ccf(W1,W2,type = "covariance", lag.max = w1.length, plot = FALSE)$acf
  cov.W1W2 = as.vector(cross.cov)[1]^2/2 + sum(as.vector(cross.cov)[-1]^2)
  return(cov.W1W2)  
}

compute.V = function(sig.modwt){
  # Compute number of levels
  nb.lev = length(sig.modwt) - 1
  
  # Initialisation of covariance matrix
  full.V = matrix(NA,nb.lev,nb.lev)
  
  # Fill matrix elements
  for (i in 1:nb.lev){
    for (j in 1:i){
      Wi = as.vector(unlist(sig.modwt[i]))
      Wj = as.vector(unlist(sig.modwt[j]))
      full.V[i,j] = compute.cov(Wi,Wj)
      full.V[j,i] = full.V[i,j]
    }
  }
  return(full.V) 
}

compute.V.spat = function(sig.modwt){
  # Compute number of levels
  nb.lev = length(sig.modwt)
  
  # Initialisation of covariance matrix
  full.V = matrix(NA,nb.lev,nb.lev)
  
  # Fill matrix elements
  for (i in 1:nb.lev){
    for (j in 1:i){
      Wi = as.vector(unlist(sig.modwt[i]))
      Wj = as.vector(unlist(sig.modwt[j]))
      full.V[i,j] = compute.cov(Wi,Wj)
      full.V[j,i] = full.V[i,j]
    }
  }
  return(full.V) 
}

hfilter=function(jscale){
  
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

spat.wavar.aniso = function(X,J1,J2,eff=0.6,compute.v = FALSE,compute.ci=FALSE){
  
  n = dim(X)[1]
  m = dim(X)[2]
  wv = matrix(NA,J1,J2)
  wv.rob = matrix(NA,J1,J2)
  wv.mp = matrix(NA,J1,J2)
  sig.ci = matrix(NA,J1,J2)
  mn = matrix(NA,J1,J2)
  nb.level = J1*J2
  CI = matrix(NA,min(J1,J2),3)
  CI.rob = matrix(NA,min(J1,J2),3)
  
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
      
      if(compute.ci==TRUE && j1==j2){
        
        k = k+1
        A = compute.cov(wv.coeff,wv.coeff)
        CI[k,] = c(wv[j1,j2],(wv[j1,j2]*length(wv.coeff)*wv[j1,j2]^2/A)*(1/c(qchisq(0.975,length(wv.coeff)*wv[j1,j2]^2/A),qchisq(0.025,length(wv.coeff)*wv[j1,j2]^2/A))))
        A = A/eff #compute.cov(psi.tuk(wv.coeff,wv.rob[j1,j2]),psi.tuk(wv.coeff,wv.rob[j1,j2]))/mean(der.psi.tuk(wv.coeff,wv.rob[j1,j2]))^2
        CI.rob[k,] = c(wv.rob[j1,j2],(wv.rob[j1,j2]*length(wv.coeff)*wv.rob[j1,j2]^2/A)*(1/c(qchisq(0.975,length(wv.coeff)*wv.rob[j1,j2]^2/A),qchisq(0.025,length(wv.coeff)*wv.rob[j1,j2]^2/A))))
        
      }
      
    }
    
  }
  
  V = NA
  V.rob = NA
  
  if(compute.v==TRUE){
    
    nstar = length(wv.coeffs[[nb.level]])
    wv.coeffs = lapply(wv.coeffs, function(d) d[1:nstar])
    wv.coeffs = do.call("rbind",wv.coeffs)
    
    #Standard
    #V = compute.V.spat(wv.coeffs)
    
    # Robust
    S.coeffs = matrix(NA,dim(wv.coeffs)[1],dim(wv.coeffs)[2])
    i = 0
    for(j1 in 1:J1){
      for(j2 in 1:J2){
        i = i+1
        S.coeffs[i,] = psi.tuk(wv.coeffs[i,],wv.rob[j1,j2])
      }
    }
    S = tcrossprod(S.coeffs[,1])
    for(i in 2:nstar){
      S = S + tcrossprod(S.coeffs[,i])
    }
    S = S/nstar
    M = matrix(0,nb.level,nb.level)
    for(k in 1:nstar){
      Mtemp = rep(NA,nb.level)
      i = 0
      for(j1 in 1:J1){
        for(j2 in 1:J2){
          i = i+1
          temp = der.psi.tuk(wv.coeffs[i,k],wv.rob[j1,j2])
          Mtemp[i] =  ifelse(is.na(temp),0,temp)
        }
      }
      M = M + diag(Mtemp)
    }
    M.temp = try(solve(-M/nstar),silent=T)
    if(class(M.temp)=="try-error"){
      M = ginv(-M/nstar)
    } else {
      M = M.temp
    }
    V.rob = M%*%S%*%t(M)/nstar
    
  }
  
  return(list(wv=wv,wv.rob=wv.rob,wv.mp=wv.mp,V=V,V.rob=V.rob,CI=CI,CI.rob=CI.rob))
  
}

spat.wavar.iso = function(X,J1,J2,eff=0.6){
  
  n = dim(X)[1]
  m = dim(X)[2]
  wv = matrix(NA,J1,J2)
  wv.rob = matrix(NA,J1,J2)
  wv.mp = matrix(NA,J1,J2)
  sig.ci = matrix(NA,J1,J2)
  mn = matrix(NA,J1,J2)
  
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
      
      wv.coeff = c(xhh)
      wv.coeff = wv.coeff[!is.na(wv.coeff)]
      sig.ci[j1,j2] = sum(acf(wv.coeff, type="covariance", plot=F)$acf^2)/2
      
      # Standard
      wv[j1,j2] = sum(wv.coeff^2,na.rm=T)/(pp1*pp2)
      
      # Robust
      wv.mp[j1,j2] = percival(wv.coeff)
      wv.rob[j1,j2] = sig.rob.bw(eff,wv.coeff)
      
      mn[j1,j2] = length(wv.coeff)
      
    }
    
  }
  
  if(dim(wv)[1]<dim(wv)[2]){
    wv = t(wv)
    wv.rob = t(wv.rob)
    wv.mp = t(wv.mp)
  }
  
  WV.temp = wv[1:min(J1,J2),1:min(J1,J2)]
  WV.temp1 = (WV.temp[lower.tri(WV.temp,diag=TRUE)] + t(WV.temp)[t(upper.tri(WV.temp,diag=TRUE))])/2
  WV.temp[lower.tri(WV.temp,diag=TRUE)] = WV.temp1
  wv[1:min(J1,J2),1:min(J1,J2)] = WV.temp
  wv = wv[lower.tri(wv, diag=TRUE)]
  
  WV.temp = wv.rob[1:min(J1,J2),1:min(J1,J2)]
  WV.temp1 = (WV.temp[lower.tri(WV.temp,diag=TRUE)] + t(WV.temp)[t(upper.tri(WV.temp,diag=TRUE))])/2
  WV.temp[lower.tri(WV.temp,diag=TRUE)] = WV.temp1
  wv.rob[1:min(J1,J2),1:min(J1,J2)] = WV.temp
  wv.rob = wv.rob[lower.tri(wv.rob, diag=TRUE)]
  
  WV.temp = wv.mp[1:min(J1,J2),1:min(J1,J2)]
  WV.temp1 = (WV.temp[lower.tri(WV.temp,diag=TRUE)] + t(WV.temp)[t(upper.tri(WV.temp,diag=TRUE))])/2
  WV.temp[lower.tri(WV.temp,diag=TRUE)] = WV.temp1
  wv.mp[1:min(J1,J2),1:min(J1,J2)] = WV.temp
  wv.mp = wv.mp[lower.tri(wv.mp, diag=TRUE)]
  
  sig.temp = sig.ci[1:min(J1,J2),1:min(J1,J2)]
  sig.temp1 = (sig.temp[lower.tri(sig.temp,diag=TRUE)] + t(sig.temp)[t(upper.tri(sig.temp,diag=TRUE))])/2
  sig.temp[lower.tri(sig.temp,diag=TRUE)] = sig.temp1
  sig.ci[1:min(J1,J2),1:min(J1,J2)] = sig.temp
  sig.ci = sig.ci[lower.tri(sig.ci, diag=TRUE)]
  
  mn.temp = mn[1:min(J1,J2),1:min(J1,J2)]
  mn.temp1 = (mn.temp[lower.tri(mn.temp,diag=TRUE)] + t(mn.temp)[t(upper.tri(mn.temp,diag=TRUE))])/2
  mn.temp[lower.tri(mn.temp,diag=TRUE)] = mn.temp1
  mn[1:min(J1,J2),1:min(J1,J2)] = mn.temp
  mn = mn[lower.tri(mn, diag=TRUE)]
  
  return(list(wv=wv,wv.rob=wv.rob,wv.mp=wv.mp,sig=sig.ci, mn=mn))
  
}

spat.wavar = function(X,J1,J2,eff=0.6,compute.v = FALSE,compute.ci=FALSE, iso = TRUE){
  
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
      sig.ci[j1,j2] = sum(acf(wv.coeff, type="covariance", plot=F)$acf^2)/2
      mn[j1,j2] = length(wv.coeff)
      
      if(compute.v==TRUE) wv.coeffs[[i]] = wv.coeff
      
      # Standard
      wv[j1,j2] = sum(wv.coeff^2,na.rm=T)/(pp1*pp2)
      
      # Robust
      wv.mp[j1,j2] = percival(wv.coeff)
      wv.rob[j1,j2] = sig.rob.bw(eff,wv.coeff)
      
      if(compute.ci==TRUE && j1==j2){
        
        k = k+1
        A = compute.cov(wv.coeff,wv.coeff)
        CI[k,] = c(wv[j1,j2],(wv[j1,j2]*length(wv.coeff)*wv[j1,j2]^2/A)*(1/c(qchisq(0.975,length(wv.coeff)*wv[j1,j2]^2/A),qchisq(0.025,length(wv.coeff)*wv[j1,j2]^2/A))))
        A = A/eff #compute.cov(psi.tuk(wv.coeff,wv.rob[j1,j2]),psi.tuk(wv.coeff,wv.rob[j1,j2]))/mean(der.psi.tuk(wv.coeff,wv.rob[j1,j2]))^2
        CI.rob[k,] = c(wv.rob[j1,j2],(wv.rob[j1,j2]*length(wv.coeff)*wv.rob[j1,j2]^2/A)*(1/c(qchisq(0.975,length(wv.coeff)*wv.rob[j1,j2]^2/A),qchisq(0.025,length(wv.coeff)*wv.rob[j1,j2]^2/A))))
        
      }
      
    }
    
  }
  
  V = NA
  V.rob = NA
  
  if(iso==TRUE){
    
    nb.level = J1*(J2+1)/2
    
    if(dim(wv)[1]<dim(wv)[2]){
      wv = t(wv)
      wv.rob = t(wv.rob)
      wv.mp = t(wv.mp)
    }
    
    WV.temp = wv[1:min(J1,J2),1:min(J1,J2)]
    WV.temp1 = (WV.temp[lower.tri(WV.temp,diag=TRUE)] + t(WV.temp)[t(upper.tri(WV.temp,diag=TRUE))])/2
    WV.temp[lower.tri(WV.temp,diag=TRUE)] = WV.temp1
    wv[1:min(J1,J2),1:min(J1,J2)] = WV.temp
    
    WV.temp = wv.rob[1:min(J1,J2),1:min(J1,J2)]
    WV.temp1 = (WV.temp[lower.tri(WV.temp,diag=TRUE)] + t(WV.temp)[t(upper.tri(WV.temp,diag=TRUE))])/2
    WV.temp[lower.tri(WV.temp,diag=TRUE)] = WV.temp1
    wv.rob[1:min(J1,J2),1:min(J1,J2)] = WV.temp
    
    WV.temp = wv.mp[1:min(J1,J2),1:min(J1,J2)]
    WV.temp1 = (WV.temp[lower.tri(WV.temp,diag=TRUE)] + t(WV.temp)[t(upper.tri(WV.temp,diag=TRUE))])/2
    WV.temp[lower.tri(WV.temp,diag=TRUE)] = WV.temp1
    wv.mp[1:min(J1,J2),1:min(J1,J2)] = WV.temp
    
    if(compute.v==TRUE){
      
      nstar = length(wv.coeffs[[(J1*J2)]])
      wv.coeffs = lapply(wv.coeffs, function(d) d[1:nstar])
      wv.coeffs = do.call("rbind",wv.coeffs)
      
      # Robust
      S.coeffs = matrix(NA,dim(wv.coeffs)[1],dim(wv.coeffs)[2])
      i = 0
      for(j1 in 1:J1){
        for(j2 in 1:J2){
          i = i+1
          S.coeffs[i,] = psi.tuk(wv.coeffs[i,],wv.rob[j1,j2])
        }
      }
      index = matrix(1:(J1*J2),J1,J2,byrow=T)
      S.low = S.coeffs[index[lower.tri(index,diag = T)],]
      S.up = S.coeffs[t(index)[lower.tri(index,diag = T)],]
      S.lower = tcrossprod(S.low[,1])
      S.upper = tcrossprod(S.up[,1])
      for(i in 2:nstar){
        S.lower = S.lower + tcrossprod(S.low[,i])
        S.upper = S.upper + tcrossprod(S.up[,i])
      }
      S = (S.lower + S.upper)/(2*nstar)
      M = matrix(0,(J1*J2),(J1*J2))
      for(k in 1:nstar){
        Mtemp = rep(NA,(J1*J2))
        i = 0
        for(j1 in 1:J1){
          for(j2 in 1:J2){
            i = i+1
            temp = der.psi.tuk(wv.coeffs[i,k],wv.rob[j1,j2])
            Mtemp[i] =  ifelse(is.na(temp),0,temp)
          }
        }
        M = M + diag(Mtemp)
      }
      M = (M[index[lower.tri(index,diag = T)],index[lower.tri(index,diag = T)]] + M[t(index)[lower.tri(index,diag = T)],t(index)[lower.tri(index,diag = T)]])/2
      M.temp = try(solve(-M/nstar),silent=T)
      if(class(M.temp)=="try-error"){
        M = ginv(-M/nstar)
      } else {
        M = M.temp
      }
      V.rob = M%*%S%*%t(M)/nstar
      
    }
    
    sig.temp = sig.ci[1:min(J1,J2),1:min(J1,J2)]
    sig.temp1 = (sig.temp[lower.tri(sig.temp,diag=TRUE)] + t(sig.temp)[t(upper.tri(sig.temp,diag=TRUE))])/2
    sig.temp[lower.tri(sig.temp,diag=TRUE)] = sig.temp1
    sig.ci[1:min(J1,J2),1:min(J1,J2)] = sig.temp
    sig.ci = sig.ci[lower.tri(sig.ci, diag=TRUE)]
    
    mn.temp = mn[1:min(J1,J2),1:min(J1,J2)]
    mn.temp1 = (mn.temp[lower.tri(mn.temp,diag=TRUE)] + t(mn.temp)[t(upper.tri(mn.temp,diag=TRUE))])/2
    mn.temp[lower.tri(mn.temp,diag=TRUE)] = mn.temp1
    mn[1:min(J1,J2),1:min(J1,J2)] = mn.temp
    mn = mn[lower.tri(mn, diag=TRUE)]
    
    wv = wv[lower.tri(wv,diag=TRUE)]
    wv.rob = wv.rob[lower.tri(wv.rob,diag=TRUE)]
    wv.mp = wv.mp[lower.tri(wv.mp,diag=TRUE)]
    
  } else {
    
    if(compute.v==TRUE){
      
      nstar = length(wv.coeffs[[nb.level]])
      wv.coeffs = lapply(wv.coeffs, function(d) d[1:nstar])
      wv.coeffs = do.call("rbind",wv.coeffs)
      
      # Robust
      S.coeffs = matrix(NA,dim(wv.coeffs)[1],dim(wv.coeffs)[2])
      i = 0
      for(j1 in 1:J1){
        for(j2 in 1:J2){
          i = i+1
          S.coeffs[i,] = psi.tuk(wv.coeffs[i,],wv.rob[j1,j2])
        }
      }
      S = tcrossprod(S.coeffs[,1])
      for(i in 2:nstar){
        S = S + tcrossprod(S.coeffs[,i])
      }
      S = S/nstar
      M = matrix(0,nb.level,nb.level)
      for(k in 1:nstar){
        Mtemp = rep(NA,nb.level)
        i = 0
        for(j1 in 1:J1){
          for(j2 in 1:J2){
            i = i+1
            temp = der.psi.tuk(wv.coeffs[i,k],wv.rob[j1,j2])
            Mtemp[i] =  ifelse(is.na(temp),0,temp)
          }
        }
        M = M + diag(Mtemp)
      }
      M.temp = try(solve(-M/nstar),silent=T)
      if(class(M.temp)=="try-error"){
        M = ginv(-M/nstar)
      } else {
        M = M.temp
      }
      V.rob = M%*%S%*%t(M)/nstar
      
    }
    
    wv = c(wv)
    wv.rob = c(wv.rob)
    wv.mp = c(wv.mp)
    sig.ci = c(sig.ci)
    mn = c(mn)
    
  }
  
  return(list(wv=wv,wv.rob=wv.rob,wv.mp=wv.mp,V=V,V.rob=V.rob,CI=CI,CI.rob=CI.rob,sig=sig.ci, mn=mn))
  
}