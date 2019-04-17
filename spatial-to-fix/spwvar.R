#' @title Wavelet Variance for Regular Lattice Spatial Data
#' 
#' @description
#' General function to compute standard and robust WV for regular lattice spatial data (isotropic and anisotropic)
#' @param X            A \code{matrix} with dimensions N x M
#' @param J1           TO DO    
#' @param J2           TO DO   
#' @param eff          TO DO   
#' @param compute.v    TO DO 
#' @param iso          TO DO 
#' @param ...       Further arguments passed to or from other methods.
#' @return A \code{list} with the structure: IT MIGHT BE A GOOD TO MATCH THE OUTPUT OF WVAR
#' \describe{
#'   \item{"wv"}{TO DO}
#'   \item{"wv.rob"}{TO DO}
#'   \item{"wv.mp"}{TO DO}
#'   \item{"V"}{TO DO}
#' }
#' @details 
#' TO DO
#' @author Roberto Molinari
#' @examples
#' # An example TO DO
spat_wavar = function(X, J1 = floor(log2(dim(X)[1])), J2 = floor(log2(dim(X)[2])), eff=0.6, compute.v = FALSE, iso = TRUE, ...){
  
  n = dim(X)[1]
  m = dim(X)[2]
  wv = matrix(NA, J1, J2)
  wv.rob = matrix(NA, J1, J2)
  wv.mp = matrix(NA, J1, J2)
  sig.ci = matrix(NA, J1, J2)
  mn = matrix(NA, J1, J2)
  CI = matrix(NA, min(J1, J2), 3)
  CI.rob = matrix(NA, min(J1, J2), 3)
  nb.level = J1*J2
  
  if(compute.v == TRUE) wv.coeffs = vector("list", nb.level)
  
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
      
      if(compute.v == TRUE) wv.coeffs[[i]] = wv.coeff
      
      # Standard
      wv[j1,j2] = sum(wv.coeff^2,na.rm=T)/(pp1*pp2)
      
      # Robust
      wv.rob[j1,j2] = sig_rob_bw(wv.coeff, eff)
      
    }
    
  }
  
  V = NA
  
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
    
    # Compute the covariance matrix of the WV
    if(compute.v==TRUE){
      
      nstar = length(wv.coeffs[[(J1*J2)]])
      wv.coeffs = lapply(wv.coeffs, function(d) d[1:nstar])
      wv.coeffs = do.call("rbind",wv.coeffs)
      
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
      V = M%*%S%*%t(M)/nstar
      
    }
    
    wv = wv[lower.tri(wv,diag=TRUE)]
    wv.rob = wv.rob[lower.tri(wv.rob,diag=TRUE)]
    wv.mp = wv.mp[lower.tri(wv.mp,diag=TRUE)]
    
  } else {
    
    # Compute the covariance matrix of the WV
    if(compute.v==TRUE){
      
      nstar = length(wv.coeffs[[nb.level]])
      wv.coeffs = lapply(wv.coeffs, function(d) d[1:nstar])
      wv.coeffs = do.call("rbind",wv.coeffs)
      
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
      V = M%*%S%*%t(M)/nstar
      
    }
    
    wv = c(wv)
    wv.rob = c(wv.rob)
    wv.mp = c(wv.mp)
    
  }
  
  return(list(wv=wv,wv.rob=wv.rob,wv.mp=wv.mp,V=V))
  
}