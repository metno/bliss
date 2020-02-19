#+ optimal interpolation 
oi_var_gridpoint_by_gridpoint<-function(i,
                                        corr="soar",
                                        pmax,
                                        y_elab=F,
                                        loocv=F,
                                        isolation_factor=7,
                                        xa_errvar_min=0.001) {
# return(c(xa,xa_errvar,o_errvar,xidi,idiv,av))
# NOTE: av is the leave-one-out CV. However, its errvar is not returned.
#       todo: figure out how to compute the leave-ione-out errvar.
# global variables: xgrid_spint, ygrid_spint, yb_spint, xb_spint
#                   VecX, VecY, dh, dh2,
  deltax<-abs(xgrid_spint[i]-VecX)
  deltay<-abs(ygrid_spint[i]-VecY)
  if (y_elab) {
    av<-NA
    idiv<-0
    xidi<-1/(1+eps2[i])
    if (length(which(deltax<(isolation_factor*dh)))==1) 
      return(c(yb_spint[i],NA,NA,xidi,idiv,av))
    if (length(which(deltay<(isolation_factor*dh)))==1) 
      return(c(yb_spint[i],NA,NA,xidi,idiv,av))
    exclude_i<-rep(T,n0)
    if (loocv) exclude_i[i]<-F
    ixa<-which( deltax<(isolation_factor*dh) & 
                deltay<(isolation_factor*dh) & 
                exclude_i )
    if (length(ixa)==0) return(c(yb_spint[i],NA,NA,xidi,idiv,av))
  } else {
    av<-NA
    idiv<-NA
    if (!any(deltax<(isolation_factor*dh))) return(c(xb_spint[i],NA,NA,0,idiv,av))
    if (!any(deltay<(isolation_factor*dh))) return(c(xb_spint[i],NA,NA,0,idiv,av))
    ixa<-which( deltax<(isolation_factor*dh) & 
                deltay<(isolation_factor*dh) )
    if (length(ixa)==0) return(c(xb_spint[i],NA,NA,0,idiv,av))
  }
  if (corr=="gaussian") {
    rloc<-exp( -0.5* (deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]) / dh2 )
  } else if (corr=="soar")  {
    distnorm<-sqrt(deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]) / dh
    rloc<-(1+distnorm)*exp(-distnorm)
    rm(distnorm)
  }
  if (length(ixa)>pmax) {
    ixb<-order(rloc, decreasing=T)[1:pmax]
    rloc<-rloc[ixb]
    ixa<-ixa[ixb]
    rm(ixb)
  }
  di<-yo_spint[ixa]-yb_spint[ixa]
  if (corr=="gaussian") {
    S<-exp(-0.5*(outer(VecY[ixa],VecY[ixa],FUN="-")**2. + 
                 outer(VecX[ixa],VecX[ixa],FUN="-")**2)/dh2)
  } else if (corr=="soar")  {
    distnorm<-sqrt(outer(VecY[ixa],VecY[ixa],FUN="-")**2. + 
                   outer(VecX[ixa],VecX[ixa],FUN="-")**2) / dh 
    S<-(1+distnorm)*exp(-distnorm)
    rm(distnorm)
  }
  #
  SRinv<-chol2inv(chol( (S+diag(x=eps2[ixa],length(ixa))) ))
  xidi<-sum(rloc*as.vector(rowSums(SRinv)))
  SRinv_di<-crossprod(SRinv,di)       
  o_errvar<-mean( di * (di-crossprod(S,SRinv_di)) )
  rm(S)
  xa_errvar<-max(xa_errvar_min,(o_errvar/ mean(eps2[ixa]))) * 
           (1-sum(as.vector(crossprod(rloc,SRinv))*rloc))
  xa<-xb_spint[i]+sum(rloc*as.vector(SRinv_di))
  if (y_elab & !loocv) {
    ii<-which(ixa==i)
    Wii<-sum(rloc*SRinv[ii,])
    idiv<-(xidi-Wii)/(1-Wii)
    av<-(xa-Wii*yo_spint[i])/(1-Wii)
  }
  return(c(xa,xa_errvar,o_errvar,xidi,idiv,av))
}

