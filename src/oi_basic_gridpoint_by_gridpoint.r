#+ optimal interpolation 
oi_basic_gridpoint_by_gridpoint<-function(i,
                                          corr="soar",
                                          pmax,
                                          isolation_factor=7,
                                          dh=10000) {
# return(c(xa,xa_errvar,o_errvar,xidi,idiv,av))
# NOTE: av is the leave-one-out CV. However, its errvar is not returned.
#       todo: figure out how to compute the leave-ione-out errvar.
# global variables: xgrid_spint, ygrid_spint, yb_spint, xb_spint
#                   VecX, VecY, dh, dh2,
  if( i%%(round(env$m_dim/10)) == 0) cat(".")
  # select the observations to use
  if ( (p <- length( aux <- which(envtmp$nn2$nn.idx[i,]!=0))) == 0) return( c(envtmp$xb[i],NA))
  ixa  <- envtmp$nn2$nn.idx[i,aux]
  dist <- envtmp$nn2$nn.dists[i,aux]
  if (corr=="gaussian") {
    rloc <- exp( -0.5* (dist*dist) / (dh*dh) )
  } else if (corr=="soar")  {
    rloc <- (1+dist/dh)*exp(-dist/dh)
  } else if (corr=="powerlaw")  {
    rloc <- 1 / (1 + 0.5*(dist*dist)/(dh*dh))
  } else if (corr=="TOAR")  {
    rloc <- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
  }
#  deltax<-abs(envtmp$x[i]-y_env$yo$x)
#  deltay<-abs(envtmp$y[i]-y_env$yo$y)
#  if ( !any( deltax < (isolation_factor*dh))) return( c(envtmp$xb[i],NA))
#  if ( !any( deltay < (isolation_factor*dh))) return( c(envtmp$xb[i],NA))
#  ixa<-which( deltax<(isolation_factor*dh) & deltay<(isolation_factor*dh) )
#  if ( length(ixa) == 0) return( c(envtmp$xb[i],NA))
#  
#  if (corr=="gaussian") {
#    rloc<-exp( -0.5* (deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]) / dh2 )
#  } else if (corr=="soar")  {
#    distnorm<-sqrt(deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]) / dh
#    rloc<-(1+distnorm)*exp(-distnorm)
#    rm(distnorm)
#  }
#  if (length(ixa)>pmax) {
#    ixb<-order(rloc, decreasing=T)[1:pmax]
#    rloc<-rloc[ixb]
#    ixa<-ixa[ixb]
#    rm(ixb)
#  }
#  p <- length(ixa)
  x <- y_env$yo$x[ixa]
  y <- y_env$yo$y[ixa]
  yo <- y_env$yo$value[ixa]
  yb <- envtmp$yb[ixa]
  di <- yo - yb
  eps2 <- envtmp$eps2[i]

  if (corr=="gaussian") {
    rloc <- exp( -0.5* (dist*dist) / (dh*dh) )
  } else if (corr=="soar")  {
    rloc <- (1+dist/dh)*exp(-dist/dh)
  } else if (corr=="powerlaw")  {
    rloc <- 1 / (1 + 0.5*(dist*dist)/(dh*dh))
  } else if (corr=="TOAR")  {
    rloc <- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
  }

  if (corr=="gaussian") {
    S<-exp(-0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
  } else if (corr=="soar")  {
    distnorm<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2) / dh 
    S<-(1+distnorm)*exp(-distnorm)
    rm(distnorm)
  } else if (corr=="powerlaw")  {
    S<-1 / (1 + 0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
  } else if (corr=="TOAR")  {
    dist<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)
    S<- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
    rm(dist)
  }
  #
  SRinv<-chol2inv(chol( (S+diag(x=eps2,p)) ))
  xidi<-sum(rloc*as.vector(rowSums(SRinv)))
  SRinv_di<-crossprod(SRinv,di)       
#  o_errvar<-mean( di * (di-crossprod(S,SRinv_di)) )
  rm(S)
#  xa_errvar<-max(xa_errvar_min,(o_errvar/ mean(eps2[ixa]))) * 
#           (1-sum(as.vector(crossprod(rloc,SRinv))*rloc))
  xa<-envtmp$xb[i] + sum( rloc * as.vector(SRinv_di))
#  if (y_elab & !loocv) {
#    ii<-which(ixa==i)
#    Wii<-sum(rloc*SRinv[ii,])
#    idiv<-(xidi-Wii)/(1-Wii)
#    av<-(xa-Wii*yo_spint[i])/(1-Wii)
#  }
#  return(c(xa,xa_errvar,o_errvar,xidi,idiv,av))
  return(c(xa,xidi))
}

