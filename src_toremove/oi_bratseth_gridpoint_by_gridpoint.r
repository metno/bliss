# select the p_i observations nearest to the i-th gridpoint
nearest<-function(xx,yy,box_o_nearest_halfwidth,pmax) {
  deltax<-abs(xx-xobs_spint)
  if (!any(deltax<(box_o_nearest_halfwidth))) return(list(ix=NA,dist=NA))
  deltay<-abs(yy-yobs_spint)
  if (!any(deltay<(box_o_nearest_halfwidth))) return(list(ix=NA,dist=NA))
  ixa<-which( deltax<(box_o_nearest_halfwidth) & 
              deltay<(box_o_nearest_halfwidth) ) 
  disth2<-deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]
  if (length(ixa)>pmax) {
    ixb<-order(disth2, decreasing=F)[1:pmax]
    ixa<-ixa[ixb]
    disth2<-disth2[ixb]
    rm(ixb)
  }
  return(list(ix=ixa,dist=sqrt(disth2)))
}

#+ OI Bratset's iterative method
oi_bratseth_gridpoint_by_gridpoint<-function(i,
                                             nSCloops=10,
                                             delta_hor=10000, #m
                                             dh=10000, #m
                                             box_o_nearest_halfwidth=200000, #m
                                             dz=NA,
                                             lafmin=NA,
                                             dh_adaptive=F,
                                             dh_adaptive_min=0,
                                             dh_adaptive_max=10000000,
                                             corr="soar",
                                             pmax=200,
                                             loocv=F) {
#------------------------------------------------------------------------------
# global variables: xgrid_spint, ygrid_spint, zgrid_spint, lafgrid_spint,
#                   xobs_spint, yobs_spint, zobs_spint, lafobs_spint,
#                   yo_spint, yb_spint, xb_spint, eps2_spint
#                   nobs
# Description:
# Bratset's iterative method OI analysis at the i-th point (at xgrid_spint[i],
#  ygrid_spint[i],...) given: (a) the set of observations yo_spint 
#  (at xobs_spint, yobs_spint,...); (b) the background
# 
# Input arguments
# i: gridpoint index (refers to vectors Xgrid_spint,...,xb_spint)
# dh: horizontal de-correlation length (m)
# box_o_nearest_halfwidth: half-width of the square box used to select the 
#                          nearest observations
# dz: vertical de-correlation length (m, NA if z is not considered)
# lafmin: land-area fraction minimum value for the de-correlation factor 
#         (0-1, NA if laf is not considered)
# dh_adaptive: logical. F if dh is the actual horizontal de-correlation length
#                       T if the actual horizontal de-correlation length is 
#                         obtained as the 10-th percentile of the distances
#                         between the i-th gridpoint and the nearest observations
#                         (up to the nearest pmax observations). In this case dh
#                         is used as a lower limit.
# corr: model for the correlation functions. 
#       "soar" second order auto-regressive (only horizontal distance is used).
#       "gaussian" gaussian.
# pmax: mximum number of nearest observations to consider
#
# return(c(xa,dh))
#
#------------------------------------------------------------------------------
  deltax<-abs(xgrid_spint[i]-xobs_spint)
  deltay<-abs(ygrid_spint[i]-yobs_spint)
  if (!is.na(dz)) {deltaz<-abs(zgrid_spint[i]-zobs_spint);dz2<-dz*dz}
  if (!is.na(lafmin)) deltalaf<-abs(lafgrid_spint[i]-lafobs_spint)
  consider_obs<-rep(T,nobs)
  res<-ifelse(exists("xb_spint"),xb_spint[i],NA)
  av<-NA
  idiv<-NA
  xidi<-0
  if (loocv) consider_obs[which(deltax<1 & deltay<1)]<-F
  if (!any(deltax<(box_o_nearest_halfwidth))) return(c(res,dh))
  if (!any(deltay<(box_o_nearest_halfwidth))) return(c(res,dh))
  ixa<-which( deltax<(box_o_nearest_halfwidth) & 
              deltay<(box_o_nearest_halfwidth) & 
              consider_obs )
  if (length(ixa)==0) return(c(res,dh))
  disth2<-deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]
  if (length(ixa)>pmax) {
    ixb<-order(disth2, decreasing=F)[1:pmax]
    ixa<-ixa[ixb]
    disth2<-disth2[ixb]
    rm(ixb)
  }
  p_i<-length(ixa)
  yb_spint_i<-yb_spint[ixa]
  xb_i<-xb_spint[i]
  # correlation matrices
  if (dh_adaptive) {
    sw<-nearest(xgrid_spint[i]-delta_hor,
                ygrid_spint[i]-delta_hor,
                box_o_nearest_halfwidth,pmax) 
    se<-nearest(xgrid_spint[i]-delta_hor,
                ygrid_spint[i]+delta_hor,
                box_o_nearest_halfwidth,pmax) 
    ne<-nearest(xgrid_spint[i]+delta_hor,
                ygrid_spint[i]+delta_hor,
                box_o_nearest_halfwidth,pmax) 
    nw<-nearest(xgrid_spint[i]-delta_hor,
                ygrid_spint[i]+delta_hor,
                box_o_nearest_halfwidth,pmax) 
    dh<-min(dh_adaptive_max,
            max(dh_adaptive_min,
            mean(c(as.numeric(quantile(sqrt(disth2),probs=0.1,na.rm=T)),
                   as.numeric(quantile(sw$dist,probs=0.1,na.rm=T)),
                   as.numeric(quantile(se$dist,probs=0.1,na.rm=T)),
                   as.numeric(quantile(nw$dist,probs=0.1,na.rm=T)),
                   as.numeric(quantile(ne$dist,probs=0.1,na.rm=T))),na.rm=T)))
  }
  dh2<-dh*dh
  if (corr=="gaussian") {
    Gi<-exp( -0.5* disth2 / dh2 )
    S<-exp(-0.5*(outer(yobs_spint[ixa],yobs_spint[ixa],FUN="-")**2. + 
                 outer(xobs_spint[ixa],xobs_spint[ixa],FUN="-")**2)/dh2)
    # adjust gaussian correlations by taking into account more geo-parameters
    if (!is.na(dz)) {
      S<-S*exp(-0.5*abs(outer(zobs_spint[ixa],zobs_spint[ixa],FUN="-"))/dz2) 
      Gi<-Gi*exp(-0.5*deltaz[ixa]/dz2)
    }
    if (!is.na(lafmin)) {
      S<-S*(1-(1-lafmin)*
         abs(outer(lafobs_spint[ixa],lafobs_spint[ixa],FUN="-")))
      Gi<-Gi*(1-(1-lafmin)*deltalaf[ixa])
    }
  } else if (corr=="soar")  {
    distnorm_loc<-sqrt(disth2) / dh
    Gi<-(1+distnorm_loc)*exp(-distnorm_loc)
    rm(distnorm_loc)
    distnorm<-sqrt(outer(yobs_spint[ixa],yobs_spint[ixa],FUN="-")**2. + 
                   outer(xobs_spint[ixa],xobs_spint[ixa],FUN="-")**2) / dh 
    S<-(1+distnorm)*exp(-distnorm)
    rm(distnorm)
  }
  SR<-S+diag(x=eps2[ixa],length(ixa))
  Minv<-1/rowSums(abs(SR))
  I_SR_Minv<-diag(length(ixa))-t( t(SR) * Minv )
  d0<-yo_spint[ixa]-yb_spint[ixa]
  di<-yo_spint[ixa]-yb_spint[ixa]
  for (sc in 1:nSCloops)  di<-as.vector(crossprod(t(I_SR_Minv),as.vector(di)))
  di<-d0+di
  xa<-xb_spint[i]+sum(Gi*(Minv*di))
  return(c(xa,dh))
}
