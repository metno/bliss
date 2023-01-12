#+ vertical profile of temperature (Frei, 2014)
tvertprof<-function(z,t0,gamma,a,h0,h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region 
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
#  a= numeric. inversion (spatial) length
#  h0= numeric. z where inversion starts [m]
#  h1i= numeric. h0+h1i is z where inversion stops [m]
#       (Frei uses h1 directly, I use an increment to h0 so to avoid ending
#        up with h1<=h0 during the optimization)
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  t<-z
  t[]<-NA
  h1<-h0+abs(h1i)
  z.le.h0<-which(z<=h0)
  z.ge.h1<-which(z>=h1)
  z.in<-which(z>h0 & z<h1)
  if (length(z.le.h0)>0)
   t[z.le.h0]<-t0-gamma*z[z.le.h0]-a 
  if (length(z.ge.h1)>0)
   t[z.ge.h1]<-t0-gamma*z[z.ge.h1] 
  if (length(z.in)>0)
   t[z.in]<-t0-gamma*z[z.in]-a/2*(1+cos(pi*(z[z.in]-h0)/(h1-h0)))
  return(t)
}

#+ cost function used for optimization of tvertprof parameter
tvertprof2opt<-function(par) {
  te<-tvertprof(z=zopt,t0=par[1],gamma=par[2],a=par[3],h0=par[4],h1i=par[5])
  return(log((mean((te-topt)**2))**0.5))
}

#+ vertical profile of temperature (linear)
tvertprof_basic<-function(z,t0,gamma) {
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  return(t0+gamma*z)
}

#+ cost function used for optimization of tvertprof parameter
tvertprofbasic2opt<-function(par) {
  te<-tvertprof_basic(z=zopt,t0=par[1],gamma=argv$gamma.standard)
  return(log((mean((te-topt)**2))**0.5))
}

#+ cost function used for optimization of vertprof parameter
vertprofbasic2opt<-function(par,vert_coord,gamma,obs) {
  pred<-tvertprof_basic(z=vert_coord,t0=par[1],gamma=gamma)
  return(log((mean((pred-obs)**2))**0.5))
}

#+ cost function used for optimization of vertprof parameter
vertprof2opt<-function(par,vert_coord,obs) {
  pred<-tvertprof(z=vert_coord,
                  t0=par[1],
                  gamma=par[2],
                  a=par[3],
                  h0=par[4],
                  h1i=par[5])
  return(log((mean((pred-obs)**2))**0.5))
}

#------------------------------------------------------------------------------
# optimal interpolation 
oi_var_gridpoint_by_gridpoint<-function(i,
                                        dh=10000, #m
                                        box_o_nearest_halfwidth=100000, #m
                                        dz=NA,
                                        lafmin=NA,
                                        dh_adaptive=F,
                                        corr="soar",
                                        pmax,
                                        fg=NA,
                                        fg_gamma=NA,
                                        fg_min=NA,
                                        fg_max=NA,
                                        return_fg_only=F,
                                        succ_corr=F,
                                        y_elab=F,
                                        loocv=F,
                                        o_errvar_min=0.001,
                                        o_errvar_max=4,
                                        xa_errvar_min=0.001,
                                        xa_errvar_max=4) {
# global variables: xgrid_spint, ygrid_spint, zgrid_spint, lafgrid_spint,
#                   xobs_spint, yobs_spint, zobs_spint, lafobs_spint,
#                   yo_spint, yb_spint, xb_spint, eps2_spint
#                   nobs
# Description:
# OI analysis at the i-th point (at xgrid_spint[i],ygrid_spint[i],...)
# given the set of observations yo_spint (at xobs_spint,yobs_spint,...)
# with or without a backgorund
# 
# Input arguments
# i: gridpoint index (refers to vectors Xgrid_spint,...,xb_spint
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
# fg: method used to compute the background (first-guess). 
#  NA (default), background available in yb_spint, xb_spint
#  "linear", background as a linear function of the vertical coordinate  
#  "Frei", background as a non-linear function of the vertical coordinate  
#  "mean", background set to the mean value of the nearest observations
# fg_gamma: vertical lapse rate used for fg="linear" and optimized for fg="Frei"
# fg_min: minimum allowed value for the background
# fg_max: maximum allowed value for the background
# succ_corr: logical.
#  F (default). the background is used as it is.
#  T. three steps of successive corrections are applied to the background.
# y_elab: logical.
#  F (default). grid.. vectors and obs.. vectors refer to the same locations
#  T. grid.. and obs.. vectors refer to different locations
# loocv. logical. Leave one out cross-validation
#  F (default). standard run, use all the observations
#  T. run in loocv mode, without considering the observations at the i-th gridpoint
# o_errvar_min. minimum allowed value for the observation error variance.
# o_errvar_max. maximum allowed value for the observation error variance.
# xa_errvar_min. minimum allowed value for the analysis error variance.
# xa_errvar_max. maximum allowed value for the analysis error variance.
#
# return(c(xa,xa_errvar,o_errvar,xidi,idiv,av,dh))
# NOTE: av is the leave-one-out CV. However, its errvar is not returned.
#       todo: figure out how to compute the leave-ione-out errvar.
#------------------------------------------------------------------------------
# select the p_i observations nearest to the i-th gridpoint
  deltax<-abs(xgrid_spint[i]-xobs_spint)
  deltay<-abs(ygrid_spint[i]-yobs_spint)
  if (!is.na(dz)) {deltaz<-abs(zgrid_spint[i]-zobs_spint);dz2<-dz*dz}
  if (!is.na(lafmin)) deltalaf<-abs(lafgrid_spint[i]-lafobs_spint)
  consider_obs<-rep(T,nobs)
  if (y_elab) {
    res<-ifelse(exists("yb_spint"),yb_spint[i],NA)
    av<-NA
    idiv<-0
    xidi<-1/(1+eps2_spint[i])
    if (loocv) consider_obs[i]<-F
    if (length(which(deltax<(box_o_nearest_halfwidth)))==1) { 
      if (return_fg_only) {return(NA)} 
      else { return(c(res,NA,NA,xidi,idiv,av,dh))}
    }
    if (length(which(deltay<(box_o_nearest_halfwidth)))==1) {
        if (return_fg_only) {return(NA)} 
        else { return(c(res,NA,NA,xidi,idiv,av,dh))}
    }
  } else {
    res<-ifelse(exists("xb_spint"),xb_spint[i],NA)
    av<-NA
    idiv<-NA
    xidi<-0
    if (loocv) consider_obs[which(deltax<1 & deltay<1)]<-F
    if (!any(deltax<(box_o_nearest_halfwidth))) {
      if (return_fg_only) {return(NA)} 
      else { return(c(res,NA,NA,xidi,idiv,av,dh))}
    }
    if (!any(deltay<(box_o_nearest_halfwidth))) {
      if (return_fg_only) {return(NA)} 
      else { return(c(res,NA,NA,xidi,idiv,av,dh))}
    }
  }
  ixa<-which( deltax<(box_o_nearest_halfwidth) & 
              deltay<(box_o_nearest_halfwidth) & 
              consider_obs )
  if (length(ixa)==0) 
    if (return_fg_only) {return(NA)} 
    else { return(c(res,NA,NA,xidi,idiv,av,dh))}
  disth2<-deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]
  if (length(ixa)>pmax) {
    ixb<-order(disth2, decreasing=F)[1:pmax]
    ixa<-ixa[ixb]
    disth2<-disth2[ixb]
    rm(ixb)
  }
  p_i<-length(ixa)
  # first-guess from nearest observations
  if (!is.na(fg)) {
    yo_mean<-mean(yo_spint[ixa])
    # Frei profile may provide unrealistic values if 
    #  it is required to extrapole a value for elevations 
    #  far from the ones used to optimize the parameters
    #  OR if the elevations used are within a narrow layer
    if (fg=="Frei") {
      zmin<-sort(zobs_spint[ixa])[min(c(length(ixa),2))]
      zmax<-sort(zobs_spint[ixa])[max(c(1,(p_i-1)))]
      if ( (zmin-zgrid_spint[i])>50 |
           (zgrid_spint[i]-zmax)>50 |
           (zmax-zmin)<25 ) fg<-"linear"
    }
    if (fg=="linear") {
      par<-c(yo_mean)
      opt<-optimize(f=vertprofbasic2opt,
                    interval=c(fg_min,fg_max),
                    vert_coord=zobs_spint[ixa],
                    gamma=fg_gamma,
                    obs=yo_spint[ixa])
      yb_spint_i<-tvertprof_basic(zobs_spint[ixa],
                                  t0=opt$minimum,
                                  gamma=fg_gamma)
      xb_i<-tvertprof_basic(zgrid_spint[i],
                            t0=opt$minimum,
                            gamma=fg_gamma)
    } else if (fg=="Frei") {
      par<-c(yo_mean,
             fg_gamma,
             5,
             zmin,
             zmax)
      opt<-optim(par,vertprof2opt,vert_coord=zobs_spint[ixa],obs=yo_spint[ixa])
      yb_spint_i<-tvertprof(z=zobs_spint[ixa],
                            t0=opt$par[1],
                            gamma=opt$par[2],
                            a=opt$par[3],
                            h0=opt$par[4],
                            h1i=opt$par[5])
      xb_i<-tvertprof(z=zgrid_spint[i],
                      t0=opt$par[1],
                      gamma=opt$par[2],
                      a=opt$par[3],
                      h0=opt$par[4],
                      h1i=opt$par[5])
    } else if (fg=="mean") {
      yb_spint_i<-rep(yo_mean,length=p_i)
      xb_i<-yo_mean
    }
  } else {
    yb_spint_i<-yb_spint[ixa]
    xb_i<-xb_spint[i]
  } # end compute first-guess
  if (!is.null(fg_min)) {
    if (!is.na(fg_min) & !is.nan(fg_min)) {
      yb_spint_i[which(yb_spint_i<fg_min)]<-fg_min
      xb_i<-max(xb_i,fg_min)
    }
  }
  if (!is.null(fg_max)) {
    if (!is.na(fg_max) & !is.nan(fg_max)) {
      yb_spint_i[which(yb_spint_i>fg_max)]<-fg_max
      xb_i<-min(xb_i,fg_max)
    }
  }
  if (return_fg_only) return(xb_i)
  # correlation matrices
  if (dh_adaptive) {
    dh<-max(dh,as.numeric(quantile(sqrt(disth2),probs=0.1)))
  }
  dh2<-dh*dh
  if (corr=="gaussian") {
    rloc<-exp( -0.5* disth2 / dh2 )
  } else if (corr=="soar")  {
    distnorm_loc<-sqrt(disth2) / dh
    rloc<-(1+distnorm_loc)*exp(-distnorm_loc)
    if (!succ_corr) rm(distnorm_loc)
  }
  if (corr=="gaussian") {
    S<-exp(-0.5*(outer(yobs_spint[ixa],yobs_spint[ixa],FUN="-")**2. + 
                 outer(xobs_spint[ixa],xobs_spint[ixa],FUN="-")**2)/dh2)
  } else if (corr=="soar")  {
    distnorm<-sqrt(outer(yobs_spint[ixa],yobs_spint[ixa],FUN="-")**2. + 
                   outer(xobs_spint[ixa],xobs_spint[ixa],FUN="-")**2) / dh 
    S<-(1+distnorm)*exp(-distnorm)
    if (!succ_corr) rm(distnorm)
  }
  # successive corrections (Barnes scheme) step
  if (succ_corr) {
    for (sc in 3:1) {
      if (corr=="gaussian") {
        S1<-S**(1/sc**2)
        rloc1<-rloc**(1/sc**2)
      } else if (corr=="soar")  {
        distnorm1<-distnorm/sc
        S1<-(1+distnorm1)*exp(-distnorm1)
        distnorm1_loc<-distnorm_loc/sc
        rloc1<-(1+distnorm1_loc)*exp(-distnorm1_loc)
      }
      yb_spint_i<-yb_spint_i+crossprod(S1,(yo_spint[ixa]-yb_spint_i)) / 
                  (rowSums(S1)+eps2_spint[ixa])
      xb_i<-xb_i+sum(rloc1*(yo_spint[ixa]-yb_spint_i))/sum(rloc1+eps2_spint[ixa])
    }
    rm(S1,rloc1)
    if (corr=="soar") rm(distnorm1,distnorm,distnorm_loc,distnorm1_loc)
  }
  # adjust gaussian correlations by taking into account more geo-parameters
  if (corr=="gaussian") {
    if (!is.na(dz)) {
      S<-S*exp(-0.5*abs(outer(zobs_spint[ixa],zobs_spint[ixa],FUN="-"))/dz2) 
      rloc<-rloc*exp(-0.5*deltaz[ixa]/dz2)
    }
    if (!is.na(lafmin)) {
      S<-S*(1-(1-lafmin)*
         abs(outer(lafobs_spint[ixa],lafobs_spint[ixa],FUN="-")))
      rloc<-rloc*(1-(1-lafmin)*deltalaf[ixa])
    }
  }
  # innovation
  di<-yo_spint[ixa]-yb_spint_i
  # OI analysis
  SRinv<-chol2inv(chol( (S+diag(x=eps2_spint[ixa],length(ixa))) ))
  xidi<-sum(rloc*as.vector(rowSums(SRinv)))
  SRinv_di<-crossprod(SRinv,di)       
  o_errvar<-min(c(o_errvar_max,
                  max(c(o_errvar_min,
                        mean(di*(di-crossprod(S,SRinv_di)))))))
  rm(S)
  xa_errvar<-min(c(xa_errvar_max,
                   max(c(xa_errvar_min,
                         (o_errvar/ mean(eps2_spint[ixa])) * 
                         (1-sum(as.vector(crossprod(rloc,SRinv))*rloc))))))
  xa<-xb_i+sum(rloc*as.vector(SRinv_di))
  if (y_elab & !loocv) {
    ii<-which(ixa==i)
    Wii<-sum(rloc*SRinv[ii,])
    idiv<-(xidi-Wii)/(1-Wii)
    av<-(xa-Wii*yo_spint[i])/(1-Wii)
  }
  # debug
  #if (yo_to_check[i]>20) {
  #png(file=paste0("png/tvert_",i,".png"),width=800,height=800)
  #plot(aaa,1:2000,ylim=range(c(zobs_spint[ixa],zgrid_spint[i])),xlim=range(c(yo_spint[ixa],yb_spint_i,xb_i,yo_to_check[i])),
  #main=paste(round(xa,1),round(xa_errvar,3),round(o_errvar,3),round(xidi,3),round(dh,1),round((yo_to_check[i]-xa)**2/(xa_errvar+o_errvar),2)))
  #points(yo_spint[ixa],zobs_spint[ixa],pch=21,bg="blue")
  #points(yb_spint_i,zobs_spint[ixa],pch=21,bg="cyan")
  #points(xb_i,zgrid_spint[i],pch=21,bg="red")
  #points(yo_to_check[i],zgrid_spint[i],pch=21,bg="red")
  #dev.off()
  #}
  return(c(xa,xa_errvar,o_errvar,xidi,idiv,av,dh))
}
