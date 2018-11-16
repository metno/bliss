#!/usr/bin/env Rscript
# --~- bliss.r -~--
# Bayesian statisticaL Interpolation for Spatial analySis
# See the software repository here: 
#..............................................................................
#Copyright and license
# Copyright (C) 2018 MET Norway. The software is licensed under GPL version 3 
# or (at your option) any later version.
# https://www.gnu.org/licenses/gpl-3.0.en.html
# 
# History:
# 26.10.2018 - Cristian Lussana. Original code.
# -----------------------------------------------------------------------------
#
rm(list=ls())
#
# -----------------------------------------------------------------------------
# Libraries
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sp"))
suppressPackageStartupMessages(library("raster"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("rgdal"))
suppressPackageStartupMessages(library("ncdf4"))
suppressPackageStartupMessages(library("dotnc"))
options(warn = 2, scipen = 999)
#options(scipen = 999)
# 
# -----------------------------------------------------------------------------
# Constants
# CRS strings (always useful for copy-and-paste)
proj4.llwgs84<-"+proj=longlat +datum=WGS84"
#proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
#proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#proj4.aeqd<-"+proj=aeqd +lat_0=59.458925 +lon_0=10.564472 +R=6371000 +datum=WGS84"
#
#..............................................................................
# Functions
# + manage fatal error
boom<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}

#+ check if two rasters match
rasters_match<-function(r1,r2) {
  return ( projection(r1)==projection(r2) &
           extent(r1)==extent(r2) &
           res(r1)[1]==res(r2)[1] & 
           res(r1)[2]==res(r2)[2])
}

#+
set_NAs_to_NULL<-function(x) {
  if (!is.null(x)) {
    if (is.na(x)) x<-NULL
  }
  x
}

#
debug_plots<-function() {
  png(file=paste0("rb_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  image(rb,breaks=c(0,seq(0.001,0.8,length=25)),
        col=c("gray",rev(rainbow(24))))
  points(VecX,VecY,pch=19,col="black",cex=0.5)
  dev.off()
  png(file=paste0("innoh_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  d<-yo_relan[ixwet]-yb[ixwet]
  hist(d,breaks=seq(-1000.025,1000,by=.05),xlim=c(-1,1))
  dev.off()
  png(file=paste0("plot_",formatC(l,width=2,flag="0"),".png"),
      width=800,height=800)
  plot(yo_relan,yb,pch=19,col="blue",xlim=c(0,1),ylim=c(0,1))
  ix9<-which(prId==9)
  points(yo_relan[ix9],yb[ix9],pch=19,col="red")
  lines(-1000:1000,-1000:1000,col="red",lwd=2)
  dev.off()
}

#
debug_01_plots<-function() {
  png(file=paste0("ra.png"),width=800,height=800)
  image(ra,breaks=c(0,seq(0.1,50,length=45)),col=c("gray",rev(rainbow(44))))
  dev.off()
  png(file=paste0("ra1.png"),width=800,height=800)
  image(ra,breaks=c(0,0.09999,1000),col=c("gray","blue"))
  points(VecX[ixwet],VecY[ixwet],pch=19,col="cyan",cex=0.5)
  points(VecX[ixdry],VecY[ixdry],pch=19,col="red",cex=0.5)
  dev.off()
}

#+ find the n-th largest element from each matrix row 
findRow <- function (x,n) {   
# order by row number then by value
  y<-t(x)
  array(y[order(col(y), y)], dim(y))[nrow(y) - (n-1), ]
}

#+ mean radial distance between an observation and its k-th closest obs
dobs_fun<-function(obs,k) {
# NOTE: k=1 will return 0 everywhere (1st closest obs is the obs itself)
  nobs<-length(obs$x)
  if (nobs<k) return(NA)
  # distance matrices
  disth<-(outer(obs$x,obs$x,FUN="-")**2.+
          outer(obs$y,obs$y,FUN="-")**2.)**0.5
  dobsRow<-findRow(x=disth,n=(nobs-k+1))
  mean(dobsRow)
}

# + replace elements of a string with date-time elements
`replaceDate`<-function(string=NULL,
                        date.str=NULL,
                        format="%Y-%m-%d %H:%M:%S") {
#------------------------------------------------------------------------------
  if (is.null(string) | is.null(date.str)) return(NULL)
  Rdate<-as.POSIXlt(str2Rdate(date.str,format=format))
  yyyy<-Rdate$year+1900
  mm<-formatC(Rdate$mon+1,width=2,flag="0")
  dd<-formatC(Rdate$mday,width=2,flag="0")
  hh<-formatC(Rdate$hour,width=2,flag="0")
  out<-gsub("yyyy",yyyy,string)
  out<-gsub("mm",formatC(mm,width=2,flag="0"),out)
  out<-gsub("dd",formatC(dd,width=2,flag="0"),out)
  out<-gsub("hh",formatC(hh,width=2,flag="0"),out)
  out
}

#------------------------------------------------------------------------------
# code specific for RR
#+
`OI_RR_fast`<-function(yo,
                       yb,
                       xb,
                       xgrid,
                       ygrid,
                       zgrid=NULL,
                       VecX,
                       VecY,
                       VecZ=NULL,
                       Dh,
                       Dz=NULL) {
#------------------------------------------------------------------------------
  if (is.null(Dz)) {
    Dz<-100000
    VecZ<-rep(0,length(VecX))
    zgrid<-rep(0,length(xgrid))
  }
  no<-length(yo)
  ng<-length(xb)
  xa<-vector(mode="numeric",length=ng)
  vec<-vector(mode="numeric",length=no)
  d<-yo-yb
  out<-.C("oi_rr_first",no=as.integer(no), 
                        innov=as.double(d),
                        SRinv=as.numeric(InvD),
                        vec=as.double(vec) ) 
  vec[1:no]<-out$vec[1:no]
  rm(out)
  out<-.C("oi_rr_fast",ng=as.integer(ng),
                       no=as.integer(no),
                       xg=as.double(xgrid),
                       yg=as.double(ygrid),
                       zg=as.double(zgrid),
                       xo=as.double(VecX),
                       yo=as.double(VecY),
                       zo=as.double(VecZ),
                       Dh=as.double(Dh),
                       Dz=as.double(Dz),
                       xb=as.double(xb),
                       vec=as.double(vec),
                       xa=as.double(xa) )
  out$xa[1:ng]
}

#+
`OI_RR_var`<-function(yo,
                      yb,
                      xb,
                      gx,
                      gy,
                      ox,
                      oy,
                      Dh,
                      eps2
                      ) {
#------------------------------------------------------------------------------
  no<-length(yo)
  ng<-length(gx)
  xa<-vector(mode="numeric",length=ng)
  ya<-vector(mode="numeric",length=no)
  xa_errvar<-vector(mode="numeric",length=ng)
  ya_errvar<-vector(mode="numeric",length=no)
  o_errvar<-0
  xa[]<-0
  ya[]<-0
  xa_errvar[]<-0
  ya_errvar[]<-0
  out<-.C("oi_rr_var",ng=as.integer(ng),
                      no=as.integer(no),
                      SRinv=as.numeric(InvD),
                      eps2=as.double(eps2),
                      Dh=as.double(Dh),
                      gx=as.double(gx),
                      gy=as.double(gy),
                      ox=as.double(ox),
                      oy=as.double(oy),
                      yo=as.double(yo),
                      yb=as.double(yb),
                      xb=as.double(xb),
                      xa=as.double(xa),
                      ya=as.double(ya),
                      xa_errvar=as.double(xa_errvar),
                      ya_errvar=as.double(ya_errvar),
                      o_errvar=as.double(o_errvar))
  return( list(xa= out$xa[1:ng],
               ya= out$ya[1:no],
               xa_errvar= out$xa_errvar[1:ng],
               ya_errvar= out$ya_errvar[1:no],
               o_errvar= out$o_errvar))
}

#
boxcox<-function(x,lambda) {
  if (lambda==0) {
    return(log(x))
  } else {
    return((x**lambda-1)/lambda)
  }
}

#+
tboxcox4pdf_apply<-function(par,
                            lambda,
                            brrinf,
                            int=400) {
  mean<-par[1]
  sd<-par[2]
  if (mean<brrinf) return(0)
  if (sd==0) return(tboxcox(mean,lambda))
  aux<-qnorm(mean=mean,sd=sd,p=seq(1/int,(1-1/int),length=(int-1)))
  if (lambda!=0) aux[aux<(-1/lambda)]<-(-1/lambda)
  mean(tboxcox(aux,lambda=lambda))
}

#+
tboxcox<-function(x,lambda) {
  if (lambda==0) {
    return(exp(x))
  } else {
    return((1+lambda*x)**(1./lambda))
  }
}

#------------------------------------------------------------------------------
# code specific for Temperature
#+ background calculations. OI is used to compute the weights
`OI_T_xb_upd`<-function(xg,
                        yg,
                        zg,
                        bg,
                        wg,
                        dg,
                        dg_new,
                        xo,
                        yo,
                        vec1,
                        tpar,
                        na,
                        Dh) {
#------------------------------------------------------------------------------
  no<-length(xo)
  ng<-length(xg)
  out<-.C("oi_t_xb_upd",ng=as.integer(ng),
                        no=as.integer(no),
                        xg=as.double(xg),
                        yg=as.double(yg),
                        zg=as.double(zg),
                        bg=as.double(bg),
                        wg=as.double(wg),
                        dg=as.double(dg),
                        dg_new=as.double(dg_new),
                        xo=as.double(xo),
                        yo=as.double(yo),
                        Dh=as.double(Dh),
                        vec1=as.double(vec1),
                        tpar=as.double(tpar),
                        na=as.double(na) )
  return(list(bg_up=out$bg[1:ng],
              wg_up=out$wg[1:ng],
              dg_up=out$dg[1:ng]))
}



#+ vertical profile of temperature (Frei, 2014)
tvertprof_Frei<-function(z,t0,gamma,a,h0,h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region 
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
#  a= numeric. Temperature contrast or inversion strength [degC] 
#  h0= numeric. z where inversion starts [m]
#  h1i= numeric. h0+h1i is z where inversion stops [m]
#       (Frei uses h1 directly, I use an increment to h0 so to avoid ending
#        up with h1<=h0 during the optimization)
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  t<-z
  t[]<-NA
  h1<-h0+h1i
  z.le.h0<-which(z<=h0)
  z.ge.h1<-which(z>=h1)
  z.in<-which(z>h0 & z<h1)
  if (length(z.le.h0)>0)
   t[z.le.h0]<-t0+gamma*z[z.le.h0]-a 
  if (length(z.ge.h1)>0)
   t[z.ge.h1]<-t0+gamma*z[z.ge.h1] 
  if (length(z.in)>0)
   t[z.in]<-t0+gamma*z[z.in]-a/2*(1+cos(pi*(z[z.in]-h0)/h1i))
  return(t)
}

#+ cost function used for optimization of tvertprof parameter
tvertprofFrei2opt<-function(par) {
  te<-tvertprof_Frei(z=zopt,t0=par[1],gamma=par[2],a=par[3],h0=par[4],h1i=par[5])
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

#+
background_incAv<-function(ixynp, # sub-region centroids
                           maxboxl=250000,
                           nmin=50,
                           refdist=50000,
                           dzmin=30,
                           eps2=0.25,
                           closeNth=4
                           ) {
#-------------------------------------------------------------------------------
  # define the extension of the rectangular box used as sub-region (local scale) 
  n<-0
  fact<-1.9
  while (n<nmin) {
    fact<-fact+0.1
    ix<-which( abs(xr[ixynp[1]]-VecX_bg)<(fact*res[1]) & 
               abs(yr[ixynp[1]]-VecY_bg)<(fact*res[2]) )
    n<-length(ix)
    if (mean(fact*res)>maxboxl) break
  }
  if (n<nmin) return()
  assign("zopt",VecZ_bg[ix],envir=.GlobalEnv)
  dz<-as.numeric(quantile(zopt,probs=0.95))-
      as.numeric(quantile(zopt,probs=0.05))
  assign("topt",yo_bg[ix],envir=.GlobalEnv)
  # if region is too flat or not enough obs, then go for a basic profile
  tpar<-vector(length=6,mode="numeric")
  tpar[]<-NA
  if (dz<dzmin) {
    par<-c(mean(topt))
    opt<-optimize(f=tvertprofbasic2opt,interval=c(argv$vmin,argv$vmax))
    tpar<-c(1,opt$minimum,argv$gamma.standard,na,na,na)
  # otherwise, fit a vertical profile of temperature
  } else {
    ui<-matrix(ncol=5,nrow=10,data=0)
    ci<-vector(length=10)
    for (i in 1:5) {
      ui[(2*i-1),i]<-1
      ui[(2*i),i]<-(-1)
    }
    ci[1:2]<-c(range(topt)[1],-range(topt)[2])
    aux<-sort(c(0.5*argv$gamma.standard,2*argv$gamma.standard))
    ci[3:4]<-c(aux[1],-aux[2])
    ci[5:6]<-c(0,-20)
    if (as.numeric(quantile(zopt,probs=0.1))==
        as.numeric(quantile(zopt,probs=0.25))) {
      aux<-as.numeric(quantile(zopt,probs=0.25))+
           mean(zopt)/10
    } else {
      aux<-as.numeric(quantile(zopt,probs=0.25))
    }
    ci[7:8]<-c(-10,-aux)
    aux<-as.numeric(quantile(zopt,probs=0.75))-as.numeric(quantile(zopt,probs=0.25))
    ci[9:10]<-c(0,-aux)
    par<-c(mean(range(topt)),
           argv$gamma.standard,
           5,
           as.numeric(quantile(zopt,probs=0.1)),
           -ci[10]/2)
    opt<-constrOptim(theta=par,
                     f=tvertprofFrei2opt,
                     ui=ui,
                     ci=ci,
                     grad=NULL)
    tpar<-c(2,opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5])
  }
  # aggregate obs to speed-up the elaboration
  raux<-raster(xmn=(xr[ixynp[1]]-fact*res[1]),xmx=(xr[ixynp[1]]+fact*res[1]),
               ymn=(yr[ixynp[1]]-fact*res[2]),ymx=(yr[ixynp[1]]+fact*res[2]),
               res=c((2*(fact*res[1])/10),(2*(fact*res[2])/10)))
  raux[]<-1:ncell(raux)
  nobs_aux<-getValues(rasterize(cbind(VecX_bg[ix],VecY_bg[ix]),
                      raux,yo_bg[ix],fun=function(x,...)length(x)))
  ii_aux<-extract(raux,cbind(VecX_bg[ix],VecY_bg[ix]))
  ixnr_aux<-which(!is.na(nobs_aux) & nobs_aux>0)
  nout_aux<-length(ixnr_aux)
  x_aux<-vector()
  y_aux<-vector()
  for (i in 1:nout_aux) {
    x_aux[i]<-mean(VecX_bg[ix][which(ii_aux==ixnr_aux[i])])
    y_aux[i]<-mean(VecY_bg[ix][which(ii_aux==ixnr_aux[i])])
  }
  #
  # (S+R) definition, inversion and sum-up rows. Just in one line of code.
  disth<-sqrt(outer(x_aux,x_aux,FUN="-")**2.+
              outer(y_aux,y_aux,FUN="-")**2.)
  # maxboxl/6 is the max limit for dh
  dh_oi<-min((maxboxl/6.),mean(findRow(x=disth,n=(length(x_aux)-closeNth+1))))
  # more than 100 bg_obs (10*10) then reduce dh_oi
  if (length(ix)>100) dh_oi<-dh_oi/2.
  vec1<-rowSums(chol2inv(chol( 
                ( eps2*diag(length(x_aux)) + exp(-0.5*(disth*disth)/(refdist*refdist)) ) )))
  rm(disth)
  if (!argv$cv_mode) {
    # grid
    ixg<-which(((min(x_aux)-xgrid)<=(7*refdist)) & 
               ((xgrid-max(x_aux))<=(7*refdist)) &
               ((min(y_aux)-ygrid)<=(7*refdist)) & 
               ((ygrid-max(y_aux))<=(7*refdist)) )
    out<-OI_T_xb_upd( xg=xgrid[ixg],
                      yg=ygrid[ixg],
                      zg=dem[ixg],
                      bg=xb[ixg],
                      wg=xw[ixg],
                      dg=xdh_oi[ixg],
                      dg_new=dh_oi,
                      xo=x_aux,
                      yo=y_aux,
                      vec1=vec1,
                      tpar=tpar,
                      na=na,
                      Dh=refdist)
    xb[ixg]<-out$bg_up
    xw[ixg]<-out$wg_up
    xdh_oi[ixg]<-out$dg_up
    assign("xb",xb,envir=.GlobalEnv)
    assign("xw",xw,envir=.GlobalEnv)
    assign("xdh_oi",xdh_oi,envir=.GlobalEnv)
  # CVmode
  } else {
    ixg<-which(((min(x_aux)-VecX_cv)<=(7*refdist)) & 
               ((VecX_cv-max(x_aux))<=(7*refdist)) &
               ((min(y_aux)-VecY_cv)<=(7*refdist)) & 
               ((VecY_cv-max(y_aux))<=(7*refdist)) )
    out<-OI_T_xb_upd( xg=VecX_cv[ixg],
                      yg=VecY_cv[ixg],
                      zg=VecZ_cv[ixg],
                      bg=xb[ixg],
                      wg=xw[ixg],
                      dg=xdh_oi[ixg],
                      dg_new=dh_oi,
                      xo=x_aux,
                      yo=y_aux,
                      vec1=vec1,
                      tpar=tpar,
                      na=na,
                      Dh=refdist)
    xb[ixg]<-out$bg_up
    xw[ixg]<-out$wg_up
    xdh_oi[ixg]<-out$dg_up
    assign("xb",xb,envir=.GlobalEnv)
    assign("xw",xw,envir=.GlobalEnv)
    assign("xdh_oi",xdh_oi,envir=.GlobalEnv)
  }
#  print(paste(length(ix),length(ixg),round(refdist,0),round(mean(fact*res),0)))
#  print(paste(length(ix),round(dh_oi,0),round(mean(fact*res),0)))
  # station points
  ixg<-which(((min(x_aux)-VecX)<=(7*refdist)) & 
             ((VecX-max(x_aux))<=(7*refdist)) &
             ((min(y_aux)-VecY)<=(7*refdist)) & 
             ((VecY-max(y_aux))<=(7*refdist)) )
  out<-OI_T_xb_upd( xg=VecX[ixg],
                    yg=VecY[ixg],
                    zg=VecZ[ixg],
                    bg=yb[ixg],
                    wg=yw[ixg],
                    dg=ydh_oi[ixg],
                    dg_new=dh_oi,
                    xo=x_aux,
                    yo=y_aux,
                    vec1=vec1,
                    tpar=tpar,
                    na=na,
                    Dh=refdist)
  yb[ixg]<-out$bg_up
  yw[ixg]<-out$wg_up
  ydh_oi[ixg]<-out$dg_up
  assign("yb",yb,envir=.GlobalEnv)
  assign("yw",yw,envir=.GlobalEnv)
  assign("ydh_oi",ydh_oi,envir=.GlobalEnv)
  #
  if (argv$debug) {
    png(file=file.path(argv$debug.dir,
                       paste0("deb_sumOfweights_",
                              formatC(ixynp[1],width=4,flag="0"),".png")),
       width=800,height=800)
    image(rlaf,xlim=c(xmn,xmx),ylim=c(ymn,ymx),breaks=c(-1000,1000),col="beige",
          main=paste(length(ix),round(refdist,0),round(mean(fact*res),0)) )
    points(VecX_bg,VecY_bg,cex=0.5)
    rect(xr[ixynp[1]]-fact*res[1],yr[ixynp[1]]-fact*res[2],
         xr[ixynp[1]]+fact*res[1],yr[ixynp[1]]+fact*res[2])
    points(xr[irx],yr[irx],pch=19,col="gray",cex=0.75)
    points(VecX_bg[ix],VecY_bg[ix],pch=19,col="blue")
    points(x_aux,y_aux,pch=19,col="cyan",cex=1.25)
    points(ixynp[2],ixynp[3],pch=19,col="red")
    dev.off()
    # grid
    ixg<-which(((min(x_aux)-xgrid)<=(2*refdist)) & 
               ((xgrid-max(x_aux))<=(2*refdist)) &
               ((min(y_aux)-ygrid)<=(2*refdist)) & 
               ((ygrid-max(y_aux))<=(2*refdist)) )
    if (tpar[1]==2) {
      xb_aux<-tvertprof_Frei(z=dem[ixg],
                             t0=tpar[2],
                             gamma=tpar[3],
                             a=tpar[4],
                             h0=tpar[5],
                             h1i=tpar[6])
    } else if (tpar[1]==1) {
      xb_aux<-tvertprof_basic(z=dem[ixg],
                              t0=tpar[2],
                              gamma=tpar[3])
    }
    png(file=file.path(argv$debug.dir,
                       paste0("deb_vertprof_",
                              formatC(ixynp[1],width=4,flag="0"),".png")),
       width=800,height=800)
    mnt<-min(c(xb[ixg],xb_aux,yo_bg[ix]))
    mxt<-max(c(xb[ixg],xb_aux,yo_bg[ix]))
    plot(xb_aux,dem[ixg],
         pch=19,col="white",main=round(tpar,4),
         xlim=c(mnt,mxt))
    points(xb[ixg],dem[ixg],pch=19,col="gray")
    points(yo_bg[ix],VecZ_bg[ix],pch=19,col="blue")
    points(xb_aux,dem[ixg],pch=19,col="black")
    dev.off()
    file<-file.path(argv$debug.dir,
                    paste0("deb_backg_",
                           formatC(ixynp[1],width=4,flag="0"),".png"))
    rdeb<-rmaster
    xb[which(xb==na)]<-NA
    rdeb[mask]<-xb
    writeRaster(rdeb, file, format="CDF",overwrite=T)
    file<-file.path(argv$debug.dir,
                    paste0("deb_weights_",
                           formatC(ixynp[1],width=4,flag="0"),".png"))
    rdeb<-rmaster
    xw[which(xw==na)]<-NA
    rdeb[mask]<-xw
    writeRaster(rdeb, file, format="CDF",overwrite=T)
    file<-file.path(argv$debug.dir,
                    paste0("deb_dhoi_",
                           formatC(ixynp[1],width=4,flag="0"),".png"))
    rdeb<-rmaster
    xdh_oi[which(xdh_oi==na)]<-NA
    rdeb[mask]<-xdh_oi
    writeRaster(rdeb, file, format="CDF",overwrite=T)
  }
}

#+ find the n-th largest element from each matrix row 
findRow <- function (x,n) {   
# order by row number then by value
  y<-t(x)
#  array(y[order(col(y), y)], dim(y))[nrow(y) - 1, ]
  array(y[order(col(y), y)], dim(y))[nrow(y) - (n-1), ]
}

#+ mean radial distance between an observation and its k-th closest obs
dobs_fun<-function(ixynp,k) {
# NOTE: k=1 will return 0 everywhere (1st closest obs is the obs itself)
  # jj, index for the stations in the box
  jj<-which(iobs==ixynp[1])
  njj<-length(jj)
  if (njj<k) return(NA)
  # distance matrices
  disth<-(outer(VecX[jj],VecX[jj],FUN="-")**2.+
          outer(VecY[jj],VecY[jj],FUN="-")**2.)**0.5
#  distz<-abs(outer(data$z[jj],data$z[jj],FUN="-"))
  dobsRow<-findRow(x=disth,n=(njj-k+1))
  mean(dobsRow)
}

# superobbing for the ith box (divide it in a n1 x n2 boxes) 
superobs<-function(ixy,         #i-index,x,y
                   res,         #resx,resy
                   n=c(1,1),
                   nmin=2) {  #nx,ny
  ixI<-which(VecI==ixy[1])
  if (n[1]==1 & n[2]==1) {
    nobs<-length(ixI)
    if (nobs>nmin) {
      nout<-1
      x_so<-mean(VecX[ixI])
      y_so<-mean(VecY[ixI])
      z_so<-mean(VecZ[ixI])
      laf_so<-mean(VecLaf[ixI])
      yo_so<-mean(yo[ixI])
      rm_me[ixI]<-T
      assign("rm_me",rm_me,envir=.GlobalEnv)
    } else {
      nout<-0
      x_so<-NULL
      y_so<-NULL
      z_so<-NULL
      laf_so<-NULL
      yo_so<-NULL
    }
  } else {
    r<-raster(xmn=ixy[2]-res[1]/2,xmx=ixy[2]+res[1]/2,
              ymn=ixy[3]-res[2]/2,ymx=ixy[3]+res[2]/2,
              res=c(res[1]/n[1],res[2]/n[2]))
    r[]<-1:ncell(r)
    nobs<-getValues(
           rasterize(cbind(VecX[ixI],VecY[ixI]),
                     r,
                     yo[ixI],
                     fun=function(x,...)length(x)) )
    ii<-extract(r,cbind(VecX[ixI],VecY[ixI]))
    ixnr<-which(!is.na(nobs) & nobs>nmin)
    nout<-length(ixnr)
    if (nout>0) {
      x_so<-vector()
      y_so<-vector()
      z_so<-vector()
      laf_so<-vector()
      yo_so<-vector()
      for (i in 1:nout) {
        x_so[i]<-mean(VecX[ixI][which(ii==ixnr[i])])
        y_so[i]<-mean(VecY[ixI][which(ii==ixnr[i])])
        z_so[i]<-mean(VecZ[ixI][which(ii==ixnr[i])])
        laf_so[i]<-mean(VecLaf[ixI][which(ii==ixnr[i])])
        yo_so[i]<-mean(yo[ixI][which(ii==ixnr[i])])
        rm_me[ixI][which(ii==ixnr[i])]<-T
      }
      assign("rm_me",rm_me,envir=.GlobalEnv)
    } else {
      x_so<-NULL
      y_so<-NULL
      z_so<-NULL
      laf_so<-NULL
      yo_so<-NULL
    }
  }
  return(list(n=nout,x=x_so,y=y_so,z=z_so,laf=laf_so,yo=yo_so))
}

#+ OI for temperature (variable Dh) 
oiIT<-function(par,eps2=0.25,dz=250,lafmn=0.5,cv=FALSE,nmaxo=20) {
# 1=x,2=y,3=z,4=laf,5=dh,6=xb,7=index
  if ((par[7]%%100000)==0) print(paste(par[7],length_tot))
  if (any(is.na(par))) {
    if (cv) {
      return(rep(NA,4))
    } else {
      return(rep(NA,2))
    }
  }
  ixo<-integer(0)
  n<-0
  while (length(ixo)<nmaxo) {
    n<-n+.25
    ixo<-which(abs(VecX-par[1])<(n*par[5]) & abs(VecY-par[2])<(n*par[5]))
  }
  disth.i<-(par[1]-VecX[ixo])*(par[1]-VecX[ixo]) + 
           (par[2]-VecY[ixo])*(par[2]-VecY[ixo])
  # indicies to the nmaxo closet obs
  ixo<-ixo[order(disth.i)[1:nmaxo]]
  x<-VecX[ixo]
  y<-VecY[ixo]
  z<-VecZ[ixo]
  laf<-VecLaf[ixo]
  inn<-innov[ixo]
  mat<-chol2inv(chol( 
                ( eps2*diag(length(ixo)) + 
   (exp(-0.5*(outer(x,x,FUN="-")**2.+outer(y,y,FUN="-")**2.)/(par[5]*par[5])
        -0.5*(outer(z,z,FUN="-")**2.)/(dz*dz)) * 
    (1-(1-lafmn)*abs(outer(laf,laf,FUN="-"))) )  ) ))
  vec1<-crossprod(mat,inn)
  vec2<-rowSums(mat)
  g_row<-exp(-0.5*((par[1]-x)*(par[1]-x)+(par[2]-y)*(par[2]-y))/(par[5]*par[5])
          -0.5*((par[3]-z)*(par[3]-z)/(dz*dz)) ) *
         (1-(1-lafmn)*abs(par[4]-laf))
  xa<-par[6] + crossprod( g_row, vec1 )
  xidi<-crossprod( g_row, vec2 )
#  print(paste(par[6],xa))
  if (!cv) return(c(xa,xidi))
  coeff<-1./(1.-tcrossprod( g_row, mat )[which(ixo==par[7])])
  xav<-yo[par[7]] + coeff * (xa-yo[par[7]]) 
  xidiv<-1 + coeff * (xidi-1) 
  return(c(xa,xidi,xav,xidiv))
}


#
#==============================================================================
# MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN -
#==============================================================================
t0<-Sys.time()
# create parser object
p <- arg_parser("bliss")
#------------------------------------------------------------------------------
# miscellaneous 
p <- add_argument(p, "--verbose",
                  help="verbose mode",
                  flag=T,
                  short="-v")
p <- add_argument(p, "--debug",
                  help="debug mode",
                  flag=T)
p <- add_argument(p, "--rrinf",
                  help="precipitation yes/no threshold",
                  type="numeric",
                  default=0.1)
#------------------------------------------------------------------------------
# cross-validation mode
p <- add_argument(p, "--cv_mode",
                  help="standard cross-validation mode",
                  flag=T)
p <- add_argument(p, "--loocv_mode",
                  help="leave-one-out cross-validation mode",
                  type="logical",
                  default=F)
p <- add_argument(p, "--idiv_instead_of_elev",
                  help="leave-one-out cross-validation mode",
                  type="logical",
                  default=F)
p <- add_argument(p, "--twostep_superobbing",
                  help="superobbing (used only if \"OI_twosteptemperature\")",
                  flag=T)
#------------------------------------------------------------------------------
# statistical interpolation mode
p <- add_argument(p, "--mode",
                  help="statistical interpolation scheme (\"OI_multiscale\",\"OI_firstguess\",\"OI_twosteptemperature\")",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# time-related variables
p <- add_argument(p, "--date_out",
                  help="date to write in output files (nc, %Y%m%d%H%M)",
                  type="character",
                  default=NULL)
p <- add_argument(p, "--date_out_fmt",
                  help="date out format",
                  type="character",
                  default="%Y%m%d%H%M")
p <- add_argument(p, "--time_bnds_string",
                  help="time bounds with respect to date_out (e.g., \"-1 day\" \"-1 min\")",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# OI_multiscale / OI_firstguess parameters
p <- add_argument(p, "--eps2",
                  help="ratio of observation to background error covariance",
                  type="numeric",
                  default=1)
p <- add_argument(p, "--Dh",
                  help="horizontal de-corellation length scale (km)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--ovarc",
                  help="observation error variance correction factor",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--ovarc.prId",
                  help="observation error variance correction factor (prId)",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--prId.exclude",
                  help="observation provider identifiers to exclude",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
p <- add_argument(p, "--prId.cv",
                  help="observation provider identifiers to reserve for cross-validation",
                  type="numeric",
                  default=NULL,
                  nargs=Inf)
# additional OI parameters two-step temperature
p <- add_argument(p, "--nmaxo",
                  help="number of closest observations to be used in the OI to adjust background values at a grid point",
                  type="numeric",
                  default=20)
p <- add_argument(p, "--dz",
                  help="background error covariance matrix, vertical decorrelation length scale (m)",
                  type="numeric",
                  default=600)
p <- add_argument(p, "--lafmin",
                  help="background error covariance matrix, land area fraction parameter",
                  type="numeric",
                  default=0.5)
# Background parameters
p <- add_argument(p, "--grid.bg",
                  help="nrow ncol (i.e. number_of_rows number_of_columns) used to define grid of sub-regions.",
                  type="integer",
                  nargs=2,
                  default=c(5,5))
p <- add_argument(p, "--vmin",
                  help="minimum allowed value [units of the variable specified]",
                  type="numeric",
                  default=-50)
p <- add_argument(p, "--vmax",
                  help="maximum allowed value [units of the variable specified]",
                  type="numeric",
                  default=40)
p <- add_argument(p, "--gamma.standard",
                  help="standard value for the moist adiabatic temperature lapse rate dT/dz (degC/m)",
                  type="numeric",
                  default=-0.0065)
p <- add_argument(p, "--n.bg",
                  help="minimum number of observations allowed in a sub-region",
                  type="numeric",
                  default=50)
p <- add_argument(p, "--dz.bg",
                  help="elevation difference. If the elevation range (=elev 95th-perc minus elev 5th-perc) is less than this value, then fit a linear profile of temperature. Otherwise, fit a non-linear profile.",
                  type="numeric",
                  default=30)
p <- add_argument(p, "--nclose.bg",
                  help="n-th closest observation to consider in the calculation of the OI de-correlation length",
                  type="numeric",
                  default=4)
p <- add_argument(p, "--maxboxl",
                  help="maximum lenght (m) of the box used to define a sub-region",
                  type="numeric",
                  default=250000)
p <- add_argument(p, "--obs.outbuffer",
                  help="distance (m) defining the \"buffer\" region outside the masked region where to consider the observation, so to reduce border effects",
                  type="numeric",
                  default=50000)
#------------------------------------------------------------------------------
# paths
p <- add_argument(p, "--path2src",
                  help="path to the shared objects (.so files)",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# input files
p <- add_argument(p, "--config.file",
                  help="configuration file",
                  type="character",
                  default=NULL,
                  short="cf")
p <- add_argument(p, "--iff_obs",
                  help="full file name for the observations (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg",
                  help="full file name for the first-guess field (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf",
                  help="full file name for the rescaling factor (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf",
                  help="full file name for the land area fraction (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask",
                  help="full file name for the mask (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_black",
                  help="blacklist file",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# output files
p <- add_argument(p, "--off_stn",
                  help="full file name for output at station locations (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_ver",
                  help="full file name for output at station locations analysis vs obs - verif format (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_ver_fg",
                  help="full file name for output at station locations fist-guess vs obs - verif format (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_ver_cv",
                  help="full file name for output at station locations cv-analysis vs obs - verif format (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--verif_var",
                  help="name of the verif variable",
                  type="character",
                  default="Precipitation")
p <- add_argument(p, "--verif_units",
                  help="units of the verif variable",
                  type="character",
                  default="mm")
p <- add_argument(p, "--off_grd",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_cvstn",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_cvver",
                  help="full file name for output at gridpoints (nc)",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# customization of the observation file
p <- add_argument(p, "--iff_obs.sep",
                  help="separator character",
                  type="character",
                  default=";")
p <- add_argument(p, "--iff_obs.x",
                  help="easting coordinate name",
                  type="character",
                  default="x")
p <- add_argument(p, "--iff_obs.y",
                  help="northing coordinate name",
                  type="character",
                  default="y")
p <- add_argument(p, "--iff_obs.z",
                  help="elevation name",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.value",
                  help="variable name",
                  type="character",
                  default="value")
p <- add_argument(p, "--iff_obs.proj4",
                  help="proj4 string for the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.sourceId",
                  help="station identifier",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.prId",
                  help="provider identifier",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_obs.dqc",
                  help="data quality control",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# Master grid definition
p <- add_argument(p, "--grid_master.x1",
                  help="easting coordinate of the first gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.xn",
                  help="easting coordinate of the last gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.y1",
                  help="northing coordinate of the first gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.yn",
                  help="northing coordinate of the last gridpoint (e.g., value returned by ncdump -c)",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.resx",
                  help="grid spacing along the easting coordinate",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.resy",
                  help="grid spacing along the northing coordinate",
                  type="character",
                  default="none")
p <- add_argument(p, "--grid_master.proj4",
                  help="proj4 string for the master grid",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# Mask file netcdf parameters
p <- add_argument(p, "--iff_mask.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_mask.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_mask.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_mask.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_mask.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# Rescaling factor file netcdf parameters
p <- add_argument(p, "--iff_rf.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.varname_lat",
                  help="name of the latitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.varname_lon",
                  help="name of the longitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_rf.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_rf.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_rf.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_rf.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# Land area fraction file netcdf parameters
p <- add_argument(p, "--iff_laf.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.varname_lat",
                  help="name of the latitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.varname_lon",
                  help="name of the longitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_laf.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_laf.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_laf.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_laf.adjfact",
                  help="correction factor (laf should be 0-1)",
                  type="numeric",
                  default=1)
#------------------------------------------------------------------------------
# Land area fraction file netcdf parameters
p <- add_argument(p, "--iff_dem.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.varname_lat",
                  help="name of the latitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.varname_lon",
                  help="name of the longitude variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_dem.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_dem.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_dem.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_dem.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_dem.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_dem.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# first-guess file netcdf parameters
p <- add_argument(p, "--iff_fg.varname",
                  help="name of the variable to read from the file",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.topdown",
                  help="turn the field upside-down",
                  type="logical",
                  default=F)
p <- add_argument(p, "--iff_fg.ndim",
                  help="number of dimensions for the variable",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.tpos",
                  help="position of the time variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.epos",
                  help="position of the ensemble_member variable among the dimensions",
                  type="numeric",
                  default=NULL)
p <- add_argument(p, "--iff_fg.names",
                  help="dimension names for the variable",
                  type="character",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--iff_fg.proj4",
                  help="proj4 string identyfing the coordinate reference system",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.t",
                  help="timestamp to read from file (defualt, read the first)",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.tfmt",
                  help="timestamp time format",
                  type="character",
                  default="none")
p <- add_argument(p, "--iff_fg.e",
                  help="label of the ensemble member to read from file (default is null)",
                  type="numeric",
                  default=NULL)
#------------------------------------------------------------------------------
# output file netcdf parameters
p <- add_argument(p, "--off_grd.grid",
                  help="grid type",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.variables",
                  help="type of variable (analysis,background,idi)",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_grd.varname",
                  help="variable name",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_grd.varlongname",
                  help="variable long name",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.standardname",
                  help="variable standard name",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_grd.varversion",
                  help="variable version",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.varunit",
                  help="variable unit",
                  type="character",
                  nargs=Inf,
                  default=NA)
p <- add_argument(p, "--off_grd.timesunit",
                  help="time unit",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.reference",
                  help="references",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.write_lonlat",
                  help="add latitude and longitude variables",
                  type="logical",
                  default=F)
p <- add_argument(p, "--off_grd.diground",
                  help="rounding digits",
                  type="numeric",
                  default=3)
p <- add_argument(p, "--off_grd.summary",
                  help="summary",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.sourcestring",
                  help="source string",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.title",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.comment",
                  help="title",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_grd.cell_methods",
                  help="title",
                  type="character",
                  default="none")
#------------------------------------------------------------------------------
# gaussian anamorphosis
p <- add_argument(p, "--transf",
                  help="transformation used in the gaussian anamorphosis (\"none\",\"Box-Cox\")",
                  type="character",
                  default="none")
p <- add_argument(p, "--transf.boxcox_lambda",
                  help="Box-Cox power parameter",
                  type="numeric",
                  default=NA)

#------------------------------------------------------------------------------
#
argv <- parse_args(p)
#
#------------------------------------------------------------------------------
# read configuration file
if (!is.na(argv$config.file)) {
  if (file.exists(argv$config.file)) {
    source(argv$config.file)
    argv_tmp<-append(argv,conf)
    names_argv_tmp<-names(argv_tmp)
    argv_def<-list()
    names_argv_def<-integer(0)
    k<-0
    for (i in 1:length(argv_tmp)) {
      if (names_argv_tmp[i] %in% names_argv_def) next
      k<-k+1
      j<-which(names_argv_tmp==names_argv_tmp[i])
      argv_def[[k]]<-argv_tmp[[j[length(j)]]]
      names_argv_def<-c(names_argv_def,names_argv_tmp[i])
    }
    names(argv_def)<-names_argv_def
    rm(argv_tmp,names_argv_tmp,names_argv_def)
    rm(argv)
    argv<-argv_def
    rm(argv_def)
  } else {
    print("WARNING: config file not found")
    print(argv$config.file)
  }
}
#
#------------------------------------------------------------------------------
# check input arguments
#
if (!(file.exists(argv$iff_obs))) boom(paste0("file not found ",argv$iff_obs))
if (argv$mode=="OI_multiscale") {
  if (argv$verbose) {
    if (!(file.exists(argv$iff_rf))) 
      print("warning: file not found",argv$iff_rf)
  }
} else if (argv$mode=="OI_firstguess") {
  if (!(file.exists(argv$iff_fg))) 
    boom(paste0("file not found ",argv$iff_fg))
} else if (argv$mode=="OI_twosteptemperature") {
  if (argv$verbose) {
    if (!(file.exists(argv$iff_dem))) 
      print("warning: file not found",argv$iff_dem)
  }
} else {
  boom("error statistical interpolation scheme undefined")
}
# define/check paths and load external functions
if ( !(file.exists(argv$path2src)) ) 
  ext<-boom("path not found")
#
# load external C functions
dyn.load(file.path(argv$path2src,"oi_rr_first.so"))
dyn.load(file.path(argv$path2src,"oi_rr_fast.so"))
dyn.load(file.path(argv$path2src,"oi_rr_var.so"))
dyn.load(file.path(argv$path2src,"oi_t_xb_upd.so"))
#
#------------------------------------------------------------------------------
# Create master grid
xmn<-as.numeric(argv$grid_master.x1)-as.numeric(argv$grid_master.resx)/2
xmx<-as.numeric(argv$grid_master.xn)+as.numeric(argv$grid_master.resx)/2
ymn<-as.numeric(argv$grid_master.y1)-as.numeric(argv$grid_master.resy)/2
ymx<-as.numeric(argv$grid_master.yn)+as.numeric(argv$grid_master.resy)/2
rmaster<-raster(extent(xmn,xmx,ymn,ymx),
                res=c(as.numeric(argv$grid_master.resx),
                      as.numeric(argv$grid_master.resy)),
                crs=argv$grid_master.proj4)
rmaster[]<-1
# use mask if provided
if (file.exists(argv$iff_mask)) {
  argv$iff_mask.epos<-set_NAs_to_NULL(argv$iff_mask.epos)
  argv$iff_mask.tpos<-set_NAs_to_NULL(argv$iff_mask.tpos)
  argv$iff_mask.e<-set_NAs_to_NULL(argv$iff_mask.e)
  if (argv$iff_mask.t=="none") argv$iff_mask.t<-nc4.getTime(argv$iff_mask)[1]
  raux<-try(read_dotnc(nc.file=argv$iff_mask,
                       nc.varname=argv$iff_mask.varname,
                       topdown=argv$iff_mask.topdown,
                       out.dim=list(ndim=argv$iff_mask.ndim,
                                    tpos=argv$iff_mask.tpos,
                                    epos=argv$iff_mask.epos,
                                    names=argv$iff_mask.names),
                       proj4=argv$iff_mask.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_mask.t,
                                      format=argv$iff_mask.tfmt,
                                      e=argv$iff_mask.e)))
  if (is.null(raux)) 
    boom("error reading the mask file")
  rmask<-raux$stack; rm(raux)
  if (!rasters_match(rmask,rmaster)) rmask<-projectRaster(rmask,rmaster)
  rmaster<-mask(rmaster,rmask); rm(rmask)
}
#
nx<-ncol(rmaster)
ny<-nrow(rmaster)
xmn<-xmin(rmaster)
xmx<-xmax(rmaster)
ymn<-ymin(rmaster)
ymx<-ymax(rmaster)
#
# extract all the cell values: zvalues[1] contains the rmaster[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
zvalues<-getValues(rmaster)
storage.mode(zvalues)<-"numeric"
xy<-xyFromCell(rmaster,1:ncell(rmaster))
x<-sort(unique(xy[,1]))
y<-sort(unique(xy[,2]),decreasing=T)
xgrid<-xy[,1]
ygrid<-xy[,2]
mask<-which(!is.na(zvalues))
ngrid<-length(mask)
# clean memory
rm(zvalues,xy)
# debug info
if (argv$verbose) {
  print("+---------------------------------------------------------------+")
  print("+ grid parameters")
  print(paste("nx ny dx dy",
    as.integer(nx),
    as.integer(ny),
    round(xres(rmaster),2),
    round(yres(rmaster),2)))
  print(paste("xmn xmx ymn ymx",
    round(xmn,2),
    round(xmx,2),
    round(ymn,2),
    round(ymx,2)))
  print(paste("# grid points (master/unmasked)=",as.integer(ngrid)))
}
#
#------------------------------------------------------------------------------
# read rescaling factor
if (file.exists(argv$iff_rf)) {
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print(paste("+ rescaling factor",argv$iff_rf))
  }
  argv$iff_rf.epos<-set_NAs_to_NULL(argv$iff_rf.epos)
  argv$iff_rf.tpos<-set_NAs_to_NULL(argv$iff_rf.tpos)
  argv$iff_rf.e<-set_NAs_to_NULL(argv$iff_rf.e)
  if (any(argv$iff_rf.names=="time")) {
    if (argv$iff_rf.t=="none") argv$iff_rf.t<-nc4.getTime(argv$iff_rf)[1]
  } else {
    argv$iff_rf.t<-NULL
  }
  raux<-try(read_dotnc(nc.file=argv$iff_rf,
                       nc.varname=argv$iff_rf.varname,
                       topdown=argv$iff_rf.topdown,
                       out.dim=list(ndim=argv$iff_rf.ndim,
                                    tpos=argv$iff_rf.tpos,
                                    epos=argv$iff_rf.epos,
                                    names=argv$iff_rf.names),
                       proj4=argv$iff_rf.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_rf.t,
                                      format=argv$iff_rf.tfmt,
                                      e=argv$iff_rf.e)))
  if (is.null(raux)) 
    boom("error reading the rescaling file")
  rrf<-raux$stack; rm(raux)
  if (!rasters_match(rrf,rmaster)) {
    if (argv$iff_rf.varname_lat=="none") {
      rrf<-projectRaster(rrf,rmaster)
      rf<-getValues(rrf)
    } else {
      raux<-try(read_dotnc(nc.file=argv$iff_rf,
                       nc.varname=argv$iff_rf.varname_lat,
                       topdown=argv$iff_rf.topdown,
                       out.dim=list(ndim=argv$iff_rf.ndim,
                                    tpos=argv$iff_rf.tpos,
                                    epos=argv$iff_rf.epos,
                                    names=argv$iff_rf.names),
                       proj4=argv$iff_rf.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_rf.t,
                                      format=argv$iff_rf.tfmt,
                                      e=argv$iff_rf.e)))
      if (is.null(raux)) 
        boom("error reading the rescaling file (lat)")
      lat<-raux$data; rm(raux)
      raux<-try(read_dotnc(nc.file=argv$iff_rf,
                           nc.varname=argv$iff_rf.varname_lon,
                           topdown=argv$iff_rf.topdown,
                           out.dim=list(ndim=argv$iff_rf.ndim,
                                        tpos=argv$iff_rf.tpos,
                                        epos=argv$iff_rf.epos,
                                        names=argv$iff_rf.names),
                           proj4=argv$iff_rf.proj4,
                           nc.proj4=list(var=NULL,
                                         att=NULL),
                           selection=list(t=argv$iff_rf.t,
                                          format=argv$iff_rf.tfmt,
                                          e=argv$iff_rf.e)))
      if (is.null(raux)) 
        boom("error reading the rescaling file (lon)")
      lon<-raux$data; rm(raux)
      coord.new<-spTransform(SpatialPoints(cbind(lon,lat),
                                           proj4string=argv$iff_rf.proj4),
                             CRS(argv$grid_master.proj4))
      #
      rf<-getValues(rrf)
      rrfagg<-rasterize(coord.new,
                        aggregate(rmaster,fact=4),
                        rf)
      rrf<-mask( crop( disaggregate(rrfagg,fact=4,method="bilinear"),
                       rmaster),
                 rmaster)
      rf<-getValues(rrf)
    }
  }
}
#
#------------------------------------------------------------------------------
# read land area fraction
if (file.exists(argv$iff_laf)) {
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print(paste("+ land area fraction",argv$iff_laf))
  }
  argv$iff_laf.epos<-set_NAs_to_NULL(argv$iff_laf.epos)
  argv$iff_laf.tpos<-set_NAs_to_NULL(argv$iff_laf.tpos)
  argv$iff_laf.e<-set_NAs_to_NULL(argv$iff_laf.e)
  if (any(argv$iff_laf.names=="time")) {
    if (argv$iff_laf.t=="none") argv$iff_laf.t<-nc4.getTime(argv$iff_laf)[1]
  } else {
    argv$iff_laf.t<-NULL
  }
  raux<-try(read_dotnc(nc.file=argv$iff_laf,
                       nc.varname=argv$iff_laf.varname,
                       topdown=argv$iff_laf.topdown,
                       out.dim=list(ndim=argv$iff_laf.ndim,
                                    tpos=argv$iff_laf.tpos,
                                    epos=argv$iff_laf.epos,
                                    names=argv$iff_laf.names),
                       proj4=argv$iff_laf.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_laf.t,
                                      format=argv$iff_laf.tfmt,
                                      e=argv$iff_laf.e)))
  if (is.null(raux)) 
    boom("error reading the land area fraction file")
  rlaf<-raux$stack; rm(raux)
  laf<-getValues(rlaf)*argv$iff_laf.adjfact
  if (!rasters_match(rlaf,rmaster)) {
    if (argv$iff_laf.varname_lat=="none") {
      rlaf<-projectRaster(rlaf,rmaster)
      laf<-getValues(rlaf)*argv$iff_laf.adjfact
      rlaf[]<-laf
    } else {
      raux<-try(read_dotnc(nc.file=argv$iff_laf,
                       nc.varname=argv$iff_laf.varname_lat,
                       topdown=argv$iff_laf.topdown,
                       out.dim=list(ndim=argv$iff_laf.ndim,
                                    tpos=argv$iff_laf.tpos,
                                    epos=argv$iff_laf.epos,
                                    names=argv$iff_laf.names),
                       proj4=argv$iff_laf.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_laf.t,
                                      format=argv$iff_laf.tfmt,
                                      e=argv$iff_laf.e)))
      if (is.null(raux)) 
        boom("error reading the land area fraction (lat)")
      lat<-raux$data; rm(raux)
      raux<-try(read_dotnc(nc.file=argv$iff_laf,
                           nc.varname=argv$iff_laf.varname_lon,
                           topdown=argv$iff_laf.topdown,
                           out.dim=list(ndim=argv$iff_laf.ndim,
                                        tpos=argv$iff_laf.tpos,
                                        epos=argv$iff_laf.epos,
                                        names=argv$iff_laf.names),
                           proj4=argv$iff_laf.proj4,
                           nc.proj4=list(var=NULL,
                                         att=NULL),
                           selection=list(t=argv$iff_laf.t,
                                          format=argv$iff_laf.tfmt,
                                          e=argv$iff_laf.e)))
      if (is.null(raux)) 
        boom("error reading the land area fraction (lon)")
      lon<-raux$data; rm(raux)
      coord.new<-spTransform(SpatialPoints(cbind(lon,lat),
                                           proj4string=argv$iff_laf.proj4),
                             CRS(argv$grid_master.proj4))
      #
      laf<-getValues(rlaf)
      rlafagg<-rasterize(coord.new,
                        aggregate(rmaster,fact=4),
                        laf)
      rlaf<-mask( crop( disaggregate(rlafagg,fact=4,method="bilinear"),
                       rmaster),
                 rmaster)
      laf<-getValues(rlaf)*argv$iff_laf.adjfact
      rlaf[]<-laf
    }
  }
}
#
#------------------------------------------------------------------------------
# read digital elevation model 
if (file.exists(argv$iff_dem)) {
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print(paste("+ digital elevation model",argv$iff_dem))
  }
  argv$iff_dem.epos<-set_NAs_to_NULL(argv$iff_dem.epos)
  argv$iff_dem.tpos<-set_NAs_to_NULL(argv$iff_dem.tpos)
  argv$iff_dem.e<-set_NAs_to_NULL(argv$iff_dem.e)
  if (any(argv$iff_dem.names=="time")) {
    if (argv$iff_dem.t=="none") argv$iff_dem.t<-nc4.getTime(argv$iff_dem)[1]
  } else {
    argv$iff_dem.t<-NULL
  }
  raux<-try(read_dotnc(nc.file=argv$iff_dem,
                       nc.varname=argv$iff_dem.varname,
                       topdown=argv$iff_dem.topdown,
                       out.dim=list(ndim=argv$iff_dem.ndim,
                                    tpos=argv$iff_dem.tpos,
                                    epos=argv$iff_dem.epos,
                                    names=argv$iff_dem.names),
                       proj4=argv$iff_dem.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_dem.t,
                                      format=argv$iff_dem.tfmt,
                                      e=argv$iff_dem.e)))
  if (is.null(raux)) 
    boom("error reading the digital elevation model file")
  rdem<-raux$stack; rm(raux)
  dem<-getValues(rdem)
  if (!rasters_match(rdem,rmaster)) {
    if (argv$iff_dem.varname_lat=="none") {
      rdem<-projectRaster(rdem,rmaster)
      dem<-getValues(rdem)
    } else {
      raux<-try(read_dotnc(nc.file=argv$iff_dem,
                       nc.varname=argv$iff_dem.varname_lat,
                       topdown=argv$iff_dem.topdown,
                       out.dim=list(ndim=argv$iff_dem.ndim,
                                    tpos=argv$iff_dem.tpos,
                                    epos=argv$iff_dem.epos,
                                    names=argv$iff_dem.names),
                       proj4=argv$iff_dem.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_dem.t,
                                      format=argv$iff_dem.tfmt,
                                      e=argv$iff_dem.e)))
      if (is.null(raux)) 
        boom("error reading digital elevation model (lat)")
      lat<-raux$data; rm(raux)
      raux<-try(read_dotnc(nc.file=argv$iff_dem,
                           nc.varname=argv$iff_dem.varname_lon,
                           topdown=argv$iff_dem.topdown,
                           out.dim=list(ndim=argv$iff_dem.ndim,
                                        tpos=argv$iff_dem.tpos,
                                        epos=argv$iff_dem.epos,
                                        names=argv$iff_dem.names),
                           proj4=argv$iff_dem.proj4,
                           nc.proj4=list(var=NULL,
                                         att=NULL),
                           selection=list(t=argv$iff_dem.t,
                                          format=argv$iff_dem.tfmt,
                                          e=argv$iff_dem.e)))
      if (is.null(raux)) 
        boom("error reading digital elevation model (lon)")
      lon<-raux$data; rm(raux)
      coord.new<-spTransform(SpatialPoints(cbind(lon,lat),
                                           proj4string=argv$iff_dem.proj4),
                             CRS(argv$grid_master.proj4))
      #
      dem<-getValues(rdem)
      rdemagg<-rasterize(coord.new,
                        aggregate(rmaster,fact=4),
                        dem)
      rdem<-mask( crop( disaggregate(rdemagg,fact=4,method="bilinear"),
                       rmaster),
                 rmaster)
      dem<-getValues(rdem)
    }
  }
}
#
#------------------------------------------------------------------------------
# read first guess
if (file.exists(argv$iff_fg)) {
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print(paste("+ first guess ",argv$iff_fg))
  }
  argv$iff_fg.epos<-set_NAs_to_NULL(argv$iff_fg.epos)
  argv$iff_fg.tpos<-set_NAs_to_NULL(argv$iff_fg.tpos)
  argv$iff_fg.e<-set_NAs_to_NULL(argv$iff_fg.e)
  if (argv$iff_fg.t=="none") argv$iff_fg.t<-nc4.getTime(argv$iff_fg)[1]
  raux<-try(read_dotnc(nc.file=argv$iff_fg,
                       nc.varname=argv$iff_fg.varname,
                       topdown=argv$iff_fg.topdown,
                       out.dim=list(ndim=argv$iff_fg.ndim,
                                    tpos=argv$iff_fg.tpos,
                                    epos=argv$iff_fg.epos,
                                    names=argv$iff_fg.names),
                       proj4=argv$iff_fg.proj4,
                       nc.proj4=list(var=NULL,
                                     att=NULL),
                       selection=list(t=argv$iff_fg.t,
                                      format=argv$iff_fg.tfmt,
                                      e=argv$iff_fg.e)))
  if (is.null(raux)) 
    boom("error reading the rescaling file")
  rfg<-raux$stack; rm(raux)
  if (!rasters_match(rfg,rmaster)) rfg<-projectRaster(rfg,rmaster)
  rfg<-mask(rfg,rmaster)
  xb0<-getValues(rfg)
  aix<-which(!is.na(xb0))
  xb<-xb0[aix]
  rm(xb0)
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
  }
}
#
#------------------------------------------------------------------------------
# Read observations
dat<-read.table(file=argv$iff_obs,
                header=T,
                sep=argv$iff_obs.sep,
                stringsAsFactors=F,
                strip.white=T)
# read x_orig,y_orig,value & set x,y 
varidxtmp<-match(c(argv$iff_obs.x,
                   argv$iff_obs.y,
                   argv$iff_obs.value),
                   names(dat))
if (any(is.na(varidxtmp))) {
  print("ERROR in the specification of the variable names")
  print(paste("        x=",argv$iff_obs.x))
  print(paste("        y=",argv$iff_obs.y))
  print(paste("    value=",argv$iff_obs.value))
  print("header of input file:")
  print(argv$argv$iff_obs)
  print(names(dat))
  quit(status=1)
}
data<-data.frame(dat[,varidxtmp])
names(data)<-c("x_orig","y_orig","value")
data$x_orig<-suppressWarnings(as.numeric(data$x_orig))
data$y_orig<-suppressWarnings(as.numeric(data$y_orig))
data$value<-suppressWarnings(as.numeric(data$value))
if (argv$iff_obs.proj4!=argv$grid_master.proj4) {
  xymaster<-spTransform(SpatialPoints(cbind(data$x_orig,data$y_orig),
                                      proj4string=CRS(argv$iff_obs.proj4)) ,
                        CRS(argv$grid_master.proj4))
  data$x<-attr(xymaster,"coords")[,1]
  data$y<-attr(xymaster,"coords")[,2]
  rm(xymaster)
} else {
  data$x<-data$x_orig
  data$y<-data$y_orig
}
# lat-lon coords are required by verif
if (argv$iff_obs.proj4!=proj4.llwgs84) {
  xyll<-spTransform(SpatialPoints(cbind(data$x_orig,data$y_orig),
                                  proj4string=CRS(argv$iff_obs.proj4)) ,
                    CRS(proj4.llwgs84))
  data$lon<-attr(xyll,"coords")[,1]
  data$lat<-attr(xyll,"coords")[,2]
  rm(xyll)
} else {
  data$lon<-data$x_orig
  data$lat<-data$y_orig
}
# read z 
if (argv$iff_obs.z!="none") {
  varidxtmp<-which(argv$iff_obs.z==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     z=",argv$iff_obs.z))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$z<-dat[,varidxtmp]
} else {
  data$z<-rep(0,length(data$x))
}
# read sourceId 
if (argv$iff_obs.sourceId!="none") {
  varidxtmp<-which(argv$iff_obs.sourceId==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     sourceId=",argv$iff_obs.sourceId))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$sourceId<-dat[,varidxtmp]
} else {
  data$sourceId<-1:length(data$x)
}
# read prId 
if (argv$iff_obs.prId!="none") {
  varidxtmp<-which(argv$iff_obs.prId==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     prId=",argv$iff_obs.prId))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$prId<-dat[,varidxtmp]
} else {
  data$prId<-rep(0,length(data$x))
}
# read dqc flag 
if (argv$iff_obs.dqc!="none") {
  varidxtmp<-which(argv$iff_obs.dqc==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     dqc=",argv$iff_obs.dqc))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$dqc<-dat[,varidxtmp]
} else {
  data$dqc<-rep(0,length(data$x))
}
# data, dataframe: x,y,z,x_orig,y_orig,value,sourceId,prId,dqc
#
if (file.exists(argv$iff_black)) {
  bstid<-read.csv(argv$iff_black,header=T,stringsAsFactors=F,strip.white=T)
} else {
  bstid<-integer(0)
}
# select only observation within the master grid
flag_in_master<-!is.na(extract(rmaster,cbind(data$x,data$y)))
flag_in_fg<-rep(T,length(data$x))
if (file.exists(argv$iff_fg)) 
  flag_in_fg<-!is.na(extract(rfg,cbind(data$x,data$y))) 
# on-the-fly dqc, used for testing
#  flag_in_fg<-!is.na(extract(rfg,cbind(data$x,data$y))) &
#              data$value > (extract(rfg,cbind(data$x,data$y))-0.2*extract(rfg,cbind(data$x,data$y)))
#CVmode
if (argv$cv_mode) {
# prId=1 MET-WMO stations
  ixcv<-which( data$dqc==0 & 
               data$prId %in% argv$prId.cv & 
               !is.na(data$value) &
               flag_in_master &
               flag_in_fg   )
  if (length(ixcv)==0) boom("ERROR \"cv_mode\" running without CV-observations ")
  VecX_cv<-data$x[ixcv]
  VecY_cv<-data$y[ixcv]
  VecXorig_cv<-data$x_orig[ixcv]
  VecYorig_cv<-data$y_orig[ixcv]
  VecLat_cv<-data$lat[ixcv]
  VecLon_cv<-data$lon[ixcv]
  VecZ_cv<-data$z[ixcv]
  VecS_cv<-data$sourceId[ixcv]
  yo_cv<-data$value[ixcv]
  if (exists("rrf")) yrf_cv<-extract(rrf,cbind(VecX_cv,VecY_cv),na.rm=T)
  data$value[ixcv]<-NA
  data$dqc[ixcv]<-999
  ncv<-length(ixcv)
}
if (any(!is.na(argv$prId.exclude))) {
  ix0<-which(data$dqc==0 & 
             !(data$prId %in% argv$prId.exclude) &
             flag_in_master &
             flag_in_fg   )
} else {
  ix0<-which(data$dqc==0 &
             flag_in_master &
             flag_in_fg   )
}
n0<-length(ix0)
# definitive station list
if (n0==0) boom("No observations. Stop here.")
VecX<-data$x[ix0]
VecY<-data$y[ix0]
VecXorig<-data$x_orig[ix0]
VecYorig<-data$y_orig[ix0]
VecLat<-data$lat[ix0]
VecLon<-data$lon[ix0]
VecZ<-data$z[ix0]
VecS<-data$sourceId[ix0]
yo<-data$value[ix0]
if (exists("rrf"))  yrf<-extract(rrf,cbind(VecX,VecY),na.rm=T)
if (exists("rlaf")) VecLaf<-extract(rlaf,cbind(VecX,VecY),na.rm=T)*argv$iff_laf.adjfact
prId<-data$prId[ix0]
ydqc.flag<-rep(0,length=n0)
rm(data)
if (!is.na(argv$rrinf)) {
  ixwet<-which(yo>=argv$rrinf)
  ixdry<-which(yo< argv$rrinf)
  nwet<-length(ixwet)
  ndry<-length(ixdry)
}
# Easy-peasy
#if (nwet==0) {
#  if (argv$verbose) print("no rain over the whole domain")
#  writeIfNoPrec()
#  if (argv$verbose) print("Success exit")
#  quit(status=0)
#}
# observation error variance correction factor
ovarc<-rep(1,n0)
if (any(!is.na(argv$ovarc.prId))) {
  for (i in 1:length(argv$ovarc.prId)) {
    if (any(prId==argv$ovarc.prId[i])) 
      ovarc[which(prId==argv$ovarc.prId[i])]<-argv$ovarc[i]
  }
}
if (argv$verbose) { 
  print("+---------------------------------------------------------------+")
  if (!is.na(argv$rrinf)) {
    print(paste("#observations (wet/dry) =",n0,"(",nwet,"/",ndry,")"))
  } else {
    print(paste("#observations =",n0))
  }
}
#
#------------------------------------------------------------------------------
# Set the OI multi-scale parameters
if (argv$mode=="OI_multiscale") {
  if (n0<100) {
    kseq<-c(2,3,4,seq(5,n0,by=5),n0)
  } else {
    kseq<-c(2,3,4,seq(5,100,by=5),seq(100,n0,by=200),n0)
  }
  kseq<-rev(unique(kseq))
  vecd_tmp<-vector()
  vecf<-vector()
  for (i in 1:length(kseq)) 
    vecd_tmp[i]<-round(dobs_fun(obs=data.frame(x=VecX,y=VecY),k=kseq[i])/1000,0)
  if (vecd_tmp[length(kseq)]>=15) {
    vecd_tmp[c(i+1,i+2,i+3)]<-c(12,8,5)
  } else if (vecd_tmp[length(kseq)]>=10) {
    vecd_tmp[c(i+1,i+3)]<-c(10,5)
  }
  vecf_tmp<-pmin(min(c(nx,ny))/5,pmax(1,round(vecd_tmp/5,0)))
  vecd<-vecd_tmp[which(!duplicated(vecf_tmp,fromLast=T))]
  kseq<-kseq[which(!duplicated(vecf_tmp,fromLast=T))]
  vecf<-round(vecd/5,0)
  vece<-rep(argv$eps2,length(vecf))
  nl<-length(vecd)
  if (length(which(vecf<=2))>0) {
    vece[which(vecf<=2)]<-rep(0.1,length(which(vecf<=2)))
  } else {
    vece[length(vece)]<-0.1
  }
  rm(vecf_tmp,vecd_tmp)
  if (argv$verbose) {
    print("+---------------------------------------------------------------+")
    print("kseq vecd vecf vece")
    print(" #    km    -    -")
    print(cbind(kseq,vecd,vecf,vece))
    print("+---------------------------------------------------------------+")
  }
}
#
#------------------------------------------------------------------------------
# compute Disth (symmetric) matrix: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
Disth<-matrix(ncol=n0,nrow=n0,data=0.)
Disth<-(outer(VecY,VecY,FUN="-")**2.+
        outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
#
#------------------------------------------------------------------------------
# ANALYSIS
if (argv$verbose) {
  print("+---------------------------------------------------------------+")
  print("Analysis")
}
#..............................................................................
# ===> OI with background  <===
if (argv$mode=="OI_firstguess") {
  yb<-extract(rfg,
              cbind(VecX,VecY),
              method="bilinear")
  D<-exp(-0.5*(Disth/argv$Dh)**2.)
  diag(D)<-diag(D)+argv$eps2
  InvD<-chol2inv(chol(D))
  if (argv$transf=="none") {
    xa<-OI_RR_fast(yo=yo,
                   yb=yb,
                   xb=xb,
                   xgrid=xgrid[aix],
                   ygrid=ygrid[aix],
                   VecX=VecX,
                   VecY=VecY,
                   Dh=argv$Dh)
  } else if (argv$transf=="Box-Cox") {
    t00<-Sys.time()
    res<-OI_RR_var(yo=boxcox(yo,argv$transf.boxcox_lambda),
                   yb=boxcox(yb,argv$transf.boxcox_lambda),
                   xb=boxcox(xb,argv$transf.boxcox_lambda),
                   gx=xgrid[aix],
                   gy=ygrid[aix],
                   ox=VecX,
                   oy=VecY,
                   Dh=argv$Dh,
                   eps2=argv$eps2)
    t11<-Sys.time()
    if (argv$verbose) print(paste("OI_RR, time=",round(t11-t00,1),attr(t11-t00,"unit")))
    t00<-Sys.time()
    xa<-apply(cbind(res$xa,sqrt(abs(res$xa_errvar))),
              MARGIN=1,
              FUN=tboxcox4pdf_apply,
                  lambda=argv$transf.boxcox_lambda,
                  brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
    t11<-Sys.time()
    if (argv$verbose) print(paste("backtransf, time=",round(t11-t00,1),attr(t11-t00,"unit")))
  } else {
    boom("transformation not defined")
  }
  ra<-rmaster
  ra[]<-NA
  ra[aix]<-xa
  rm(rfg,xb,xa,aix,xgrid,ygrid)
  if (argv$loocv_mode) {
    W<-tcrossprod((D-argv$eps2*diag(n0)),InvD)
    if (argv$transf=="none") {
      ya<-OI_RR_fast(yo=yo,
                     yb=yb,
                     xb=yb,
                     xgrid=VecX,
                     ygrid=VecY,
                     VecX=VecX,
                     VecY=VecY,
                     Dh=argv$Dh)
      yav<-yo + 1./(1.-diag(W)) * (ya-yo)
    } else if (argv$transf=="Box-Cox") {
      ya<-apply(cbind(res$ya,sqrt(abs(res$ya_errvar))),
                MARGIN=1,
                FUN=tboxcox4pdf_apply,
                    lambda=argv$transf.boxcox_lambda,
                    brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
      yavbc<-boxcox(yo,argv$transf.boxcox_lambda) + 
             1./(1.-diag(W)) * (res$ya-boxcox(yo,argv$transf.boxcox_lambda))
      # this is from Lussana et al (2010), Eq.(19) 
      yavbc_errvar<-res$o_errvar/argv$eps2*(diag(InvD)-res$o_errvar)
      yav<-apply(cbind(yavbc,
                       sqrt(abs(yavbc_errvar))),
                       MARGIN=1,
                       FUN=tboxcox4pdf_apply,
                           lambda=argv$transf.boxcox_lambda,
                           brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
    }
  # loo-crossvalidation not required
  } else {
    if (argv$transf=="none") {
      ya<-extract(ra,cbind(VecX,VecY),method="bilinear")
    } else if (argv$transf=="Box-Cox") {
      ya<-apply(cbind(res$ya,sqrt(abs(res$ya_errvar))),
                MARGIN=1,
                FUN=tboxcox4pdf_apply,
                    lambda=argv$transf.boxcox_lambda,
                    brrinf=boxcox(argv$rrinf,argv$transf.boxcox_lambda))
    }
    yav<-rep(-9999,length(ya))
  }
  if (argv$idiv_instead_of_elev) {
    if (!exists("W")) W<-tcrossprod((D-argv$eps2*diag(n0)),InvD)
    # this is the cross-validation integral data influence ("yidiv")
    elev_for_verif<-rep(1,n0) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,n0))
  } else {
    elev_for_verif<-VecZ
  }
#..............................................................................
# ===>  OI multiscale (without background)   <===
} else if (argv$mode=="OI_multiscale") {
  # multi-scale OI operates on relative anomalies (if rescaling factor is available)
  if (exists("yrf")) {
    if (any(yrf==0)) yrf[which(yrf==0)]<-1
    yo_relan<-yo/yrf
  #  arrinf<-mean(argv$rrinf/rf[mask])
  } else {
    yo_relan<-yo
  #  arrinf<-mean(argv$rrinf/rf[mask])
  }
  # multi-scale OI
  for (l in 1:nl) {
    if (argv$verbose) 
      print(paste("scale # ",formatC(l,width=2,flag="0"),
                  " of ",nl,
                  " (",formatC(vecd[l],width=4,flag="0"),"km ",
                  "fact=",formatC(vecf[l],width=3,flag="0"),")",sep=""))
    D<-exp(-0.5*(Disth/vecd[l])**2.)
    diag(D)<-diag(D)+vece[l]*ovarc
    InvD<-chol2inv(chol(D))
    if (l==nl | vecf[l]==1) {
      r<-rmaster
    } else {
      r<-aggregate(rmaster, fact=vecf[l], expand=T, na.rm=T)
    }
    zvalues.l<-getValues(r)
    storage.mode(zvalues.l)<-"numeric"
    xy.l<-xyFromCell(r,1:ncell(r))
    x.l<-sort(unique(xy.l[,1]))
    y.l<-sort(unique(xy.l[,2]),decreasing=T)
    mask.l<-which(!is.na(zvalues.l))
    zgrid.l<-zvalues.l[mask.l]
    xgrid.l<-xy.l[mask.l,1]
    ygrid.l<-xy.l[mask.l,2]
    rm(xy.l,zvalues.l)
    if (l==1) {
      yb<-rep(mean(yo_relan),length=n0)
      xb<-rep(mean(yo_relan),length=length(xgrid.l))
    } else {
      rb<-resample(ra,r,method="bilinear")
      xb<-getValues(rb)[mask.l]
      count<-0
      while (any(is.na(xb))) {
        count<-count+1
        buffer_length<-round(vecd[l]/(10-(count-1)),0)*1000
        if (!is.finite(buffer_length)) break 
        ib<-which(is.na(xb))
        aux<-extract(rb,cbind(xgrid.l[ib],ygrid.l[ib]),
                     na.rm=T,
                     buffer=buffer_length)
        for (ll in 1:length(aux)) xb[ib[ll]]<-mean(aux[[ll]],na.rm=T)
        rb[mask.l]<-xb
        rm(aux,ib)
      }
      yb<-extract(rb,cbind(VecX,VecY),method="bilinear")
#      debug_plots()
      rm(rb)
    }
    if (any(is.na(xb))) print("xb is NA")
    if (any(is.na(yb))) print("yb is NA")
    xa.l<-OI_RR_fast(yo=yo_relan,
                     yb=yb,
                     xb=xb,
                     xgrid=xgrid.l,
                     ygrid=ygrid.l,
                     VecX=VecX,
                     VecY=VecY,
                     Dh=vecd[l]) 
    xidiw.l<-OI_RR_fast(yo=yo[ixwet],
                        yb=rep(0,nwet),
                        xb=rep(0,length(xgrid.l)),
                        xgrid=xgrid.l,
                        ygrid=ygrid.l,
                        VecX=VecX[ixwet],
                        VecY=VecY[ixwet],
                        Dh=vecd[l]) 
    xidid.l<-OI_RR_fast(yo=yo[ixdry],
                        yb=rep(0,ndry),
                        xb=rep(0,length(xgrid.l)),
                        xgrid=xgrid.l,
                        ygrid=ygrid.l,
                        VecX=VecX[ixdry],
                        VecY=VecY[ixdry],
                        Dh=vecd[l]) 
    if (any(xidid.l>xidiw.l)) xa.l[which(xidid.l>xidiw.l)]<-0
  #  if (any(xidid.l<xidiw.l & xa.l<arrinf)) 
  #    xa.l[which(xidid.l<xidiw.l & xa.l<arrinf)]<-arrinf
    ra<-r
    ra[]<-NA
    ra[mask.l]<-xa.l
  } # end of multi-scale OI
  # back to precipitation values (from relative anomalies)
  if (exists("rf")) xa<-round(xa.l*rf[mask.l],2)
  ra[mask.l]<-xa
  rm(mask.l,xa)
  ya<-extract(ra,cbind(VecX,VecY),method="bilinear")
  yb<-rep(-9999,length(ya))
  yav<-rep(-9999,length(ya))
#..............................................................................
# ===>  OI two-step spatial interpolation (without background)   <===
} else if (argv$mode=="OI_twosteptemperature") {
  mask<-which(!is.na(xgrid) & 
              !is.na(ygrid) &
              !is.na(getValues(rmaster)) &
              !is.na(dem) &
              !is.na(laf))
  xgrid<-xgrid[mask]
  ygrid<-ygrid[mask]
  dem<-dem[mask]
  laf<-laf[mask]
  # superobbing for the background field
  if (argv$twostep_superobbing) {
    if (argv$verbose) {
      print("+----------------------------+")
      print("superobbing for the bg field")
      t00<-Sys.time()
    }
    # remove_me vector, to identify obs that have been used for the so
    rm_me<-vector(mode="logical",length=n0)
    rm_me[]<-F
    # identify (1x1) boxes with more than 3 obs within them
    r<-rmaster
    r[]<-1:ncell(r)
    VecI<-extract(r,cbind(VecX,VecY))
    nobs<-getValues(rasterize(cbind(VecX,VecY),r,yo,fun=function(x,...)length(x)))
    ix<-which(!is.na(nobs) & nobs>3)
    if (length(ix)>10) {
      box_ixy<-cbind(ix,xy[ix,1],xy[ix,2])
      rm(r,xy,ix)
      # initialize bg vectors
      VecX_bg<-integer(0)
      VecY_bg<-integer(0)
      VecZ_bg<-integer(0)
      yo_bg<-integer(0)
      VecLaf_bg<-integer(0)
      for (i in 1:length(box_ixy[,1])) {
        so_i<-superobs(ixy=box_ixy[i,],res=c(dx,dy),n=c(2,2))
        if (so_i$n==0) next
        VecX_bg<-c(VecX_bg,so_i$x)
        VecY_bg<-c(VecY_bg,so_i$y)
        VecZ_bg<-c(VecZ_bg,so_i$z)
        yo_bg<-c(yo_bg,so_i$yo)
        VecLaf_bg<-c(VecLaf_bg,so_i$laf)
        rm(so_i)
      }
      ix<-which(!rm_me)
      if (length(ix)>0) {
        VecX_bg<-c(VecX_bg,VecX[ix])
        VecY_bg<-c(VecY_bg,VecY[ix])
        VecZ_bg<-c(VecZ_bg,VecZ[ix])
        yo_bg<-c(yo_bg,yo[ix])
        VecLaf_bg<-c(VecLaf_bg,VecLaf[ix])
      }
      rm(ix,rm_me)
      nbg<-length(VecX_bg)
      if (argv$verbose) {
        print(paste("# observations (after superobbing) =",nbg))
        t11<-Sys.time()
        print(t11-t00)
      }
    } else {
      print(paste("not enough observations for the",
                  "superobbing to make a difference"))
      VecX_bg<-VecX
      VecY_bg<-VecY
      VecZ_bg<-VecZ
      yo_bg<-yo
      VecLaf_bg<-VecLaf
    }
  } else {
    VecX_bg<-VecX
    VecY_bg<-VecY
    VecZ_bg<-VecZ
    yo_bg<-yo
    VecLaf_bg<-VecLaf
  }
  # Regional backgrounds
  if (argv$verbose) {
    print("+----------------------------+")
    print("two-step interpolation: regional bg field")
    t00<-Sys.time()
  }
  # Spatial scale definitions. 
  #  Regional=whole domain
  #  Sub-regional (or local scale)=dozens of observations (10-100)
  #  small-scale=few observations (1-10)
  #  sub-grid scale=not observed
  # define grid used to compute local backgrounds 
  #  (grid nodes are candidates for the centroids)
  r<-raster(extent(xmn,xmx,ymn,ymx),
            ncol=argv$grid.bg[2],
            nrow=argv$grid.bg[1],
            crs=argv$grid_master.proj4)
  res<-res(r)
  mures<-round(mean(res),0)
  xy<-xyFromCell(r,1:ncell(r))
  xr<-xy[,1]
  yr<-xy[,2]
  ir<-1:ncell(r)
  r[]<-1:ncell(r)
  # attribute each station to a grid box
  VecI_bg<-extract(r,cbind(VecX_bg,VecY_bg))
  # selection of centroids
  VecI4w<-extract(rmaster,
                  cbind(xr,yr),
                  buffer=(2*argv$obs.outbuffer),
                  na.rm=T,
                  fun=mean)
  ir_in<-unique(VecI_bg[which(!is.na(VecI_bg))])
  irx<-which( (ir %in% ir_in) & (!is.na(VecI4w)))
  if (argv$verbose) {
    print(paste("the master grid has been divided in",
                argv$grid.bg[1],"x",argv$grid.bg[2],"boxes"))
    print(paste("sub-regional area extensions (length x (m),length y (m))=",
          round(res[1],0),round(res[2],0)))
    print(paste("reference (horizontal) lenght scale to weight the sub-regional backgrounds (m)=",round(mures,0)))
    print(paste("# sub-regional centroids",length(irx)))
  }
  # count the number of observations in each box
  rnobs<-rasterize(cbind(VecX_bg,VecY_bg),r,yo,fun=function(x,...)length(x))
  nr<-getValues(rnobs)
  # create the 4D array for the function call ...
  # ... via apply centroids=(xr[irx],yr[irx])
  ixyn<-cbind(ir[irx],xr[irx],yr[irx],nr[irx])
  # initialization
  # xb=background at grid points; xw=weights; xdh_oi=
  b_ok<-F
  na<--999.
  for (i in 1:10) {
    if (!argv$cv_mode) {
      xb<-xgrid
      xw<-xgrid
      xdh_oi<-xgrid
      xb[]<-na
      xw[]<-na
      xdh_oi[]<-na
    # CVmode
    } else {
      xb<-VecX_cv
      xw<-VecX_cv
      xdh_oi<-VecX_cv
      xb[]<-na
      xw[]<-na
      xdh_oi[]<-na
    }
    yb<-VecX
    yw<-VecX
    ydh_oi<-VecX
    yb[]<-na
    yw[]<-na
    ydh_oi[]<-na
    out<-apply(ixyn,
               FUN=background_incAv,
               MARGIN=1,
               nmin=argv$n.bg,
               refdist=mures, # Dh used to compute IDI
               maxboxl=(i*argv$maxboxl),
               dzmin=argv$dz.bg,
               eps2=0.1,
               closeNth=argv$nclose.bg)
    if (!(any(yb==na) | ((!argv$cv_mode) & (any(xb==na))))) {
      b_ok<-T
      print(paste("maxboxl=",(i*argv$maxboxl),"m"))
      break
    }
  }
  if (!b_ok) {
    if (any(yb==na)) 
      boom("ERROR: found NAs in background values at station locations (try to increase \"maxboxl\"")
    if (!argv$cv_mode) {
      if (any(xb==na)) boom("ERROR: found NAs in background values at grid points")
    } 
  }
  # smooth out dh
  #save.image("img0.RData")
  #r<-rmaster
  #r[mask]<-xdh_oi
  #afact<-round((mures/mean(res(rmaster)))/2.,0)
  #print(afact)
  #rf<-focal(r, w=matrix(1,21,21), fun=mean,na.rm=T)
  #xdh_oi<-getValues(rf)[mask]
  #ydh_oi_tmp<-extract(rf,cbind(VecX,VecY))
  #ydh_oi[which(!is.na(ydh_oi_tmp))]<-ydh_oi_tmp[which(!is.na(ydh_oi_tmp))]
  #rm(ydh_oi_tmp,rf,r)
  if (argv$verbose) {
    if (!argv$cv_mode) {
      print(paste("RMS(yo-yb)=", round(sqrt(mean((yo-yb)**2)),2), "degC" ))
    } else {
      print(paste("RMS(yo-yb)=", round(sqrt(mean((yo_cv-xb)**2)),2), "degC" ))
    }
    t11<-Sys.time()
    print(t11-t00)
  } 
  # Analysis, OI
  if (argv$verbose) {
    print("Analysis")
    t00<-Sys.time()
  }
  innov<-yo-yb
  # standar mode (no cv)
  if (!argv$cv_mode) {
#  return(c(xa,xidi,xav,xidiv))
    ngrid<-length(xgrid)
    length_tot<-ngrid
    xout<-apply(cbind(xgrid,ygrid,dem,laf,xdh_oi,xb,1:ngrid),
                FUN=oiIT,
                MARGIN=1,
                eps2=argv$eps2,
                dz=argv$dz,
                lafmn=argv$lafmin,
                nmaxo=argv$nmaxo,
                cv=FALSE)
    if (argv$debug) {
      save.image("img.RData")
      file<-file.path(argv$debug.dir,"deb_xa.png")
      rdeb<-rmaster
      rdeb[mask]<-xout[1,]
      writeRaster(rdeb, file, format="CDF",overwrite=T)
    }
    if (argv$verbose) print(paste("grid points, time",(Sys.time()-t00)))
    if (argv$verbose) t000<-Sys.time()
    for (i in 1:length(argv$off_grd.variables)) {
      if (!exists("r.list")) r.list<-list()
      if (argv$off_grd.variables[i]=="analysis") {
        ra<-rmaster
        ra[]<-NA
        ra[mask]<-xout[1,]
        r.list[[i]]<-matrix(data=getValues(ra),
                            ncol=length(y),
                            nrow=length(x))
      } else if (argv$off_grd.variables[i]=="background") {
        r<-rmaster
        r[]<-NA
        r[mask]<-xb
        r.list[[i]]<-matrix(data=getValues(r),
                            ncol=length(y),
                            nrow=length(x))
        rm(r)
      } else if (argv$off_grd.variables[i]=="idi") {
        r<-rmaster
        r[]<-NA
        r[mask]<-xout[2,]
        r.list[[i]]<-matrix(data=getValues(r),
                            ncol=length(y),
                            nrow=length(x))
        rm(r)
      } else if (argv$off_grd.variables[i]=="dh") {
        r<-rmaster
        r[]<-NA
        r[mask]<-xdh_oi
        r.list[[i]]<-matrix(data=getValues(r),
                            ncol=length(y),
                            nrow=length(x))
        rm(r)
      }
    }
    rm(xout)
    # station points
#  return(c(xa,xidi,xav,xidiv))
    yout<-apply(cbind(VecX,VecY,VecZ,VecLaf,ydh_oi,yb,1:n0),
                FUN=oiIT,
                MARGIN=1,
                eps2=argv$eps2,
                dz=argv$dz,
                lafmn=argv$lafmin,
                nmaxo=argv$nmaxo,
                cv=TRUE)
    #
    if (argv$debug) save.image("img1.RData")
    ya<-yout[1,]
    yav<-yout[3,]
    if (argv$verbose) {
      print(paste("  RMS(yo-ya)=", 
                  round(sqrt(mean((yo-yout[1,])**2)),2), "degC" ))
      print(paste("RMS(yo-ycva)=",
                  round(sqrt(mean((yo-yout[3,])**2)),2), "degC" ))
      t11<-Sys.time()
      print(t11-t000)
    }
  # CVmode
  } else {
    length_tot<-length(VecX_cv)
    yout<-apply(cbind(VecX_cv,VecY_cv,VecZ_cv,VecLaf_cv,xdh_oi,xb,1:length_tot),
                FUN=oiIT,
                MARGIN=1,
                eps2=argv$eps2,
                dz=argv$dz,
                lafmn=argv$lafmin,
                nmaxo=argv$nmaxo,
                cv=FALSE)
    if (argv$debug) save.image("img1.RData")
    if (argv$verbose) {
      print(paste("RMS(yo-ycva)=", round(sqrt(mean((yo_cv-yout[1,])**2)),2), "degC" ))
      t11<-Sys.time()
      print(t11-t000)
    }
  }
  elev_for_verif<-VecZ
} # end if for the selection among OIs
#
#------------------------------------------------------------------------------
# set dry regions to no-rain
jump<-T
if (!jump) {
  print("+---------------------------------------------------------------+")
  print("adjust for unlikely wet regions left by multi OI")
  Dh<-vecd[which(kseq==10 & !is.na(kseq))]
  D<-exp(-0.5*(Disth[ixwet,ixwet]/Dh)**2.)
  diag(D)<-diag(D)+ovarc[ixwet]
  InvD<-chol2inv(chol(D))
  xidiw.l<-OI_RR_fast(yo=rep(1,nwet),
                      yb=rep(0,nwet),
                      xb=rep(0,length(xgrid.l)),
                      xgrid=xgrid.l,
                      ygrid=ygrid.l,
                      VecX=VecX[ixwet],
                      VecY=VecY[ixwet],
                      Dh=Dh) 
  D<-exp(-0.5*(Disth[ixdry,ixdry]/Dh)**2.)
  diag(D)<-diag(D)+ovarc[ixdry]
  InvD<-chol2inv(chol(D))
  xidid.l<-OI_RR_fast(yo=rep(1,ndry),
                      yb=rep(0,ndry),
                      xb=rep(0,length(xgrid.l)),
                      xgrid=xgrid.l,
                      ygrid=ygrid.l,
                      VecX=VecX[ixdry],
                      VecY=VecY[ixdry],
                      Dh=Dh)
  ixd<-which(((0.6*xidid.l)>xidiw.l) & xa<1)
  if (length(ixd)>0) xa[ixd]<-0
  #debug_01_plots()
  ra[mask.l]<-xa
  rm(r,xb,yb,xa.l)
  #
  print("remove wet regions with no obs in them")
  xa_aux<-getValues(ra)
  xa_aux[which(!is.na(xa_aux) & xa_aux<argv$rrinf)]<-NA
  ra[]<-xa_aux
  rclump<-clump(ra)
  oclump<-extract(rclump,cbind(VecX[ixwet],VecY[ixwet]))
  fr<-freq(rclump)
  # remove clumps of YESprec cells less than (4x4)km^2 or not including wet obs
  ix<-which(!is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=16) | !(fr[,1] %in% oclump)) )
  xa[which(getValues(rclump)[mask.l] %in% fr[ix,1])]<-0
  rm(xa_aux,rclump,oclump,fr,ix)
  #
  print("substitute dry regions with no obs in them with wet regions")
  ra[mask.l]<-xa
  xa_aux<-getValues(ra)
  xa_aux[which(!is.na(xa_aux) & xa_aux>=argv$rrinf)]<-NA
  xa_aux[which(!is.na(xa_aux))]<-1
  ra[]<-xa_aux
  rclump<-clump(ra)
  oclump<-extract(rclump,cbind(VecX[ixdry],VecY[ixdry]))
  fr<-freq(rclump)
  # remove clumps of YESprec cells less than (4x4)km^2 or not including wet obs
  ix<-which(!is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=16) | !(fr[,1] %in% oclump)) )
  xa[which((getValues(rclump)[mask.l] %in% fr[ix,1]) & (xidiw.l>xidid.l))]<-argv$rrinf
  rm(xa_aux,rclump,oclump,fr,ix,xidiw.l,xidid.l)
  ## final field
}
# 
#CVmode
if (argv$cv_mode) {
  ya_cv<-extract(ra,cbind(VecX_cv,VecY_cv),method="bilinear")
  if (argv$idiv_instead_of_elev) {
    Disth<-matrix(ncol=ncv,nrow=ncv,data=0.)
    Disth<-(outer(VecY_cv,VecY_cv,FUN="-")**2.+
            outer(VecX_cv,VecX_cv,FUN="-")**2.)**0.5/1000.
    D<-exp(-0.5*(Disth/argv$Dh)**2.)
    diag(D)<-diag(D)+argv$eps2
    InvD<-chol2inv(chol(D))
    W<-tcrossprod((D-argv$eps2*diag(ncv)),InvD)
    # this is the cross-validation integral data influence ("yidiv")
    elev_for_verif_cv<-rep(1,ncv) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,ncv))
  } else {
    elev_for_verif_cv<-VecZ_cv
  }
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (argv$verbose) print("++ Output")
# Station Points
# CVmode
if (argv$cv_mode) {
  cat(paste0("# variable: ",argv$verif_var,"\n"),file=argv$off_cvver,append=F)
  cat(paste0("# units: $",argv$verif_units,"$\n"),file=argv$off_cvver,append=T)
  cat("unixtime     leadtime location  lat     lon      altitude obs      fcst\n",
      file=argv$off_cvver,append=T)
  ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(yav))
  cat(paste( as.numeric(as.POSIXct(argv$date_out, format=argv$date_out_fmt))," ",
             rep(0,length(VecS_cv))," ",
             VecS_cv," ",
             VecLat_cv," ",
             VecLon_cv," ",
             elev_for_verif_cv," ",
             round(yo_cv,2)," ",
             round(ya_cv,2),"\n",
             sep=""),
      file=argv$off_cvver,append=T)
  print(paste("data saved on file",argv$off_cvver))
  cat("date;sourceId;x;y;z;yo;yb;ya;yav;dqc;\n",
      file=argv$off_cvstn,append=F)
  cat(paste(argv$date_out,
            formatC(VecS_cv,format="f",digits=0),
            formatC(VecX_cv,format="f",digits=0),
            formatC(VecY_cv,format="f",digits=0),
            formatC(VecZ_cv,format="f",digits=0),
            formatC(yo_cv,format="f",digits=1),
            formatC(yb_cv,format="f",digits=1),
            formatC(ya_cv,format="f",digits=1),
            formatC(ya_cv,format="f",digits=1),
            rep(0,length(VecS_cv)),
            "\n",sep=";"),
      file=argv$off_cvstn,append=T)
  print(paste("data saved on file",argv$off_cvstn))
# non-CVmode
} else {
  cat(paste0("# variable: ",argv$verif_var,"\n"),file=argv$off_ver,append=F)
  cat(paste0("# units: $",argv$verif_units,"$\n"),file=argv$off_ver,append=T)
  cat("unixtime     leadtime location  lat     lon      altitude obs      fcst\n",
      file=argv$off_ver,append=T)
  if (argv$off_ver_fg!="none" & argv$off_ver_cv!="none") {
    ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(ya) & !is.na(yb) & !is.na(yav))
  } else if (argv$off_ver_fg!="none") {
    ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(ya) & !is.na(yb))
  } else if (argv$off_ver_cv!="none") {
    ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(ya) & !is.na(yav))
  } else {
    ix<-which(ydqc.flag<=0 & !is.na(yo) & !is.na(ya))
  }
  cat(paste( as.numeric(as.POSIXct(argv$date_out, format=argv$date_out_fmt))," ",
             rep(0,length(ix))," ",
             VecS[ix]," ",
             VecLat[ix]," ",
             VecLon[ix]," ",
             elev_for_verif[ix]," ",
             round(yo[ix],1)," ",
             round(ya[ix],1),"\n",
             sep=""),
      file=argv$off_ver,append=T)
  print(paste("data saved on file",argv$off_ver))
  if (argv$off_ver_fg!="none") {
    cat(paste0("# variable: ",argv$verif_var,"\n"),file=argv$off_ver_fg,append=F)
    cat(paste0("# units: $",argv$verif_units,"$\n"),file=argv$off_ver_fg,append=T)
    cat("unixtime     leadtime location  lat     lon      altitude obs      fcst\n",
        file=argv$off_ver_fg,append=T)
    cat(paste( as.numeric(as.POSIXct(argv$date_out, format=argv$date_out_fmt))," ",
               rep(0,length(ix))," ",
               VecS[ix]," ",
               VecLat[ix]," ",
               VecLon[ix]," ",
               elev_for_verif[ix]," ",
               round(yo[ix],1)," ",
               round(yb[ix],1),"\n",
               sep=""),
        file=argv$off_ver_fg,append=T)
    print(paste("data saved on file",argv$off_ver_fg))
  }
  if (argv$off_ver_cv!="none") {
    cat(paste0("# variable: ",argv$verif_var,"\n"),file=argv$off_ver_cv,append=F)
    cat(paste0("# units: $",argv$verif_units,"$\n"),file=argv$off_ver_cv,append=T)
    cat("unixtime     leadtime location  lat     lon      altitude obs      fcst\n",
        file=argv$off_ver_cv,append=T)
    cat(paste( as.numeric(as.POSIXct(argv$date_out, format=argv$date_out_fmt))," ",
               rep(0,length(ix))," ",
               VecS[ix]," ",
               VecLat[ix]," ",
               VecLon[ix]," ",
               elev_for_verif[ix]," ",
               round(yo[ix],1)," ",
               round(yav[ix],1),"\n",
               sep=""),
        file=argv$off_ver_cv,append=T)
    print(paste("data saved on file",argv$off_ver_cv))
  }
  cat("date;sourceId;x;y;z;yo;yb;ya;yav;dqc;\n",
      file=argv$off_stn,append=F)
  digxy<-0
  if (argv$grid_master.proj4==proj4.llwgs84) digxy<-6
  cat(paste(argv$date_out,
            formatC(VecS,format="f",digits=0),
            formatC(VecX,format="f",digits=digxy),
            formatC(VecY,format="f",digits=digxy),
            formatC(VecZ,format="f",digits=0),
            formatC(yo,format="f",digits=1),
            formatC(yb,format="f",digits=1),
            formatC(ya,format="f",digits=1),
            formatC(yav,format="f",digits=1),
            formatC(ydqc.flag,format="f",digits=0),
            "\n",sep=";"),
      file=argv$off_stn,append=T)
  print(paste("data saved on file",argv$off_stn))
  # grid
  if (!exists("r.list")) {
    r.list<-list()
    r.list[[1]]<-matrix(data=getValues(ra),
                        ncol=length(y),
                        nrow=length(x))
  }
  # define time for output
  tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                    format="%Y%m%d%H%M",tz="GMT")
  time_bnds<-array(format(rev(seq(strptime(argv$date_out,argv$date_out_fmt),
                                           length=2,by=argv$time_bnds_string)),
                   format="%Y%m%d%H%M",tz="GMT"),dim=c(1,2))
  out<-write_dotnc(grid.list=r.list,
                   file.name=argv$off_grd,
                   grid.type=argv$off_grd.grid,
                   x=x,
                   y=y,
                   var.name=argv$off_grd.varname,
                   var.longname=argv$off_grd.varlongname,
                   var.standardname=argv$off_grd.varstandardname,
                   var.version=argv$off_grd.varversion,
                   var.unit=argv$off_grd.varunit,
                   times=tstamp_nc,
                   times.unit=argv$off_grd.timesunit,
                   reference=argv$off_grd.reference,
                   proj4.string=argv$grid_master.proj4,
                   lonlat.out=argv$off_grd.write_lonlat,
                   round.dig=argv$off_grd.diground,
                   summary=argv$off_grd.summary,
                   source.string=argv$off_grd.sourcestring,
                   title=argv$off_grd.title,
                   comment=argv$off_grd.comment,
                   atts.var.add=NULL,
                   var.cell_methods=argv$off_grd.cell_methods,
                   time_bnds=time_bnds,
                   cf_1.7=T)
  print(paste("data saved on file",argv$off_grd))
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# normal exit
t1<-Sys.time()
print(paste("Normal exit, time=",round(t1-t0,1),attr(t1-t0,"unit")))
quit(status=0)
