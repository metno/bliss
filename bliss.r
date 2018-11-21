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
    if (length(x)>1) {
      ix<-which(is.na(x))
      x[ix]<-NULL
    } else {
      if (is.na(x)) x<-NULL
    }
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
  # maxboxl/4.5 is the max limit for dh
  # NOTE: maxboxl/4.5, seNorge_2018 (maxboxl=250km so maxboxl/4.5=55km)
  dh_oi<-min((maxboxl/4.5),mean(findRow(x=disth,n=(length(x_aux)-closeNth+1))))
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

#---------------------------------------------------------------------------------------------------
#+ 
obsop_LapseRateConst<-function(s,m,xm,ym,zm,xs,ys,zs,
                               ifield,oval,mMinElevDiff,
                               mMinGradient,mMaxGradient,
                               mSearchRadius) {
  .C("obsop_LapseRateConst",s=as.integer(s),
                            m=as.integer(m),
                            xm=as.double(xm),
                            ym=as.double(ym),
                            zm=as.double(zm),
                            xs=as.double(xs),
                            ys=as.double(ys),
                            zs=as.double(zs),
                            ifield=as.double(ifield),
                            oval=as.double(oval),
                            mMinElevDiff=as.double(mMinElevDiff),
                            mMinGradient=as.double(mMinGradient),
                            mMaxGradient=as.double(mMaxGradient),
                            mSearchRadius=as.double(mSearchRadius))
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
                  help="statistical interpolation scheme (\"OI_multiscale\",\"OI_firstguess\",\"OI_twosteptemperature\",\"hyletkf\",\"letkf\")",
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
                  help="maximum length (m) of the box used to define a sub-region",
                  type="numeric",
                  default=250000)
p <- add_argument(p, "--obs.outbuffer",
                  help="distance (m) defining the \"buffer\" region outside the masked region where to consider the observation, so to reduce border effects",
                  type="numeric",
                  default=50000)
#------------------------------------------------------------------------------
# hyletkf
p <- add_argument(p, "--hyletkf.eps2_prec_default",
                  help="LETKF default value for the ratio between obs error variance and backg error variance in case of precipitation",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--hyletkf.eps2_noprec_default",
                  help="LETKF default value for the ratio between obs error variance and backg error variance in case of no-precipitation",
                  type="numeric",
                  default=0.05)
p <- add_argument(p, "--hyletkf.Dh",
                  help="horizontal de-correlation length for the LETKF localization",
                  type="numeric",
                  default=10)
p <- add_argument(p, "--hyletkf.pmax",
                  help="maximum number of observations in the neighbourhood of a gridpoint for LETKF",
                  type="numeric",
                  default=200)
p <- add_argument(p, "--hyletkf.sigma2_min",
                  help="minimum allowed background variance (in the transformed space) for LETKF",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--hyletkf.eps2_prId",
                  help="provider identifier corresponding to the specified ratio between obs error variance and backg error variance (LETKF)",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--hyletkf.eps2",
                  help="observation-provider dependent ratio between obs error variance and backg error variance (LETKF)",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--hyletkf.Dh_oi",
                  help="horizontal de-correlation length for the OI step of the HyLETKF (km)",
                  type="numeric",
                  default=10)
p <- add_argument(p, "--hyletkf.eps2_oi",
                  help="ratio between obs error variance and backg error variance for the OI step of the HyLETKF",
                  type="numeric",
                  default=0.1)
p <- add_argument(p, "--hyletkf.rloc_min",
                  help="use an observation in LETKF only if the localization weight is greater than the specified minimum value",
                  type="numeric",
                  default=0.0013)
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
p <- add_argument(p, "--iff_dem",
                  help="full file name for the digital elevation model (nc)",
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
p <- add_argument(p, "--off_vernc_a",
                  help="full file name for output at station locations, verif-type analyses versus observations (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_vernc_b",
                  help="full file name for output at station locations, verif-type background versus observations (nc)",
                  type="character",
                  default="none")
# output cv-mode
p <- add_argument(p, "--off_cvstn",
                  help="full file name for output at station locations (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_cvver",
                  help="full file name for output at station location, verif-type (txt)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_cvvernc_a",
                  help="full file name for output at station locations, verif-type analyses versus observations (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_cvvernc_b",
                  help="full file name for output at station locations, verif-type background versus observations (nc)",
                  type="character",
                  default="none")
p <- add_argument(p, "--off_vernc.varname",
                  help="name of the verif variable for off_vernc files",
                  type="character",
                  default="Precipitation")
p <- add_argument(p, "--off_vernc.stdvarname",
                  help="standard name of the verif variable for off_vernc files",
                  type="character",
                  default="precipitation_amount")
p <- add_argument(p, "--off_vernc.varunits",
                  help="units of the verif variable for off_vernc files",
                  type="character",
                  default="mm")
p <- add_argument(p, "--off_vernc.thresholds",
                  help="thresholds for the off_vernc files",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
p <- add_argument(p, "--off_vernc.quantile",
                  help="quantiles for the off_vernc files",
                  type="numeric",
                  nargs=Inf,
                  default=NULL)
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
p <- add_argument(p, "--iff_rf.adjfact",
                  help="correction factor",
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_rf.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
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
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_laf.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
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
p <- add_argument(p, "--iff_dem.adjfact",
                  help="correction factor",
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_dem.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
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
p <- add_argument(p, "--iff_fg.adjfact",
                  help="correction factor",
                  type="character",
                  default="1")
p <- add_argument(p, "--iff_fg.adjval",
                  help="adjustment value",
                  type="character",
                  default="0")
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
} else if (argv$mode=="hyletkf") {
  print("Hybrid Local Ensemble Transform Kalman Filter")
} else if (argv$mode=="letkf") {
  print("Local Ensemble Transform Kalman Filter")
} else {
  boom("error statistical interpolation scheme undefined")
}
# define/check paths and load external functions
if ( !(file.exists(argv$path2src)) ) 
  ext<-boom("path not found")
#
argv$iff_rf.adjfact<-as.numeric(gsub("_","-",argv$iff_rf.adjfact))
argv$iff_rf.adjval<-as.numeric(gsub("_","-",argv$iff_rf.adjval))
argv$iff_laf.adjfact<-as.numeric(gsub("_","-",argv$iff_laf.adjfact))
argv$iff_laf.adjval<-as.numeric(gsub("_","-",argv$iff_laf.adjval))
argv$iff_dem.adjfact<-as.numeric(gsub("_","-",argv$iff_dem.adjfact))
argv$iff_dem.adjval<-as.numeric(gsub("_","-",argv$iff_dem.adjval))
argv$iff_fg.adjfact<-as.numeric(gsub("_","-",argv$iff_fg.adjfact))
argv$iff_fg.adjval<-as.numeric(gsub("_","-",argv$iff_fg.adjval))
# load external C functions
dyn.load(file.path(argv$path2src,"oi_rr_first.so"))
dyn.load(file.path(argv$path2src,"oi_rr_fast.so"))
dyn.load(file.path(argv$path2src,"oi_rr_var.so"))
dyn.load(file.path(argv$path2src,"oi_t_xb_upd.so"))
dyn.load(file.path(argv$path2src,"obsop_LapseRateConst.so"))
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
  if (argv$iff_rf.adjfact!=1 | argv$iff_rf.adjval!=0) 
    raux$stack[]<-getValues(raux$stack)*argv$iff_rf.adjfact+argv$iff_rf.adjval
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
  if (argv$iff_laf.adjfact!=1 | argv$iff_laf.adjval!=0) 
    raux$stack[]<-getValues(raux$stack)*argv$iff_laf.adjfact+argv$iff_laf.adjval
  rlaf<-raux$stack; rm(raux)
  laf<-getValues(rlaf)
  if (!rasters_match(rlaf,rmaster)) {
    if (argv$iff_laf.varname_lat=="none") {
      rlaf<-projectRaster(rlaf,rmaster)
      laf<-getValues(rlaf)
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
      laf<-getValues(rlaf)
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
  if (argv$iff_dem.adjfact!=1 | argv$iff_dem.adjval!=0) 
    raux$stack[]<-getValues(raux$stack)*argv$iff_dem.adjfact+argv$iff_dem.adjval
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
  # First-guess is not an ensemble
  if (is.null(argv$iff_fg.epos)) {
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
                                        e=argv$iff_fg.e[e])))
    if (is.null(raux)) 
      boom("error reading the first-guess file")
    if (argv$iff_fg.adjfact!=1 | argv$iff_fg.adjval!=0) 
      raux$stack[]<-getValues(raux$stack)*argv$iff_fg.adjfact+argv$iff_fg.adjval
    rfg<-raux$stack; rm(raux)
    if (!rasters_match(rfg,rmaster)) rfg<-projectRaster(rfg,rmaster)
    rfg<-mask(rfg,rmaster)
    xb0<-getValues(rfg)
    aix<-which(!is.na(xb0))
    xb<-xb0[aix]
    rm(xb0)
  # First-guess is an ensemble
  } else {
    nens<-0 # number of ensemble members having at least one value different from NA
    for (e in 1:length(argv$iff_fg.e)) {
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
                                          e=argv$iff_fg.e[e])))
      if (is.null(raux)) 
        boom("error reading the first-guess file")
      if (argv$iff_fg.adjfact!=1 | argv$iff_fg.adjval!=0) 
        raux$stack[]<-getValues(raux$stack)*argv$iff_fg.adjfact+argv$iff_fg.adjval
      rfg<-raux$stack; rm(raux)
      if (!rasters_match(rfg,rmaster)) rfg<-projectRaster(rfg,rmaster)
      rfg<-mask(rfg,rmaster)
      xb0<-getValues(rfg)
      aix<-which(!is.na(xb0))
      if (length(aix)>0) {
        if (!exists("xb1")) xb1<-array(data=NA,
                                       dim=c(length(aix),length(argv$iff_fg.e)))
        if (length(aix)!=dim(xb1)[1]) boom("ERROR while reading the background file")
        xb1[,e]<-xb0[aix]
        nens<-nens+1
        if (!exists("valens")) valens<-integer(0)
        valens<-c(valens,argv$iff_fg.e[e])
      }
      rm(xb0)
    }
    if (nens==0) boom("ERROR while reading the background file")
    xb<-array(data=NA,dim=c(dim(xb1)[1],nens))
    e<-0
    for (i in 1:length(argv$iff_fg.e)) {
      if (any(!is.na(xb1[,i]))) {
        e<-e+1
        xb[,e]<-xb1[,i]
      }
    }
    rm(xb1)
  }
  if (argv$verbose) {
    print(paste("# grid points (not NAs)=",length(aix)))
    print("+...............................................................+")
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
  if (!is.na(argv$rrinf)) {
    ixwet_cv<-which(yo_cv>=argv$rrinf)
    ixdry_cv<-which(yo_cv< argv$rrinf)
    nwet_cv<-length(ixwet_cv)
    ndry_cv<-length(ixdry_cv)
  }
}
if (any(!is.na(argv$prId.exclude))) {
  ix0<-which(data$dqc==0 & 
             !(data$prId %in% argv$prId.exclude) &
             flag_in_master &
             flag_in_fg &
             !is.na(data$value) &
             !is.nan(data$value) )
} else {
  ix0<-which(data$dqc==0 &
             flag_in_master &
             flag_in_fg  &
             !is.na(data$value) &
             !is.nan(data$value) )
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
if (exists("rlaf")) VecLaf<-extract(rlaf,cbind(VecX,VecY),na.rm=T)
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
    if (argv$cv_mode) {
      print(paste("#cv-observations (wet/dry) =",ncv,"(",nwet_cv,"/",ndry_cv,")"))
    }
  } else {
    print(paste("#observations =",n0))
  }
  print("+...............................................................+")
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
    print(paste("reference (horizontal) length scale to weight the sub-regional backgrounds (m)=",round(mures,0)))
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
# END of OI two-step spatial interpolation (without background)
#..............................................................................
# ===>  Hybrid Local Ensemble Transform Kalman Filter  <===
} else if (argv$mode=="hyletkf") {
# xb, aix,nens
# yo, VecX,VecY, prId
  # save original vectors for future use
  VecX_orig<-VecX
  VecY_orig<-VecY
  prId_orig<-prId
  yo_orig<-yo
  xb_orig<-xb
  if (argv$cv_mode) {
    xgrid<-VecX_cv
    ygrid<-VecY_cv
  } else {
    xgrid_orig<-xgrid
    ygrid_orig<-ygrid
    # define grid-space
    # note: xb is defined only over "aix" points
    xgrid<-xgrid[aix]
    ygrid<-ygrid[aix]
  }
  # define obs-space, background
  r<-rmaster
  yb0<-array(data=NA,dim=c(n0,nens))
  for (e in 1:nens) {
    r[]<-NA
    r[aix]<-xb[,e]
    yb0[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
    auxx<-which(yb0[,e]<0 | is.na(yb0[,e]) | is.nan(yb0[,e]))
    if (length(auxx)>0) yb0[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
    rm(auxx)
  }
  rm(r)
  ix<-which(apply(yb0,MAR=1,
                  FUN=function(x){length(which(!is.na(x) & !is.nan(x)))})
            ==nens)
  n1<-length(ix)
  yb<-array(data=NA,dim=c(n1,nens))
  for (e in 1:nens) yb[,e]<-yb0[ix,e]
  rm(yb0)
  # background, probability of rain
  yb_pwet<-apply(yb,MAR=1,FUN=function(x){length(which(x>=argv$rrinf))/nens})
  # define obs-space, observations
  ix_orig<-ix
  yo<-yo[ix]
  VecX<-VecX[ix]
  VecY<-VecY[ix]
  prId<-prId[ix] 
  rm(ix)
  # yo_wet, 1=wet 0=dry
  yo_wet<-yo
  yo_wet[]<-0
  yo_wet[which(yo>argv$rrinf)]<-1
  # define the inverse of the diagonal observation error matrix
  diagRinv<-rep(argv$hyletkf.eps2_prec_default,n1)
  for (i in 1:length(argv$hyletkf.eps2_prId)) {
    ix<-which(prId==argv$hyletkf.eps2_prId[i])
    if (length(ix)>0) diagRinv[ix]<-argv$hyletkf.eps2[i]
    rm(ix)
  }
  diagRinv[which(yo<argv$rrinf)]<-argv$hyletkf.eps2_noprec_default
  if (any(diagRinv==0 | is.na(diagRinv) | is.nan(diagRinv))) 
    boom("ERROR in the definition of Rinv")
  diagRinv<-1./diagRinv
  # if CVmode, then adjust the background
  if (argv$cv_mode) {
    xb<-array(data=NA,dim=c(ncv,nens))
    r<-rmaster
    for (e in 1:nens) {
      r[]<-NA
      r[aix]<-xb_orig[,e]
      xb[,e]<-extract(r,cbind(VecX_cv,VecY_cv),method="bilinear")
      auxx<-which(xb[,e]<0 | is.na(xb[,e]) | is.nan(xb[,e]))
      if (length(auxx)>0) 
        xb[auxx,e]<-extract(r,cbind(VecX_cv[auxx],VecY_cv[auxx]))
      rm(auxx)
    }
    rm(r)
  }
  # Gaussian anamorphosis  
  if (argv$transf=="Box-Cox") {
    yo<-boxcox(yo,argv$transf.boxcox_lambda)
    yb<-boxcox(yb,argv$transf.boxcox_lambda)
    xb<-boxcox(xb,argv$transf.boxcox_lambda)
  }
  # ensemble mean
  xbm<-apply(xb,MAR=1,FUN=mean)
  ybm<-apply(yb,MAR=1,FUN=mean)
  # ensemble perturbations
  Xb<-xb-xbm
  Yb<-yb-ybm
  # ensemble variances
  ybvar<-apply(Yb,MAR=1,FUN=mean)**2
  # intialize analysis vector
  xa<-xb
  xa[]<-NA
  # letkf: loop over gridpoints
  if (argv$verbose) t00<-Sys.time()
  for (i in 1:length(xgrid)) {
    dist<-(((xgrid[i]-VecX)**2+(ygrid[i]-VecY)**2 )**0.5)/1000.
    rloc<-exp(-0.5*(dist**2./argv$hyletkf.Dh**2))
    sel<-which(rloc>argv$hyletkf.rloc_min)
    # i-th gridpoint analysis 
    if (length(sel)>0) {
      if (length(sel)>argv$hyletkf.pmax) sel<-order(dist)[1:argv$hyletkf.pmax]
      sel_wet<-sel[which(yb_pwet[sel]>=0.5)]
      if (length(sel_wet)>0) {
        sigma2<-max(c(mean(ybvar[sel]),argv$hyletkf.sigma2_min))
      } else {
        sigma2<-argv$hyletkf.sigma2_min
      }
      Yb.i<-Yb[sel,,drop=F]
      d.i<-yo[sel]-ybm[sel]
      C.i<-t(1./sigma2*diagRinv[sel]*rloc[sel]*Yb.i)
      C1.i<-crossprod(t(C.i),Yb.i)
      Cd.i<-crossprod(t(C.i),d.i)
      rm(C.i)
      diag(C1.i)<-diag(C1.i) + (nens-1)
      Pa.i<-chol2inv(chol(C1.i))
      rm(C1.i)
      # Pa*(k-1) is the matrix for which we want to square root:
      a.eig <- eigen(Pa.i*(nens-1),symmetric=T)
      Wa <- tcrossprod(
             tcrossprod( a.eig$vectors, diag(sqrt(a.eig$values))),
             t(solve(a.eig$vectors)) )
      waa<-crossprod(Pa.i, Cd.i )
      rm(Cd.i)
      W<-Wa+as.vector(waa)
      xa[i,]<-xbm[i]+Xb[i,] %*% W
      rm(W,Wa,waa,a.eig)
    # i-th gridpoint is isolated
    } else {
      xa[i,]<-xb[i,]
    }
  } # end: letkf loop over gridpoints
  if (argv$verbose) {
    t11<-Sys.time()
    print(paste("letkf time=",round(t11-t00,1),attr(t11-t00,"unit")))
  }
  # OI
  sel_oi<-which(yb_pwet<0.5 & yo_wet==1)  
  if (length(sel_oi)>0) {
    if (argv$verbose) t00<-Sys.time()
    noi<-length(sel_oi)
    Disth<-matrix(ncol=noi,nrow=noi,data=0.)
    Disth<-(outer(VecY[sel_oi],VecY[sel_oi],FUN="-")**2.+
            outer(VecX[sel_oi],VecX[sel_oi],FUN="-")**2.)**0.5/1000.
    D<-exp(-0.5*(Disth/argv$hyletkf.Dh_oi)**2.)
    diag(D)<-diag(D)+argv$hyletkf.eps2_oi
    InvD<-chol2inv(chol(D))
    xam<-OI_RR_fast(yo=yo[sel_oi],
                    yb=ybm[sel_oi],
                    xb=rowMeans(xa),
                    xgrid=xgrid,
                    ygrid=ygrid,
                    VecX=VecX[sel_oi],
                    VecY=VecY[sel_oi],
                    Dh=argv$hyletkf.Dh_oi) 
    xa<-xa-rowMeans(xa)+xam
    rm(xam)
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("oi time=",round(t11-t00,1),
                             attr(t11-t00,"unit")))
    }
  }
  # Back-transformation
  if (argv$transf=="Box-Cox") {
    xa<-tboxcox(xa,argv$transf.boxcox_lambda)
    xb<-tboxcox(xb,argv$transf.boxcox_lambda)
  }
  if (!argv$cv_mode & 
      ("idi" %in% argv$off_grd.variables | 
       argv$idiv_instead_of_elev)) {
    if (argv$verbose) t00<-Sys.time()
    D<-exp(-0.5*((outer(VecY,VecY,FUN="-")**2.+
                  outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
                 /argv$hyletkf.Dh_oi)**2.)
    diag(D)<-diag(D)+argv$hyletkf.eps2_oi
    InvD<-chol2inv(chol(D))
    rm(D)
    xidi<-OI_RR_fast(yo=rep(1,length(VecX)),
                     yb=rep(0,length(VecX)),
                     xb=rep(0,length(xgrid)),
                     xgrid=xgrid,
                     ygrid=ygrid,
                     VecX=VecX,
                     VecY=VecY,
                     Dh=argv$hyletkf.Dh_oi)
    if (argv$idiv_instead_of_elev) {
      G<-exp(-0.5*((outer(VecY_orig,VecY,FUN="-")**2.+
                    outer(VecX_orig,VecX,FUN="-")**2.)**0.5/1000.
                    /argv$hyletkf.Dh_oi)**2.)
      W<-tcrossprod(G,InvD)
      rm(G,InvD)
      # this is the cross-validation integral data influence ("yidiv")
      elev_for_verif<-rep(1,n0) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,n0))
      rm(W)
    }
    # InvD could be used in the following
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("idi time=",round(t11-t00,1),
                              attr(t11-t00,"unit")))
    }
  }
  if (!argv$idiv_instead_of_elev) elev_for_verif<-VecZ
  # back to the original vectors
  yo<-yo_orig
  prId<-prId_orig
  VecX<-VecX_orig
  VecY<-VecY_orig
  rm(yo_orig,prId_orig,VecX_orig,VecY_orig)
  # prepare point / gridded output
  if (!argv$cv_mode) {
    ya<-array(data=NA,dim=c(n0,length(argv$iff_fg.e)))
    yb<-array(data=NA,dim=c(n0,length(argv$iff_fg.e)))
    for (j in 1:length(argv$off_grd.variables)) {
      if (!exists("r.list")) r.list<-list()
      if (argv$off_grd.variables[j]=="analysis") {
        r<-rmaster
        r[]<-NA
        grid<-array(data=NA,dim=c(length(x),length(y),length(argv$iff_fg.e),1))
        for (e in 1:length(argv$iff_fg.e)) {
          if (argv$iff_fg.e[e] %in% valens) {
            i<-which(valens==argv$iff_fg.e[e])
            r[aix]<-xa[,i]
            grid[,,e,1]<-matrix(data=getValues(r),
                                ncol=length(y),
                                nrow=length(x))
            ya[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
            auxx<-which(ya[,e]<0 | is.na(ya[,e]) | is.nan(ya[,e]))
            if (length(auxx)>0) 
              ya[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
            rm(auxx)
          }
        }
        rm(r)
        r.list[[j]]<-grid
        rm(grid)
      } else if (argv$off_grd.variables[j]=="background") {
        r<-rmaster
        r[]<-NA
        grid<-array(data=NA,dim=c(length(x),length(y),length(argv$iff_fg.e),1))
        for (e in 1:length(argv$iff_fg.e)) {
          if (argv$iff_fg.e[e] %in% valens) {
            i<-which(valens==argv$iff_fg.e[e])
            r[aix]<-xb[,i]
            grid[,,e,1]<-matrix(data=getValues(r),
                                ncol=length(y),
                                nrow=length(x))
            yb[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
            auxx<-which(yb[,e]<0 | is.na(yb[,e]) | is.nan(yb[,e]))
            if (length(auxx)>0) 
              yb[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
            rm(auxx)
          }
        }
        rm(r)
        r.list[[j]]<-grid
        rm(grid)
      } else if (argv$off_grd.variables[j]=="idi") {
        r<-rmaster
        r[]<-NA
        r[aix]<-xidi
        r.list[[j]]<-matrix(data=getValues(r),
                            ncol=length(y),
                            nrow=length(x))
        rm(r)
      }
    }
  } # end prepare for gridded output
# END of Hybrid Local Ensemble Transform Kalman Filter
#..............................................................................
# ===>  Local Ensemble Transform Kalman Filter  <===
} else if (argv$mode=="letkf") {
  if (argv$verbose) t00<-Sys.time()
  # save original vectors for future use
  VecX_orig<-VecX
  VecY_orig<-VecY
  VecZ_orig<-VecZ
  prId_orig<-prId
  yo_orig<-yo
  xb_orig<-xb
  aix<-which(!is.na(xgrid) & !is.na(ygrid) & !is.na(dem))
  if (argv$cv_mode) {
    xgrid<-VecX_cv
    ygrid<-VecY_cv
  } else {
    xgrid_orig<-xgrid
    ygrid_orig<-ygrid
    dem_orig<-dem
    # define grid-space
    # note: xb is defined only over "aix" points
    xgrid<-xgrid[aix]
    ygrid<-ygrid[aix]
    dem<-dem[aix]
  }
  # define obs-space, background
  r<-rmaster
  yb0<-array(data=NA,dim=c(n0,nens))
  argv$letkf.obsop<-"obsop_LapseRateConst"
  argv$letkf.obsop.mMinElevDiff<-30
  argv$letkf.obsop.mMinGradient<-(-0.01)
  argv$letkf.obsop.mMaxGradient<-0.01
  argv$letkf.obsop.mSearchRadius<-10
  for (e in 1:nens) {
    if (argv$letkf.obsop=="obsop_LapseRateConst") {
      yb0[,e]<-obsop_LapseRateConst(s=n0,
                                    m=length(aix),
                                    xm=xgrid_orig[aix],
                                    ym=ygrid_orig[aix],
                                    zm=dem_orig[aix],
                                    xs=VecX,
                                    ys=VecY,
                                    zs=VecZ,
                                    ifield=xb[,e],
                                    oval=rep(0,n0),
                                    mMinElevDiff=argv$letkf.obsop.mMinElevDiff,
                                    mMinGradient=argv$letkf.obsop.mMinGradient,
                                    mMaxGradient=argv$letkf.obsop.mMaxGradient,
                                    mSearchRadius=argv$letkf.obsop.mSearchRadius)$oval
    }
  }
  ix<-which(apply(yb0,MAR=1,
                  FUN=function(x){length(which(!is.na(x) & !is.nan(x)))})
            ==nens)
  n1<-length(ix)
  yb<-array(data=NA,dim=c(n1,nens))
  for (e in 1:nens) yb[,e]<-yb0[ix,e]
  rm(yb0)
  # define obs-space, observations
  ix_orig<-ix
  yo<-yo[ix]
  VecX<-VecX[ix]
  VecY<-VecY[ix]
  VecZ<-VecZ[ix]
  prId<-prId[ix] 
  rm(ix)
  # define the inverse of the diagonal observation error matrix
  argv$letkf.eps2_default<-0.5
  argv$letkf.eps2_prId<-NA
  diagRinv<-rep(argv$letkf.eps2_default,n1)
  for (i in 1:length(argv$letkf.eps2_prId)) {
    ix<-which(prId==argv$letkf.eps2_prId[i])
    if (length(ix)>0) diagRinv[ix]<-argv$letkf.eps2[i]
    rm(ix)
  }
  if (any(diagRinv==0 | is.na(diagRinv) | is.nan(diagRinv))) 
    boom("ERROR in the definition of Rinv")
  diagRinv<-1./diagRinv
  # if CVmode, then adjust the background
  if (argv$cv_mode) {
    xb<-array(data=NA,dim=c(ncv,nens))
    for (e in 1:nens) {
      if (argv$letkf.obsop=="obsop_LapseRateConst") {
        xb[,e]<-obsop_LapseRateConst(s=n0,
                                      m=length(aix),
                                      xm=xgrid_orig[aix],
                                      ym=ygrid_orig[aix],
                                      zm=dem_orig[aix],
                                      xs=VecX_cv,
                                      ys=VecY_cv,
                                      zs=VecZ_cv,
                                      ifield=xb_orig[,e],
                                      oval=rep(0,n0),
                                      mMinElevDiff=argv$letkf.obsop.mMinElevDiff,
                                      mMinGradient=argv$letkf.obsop.mMinGradient,
                                      mMaxGradient=argv$letkf.obsop.mMaxGradient,
                                      mSearchRadius=argv$letkf.obsop.mSearchRadius)$oval
      }
    }
  }
  # ensemble mean
  xbm<-apply(xb,MAR=1,FUN=mean)
  ybm<-apply(yb,MAR=1,FUN=mean)
  # ensemble perturbations
  Xb<-xb-xbm
  Yb<-yb-ybm
  # ensemble variances
  ybvar<-apply(Yb,MAR=1,FUN=mean)**2
  # intialize analysis vector
  xa<-xb
  xa[]<-NA
  # letkf: loop over gridpoints
  argv$letkf.Dh<-10
  argv$letkf.rloc_min<-0.0013
  argv$letkf.pmax<-200
  argv$letkf.sigma2_min<-0.1
  if (argv$verbose) t00<-Sys.time()
  for (i in 1:length(xgrid)) {
    dist<-(((xgrid[i]-VecX)**2+(ygrid[i]-VecY)**2 )**0.5)/1000.
    rloc<-exp(-0.5*(dist**2./argv$letkf.Dh**2))
    sel<-which(rloc>argv$letkf.rloc_min)
    # i-th gridpoint analysis 
    if (length(sel)>0) {
      if (length(sel)>argv$letkf.pmax) sel<-order(dist)[1:argv$letkf.pmax]
      sigma2<-max(c(mean(ybvar[sel]),argv$letkf.sigma2_min))
      Yb.i<-Yb[sel,,drop=F]
      d.i<-yo[sel]-ybm[sel]
      C.i<-t(1./sigma2*diagRinv[sel]*rloc[sel]*Yb.i)
      C1.i<-crossprod(t(C.i),Yb.i)
      Cd.i<-crossprod(t(C.i),d.i)
      rm(C.i)
      diag(C1.i)<-diag(C1.i) + (nens-1)
      Pa.i<-chol2inv(chol(C1.i))
      rm(C1.i)
      # Pa*(k-1) is the matrix for which we want to square root:
      a.eig <- eigen(Pa.i*(nens-1),symmetric=T)
      Wa <- tcrossprod(
             tcrossprod( a.eig$vectors, diag(sqrt(a.eig$values))),
             t(solve(a.eig$vectors)) )
      waa<-crossprod(Pa.i, Cd.i )
      rm(Cd.i)
      W<-Wa+as.vector(waa)
      xa[i,]<-xbm[i]+Xb[i,] %*% W
      rm(W,Wa,waa,a.eig)
    # i-th gridpoint is isolated
    } else {
      xa[i,]<-xb[i,]
    }
  } # end: letkf loop over gridpoints
  if (argv$verbose) {
    t11<-Sys.time()
    print(paste("letkf time=",round(t11-t00,1),attr(t11-t00,"unit")))
  }
  if (!argv$cv_mode & 
      ("idi" %in% argv$off_grd.variables | 
       argv$idiv_instead_of_elev)) {
    if (argv$verbose) t00<-Sys.time()
    D<-exp(-0.5*((outer(VecY,VecY,FUN="-")**2.+
                  outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
                 /argv$hyletkf.Dh_oi)**2.)
    diag(D)<-diag(D)+argv$hyletkf.eps2_oi
    InvD<-chol2inv(chol(D))
    rm(D)
    xidi<-OI_RR_fast(yo=rep(1,length(VecX)),
                     yb=rep(0,length(VecX)),
                     xb=rep(0,length(xgrid)),
                     xgrid=xgrid,
                     ygrid=ygrid,
                     VecX=VecX,
                     VecY=VecY,
                     Dh=argv$hyletkf.Dh_oi)
    if (argv$idiv_instead_of_elev) {
      G<-exp(-0.5*((outer(VecY_orig,VecY,FUN="-")**2.+
                    outer(VecX_orig,VecX,FUN="-")**2.)**0.5/1000.
                    /argv$hyletkf.Dh_oi)**2.)
      W<-tcrossprod(G,InvD)
      rm(G,InvD)
      # this is the cross-validation integral data influence ("yidiv")
      elev_for_verif<-rep(1,n0) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,n0))
      rm(W)
    }
    # InvD could be used in the following
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("idi time=",round(t11-t00,1),
                              attr(t11-t00,"unit")))
    }
  }
  if (!argv$idiv_instead_of_elev) elev_for_verif<-VecZ
  # back to the original vectors
  yo<-yo_orig
  prId<-prId_orig
  VecX<-VecX_orig
  VecY<-VecY_orig
  rm(yo_orig,prId_orig,VecX_orig,VecY_orig)
  # prepare point / gridded output
  if (!argv$cv_mode) {
    ya<-array(data=NA,dim=c(n0,length(argv$iff_fg.e)))
    yb<-array(data=NA,dim=c(n0,length(argv$iff_fg.e)))
    for (j in 1:length(argv$off_grd.variables)) {
      if (!exists("r.list")) r.list<-list()
      if (argv$off_grd.variables[j]=="analysis") {
        r<-rmaster
        r[]<-NA
        grid<-array(data=NA,dim=c(length(x),length(y),length(argv$iff_fg.e),1))
        for (e in 1:length(argv$iff_fg.e)) {
          if (argv$iff_fg.e[e] %in% valens) {
            i<-which(valens==argv$iff_fg.e[e])
            r[aix]<-xa[,i]
            grid[,,e,1]<-matrix(data=getValues(r),
                                ncol=length(y),
                                nrow=length(x))
            ya[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
            auxx<-which(ya[,e]<0 | is.na(ya[,e]) | is.nan(ya[,e]))
            if (length(auxx)>0) 
              ya[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
            rm(auxx)
          }
        }
        rm(r)
        r.list[[j]]<-grid
        rm(grid)
      } else if (argv$off_grd.variables[j]=="background") {
        r<-rmaster
        r[]<-NA
        grid<-array(data=NA,dim=c(length(x),length(y),length(argv$iff_fg.e),1))
        for (e in 1:length(argv$iff_fg.e)) {
          if (argv$iff_fg.e[e] %in% valens) {
            i<-which(valens==argv$iff_fg.e[e])
            r[aix]<-xb[,i]
            grid[,,e,1]<-matrix(data=getValues(r),
                                ncol=length(y),
                                nrow=length(x))
            yb[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
            auxx<-which(yb[,e]<0 | is.na(yb[,e]) | is.nan(yb[,e]))
            if (length(auxx)>0) 
              yb[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
            rm(auxx)
          }
        }
        rm(r)
        r.list[[j]]<-grid
        rm(grid)
      } else if (argv$off_grd.variables[j]=="idi") {
        r<-rmaster
        r[]<-NA
        r[aix]<-xidi
        r.list[[j]]<-matrix(data=getValues(r),
                            ncol=length(y),
                            nrow=length(x))
        rm(r)
      }
    }
  } # end prepare for gridded output
# END of Local Ensemble Transform Kalman Filter
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
  if (argv$verbose) t00<-Sys.time() 
  #CVmode for deterministic analysis
  if (is.null(argv$iff_fg.epos)) {
    ya_cv<-extract(ra,cbind(VecX_cv,VecY_cv),method="bilinear")
    Dh<-argv$Dh
  #CVmode for ensemble analysis
  } else {
    ya_cv<-array(data=NA,dim=c(ncv,length(argv$iff_fg.e)))
    yb_cv<-array(data=NA,dim=c(ncv,length(argv$iff_fg.e)))
    for (e in 1:length(argv$iff_fg.e)) {
      if (argv$iff_fg.e[e] %in% valens) {
        i<-which(valens==argv$iff_fg.e[e])
        ya_cv[,e]<-xa[,i]
        yb_cv[,e]<-xb[,i]
      }
    }
    if (argv$mode=="hyletkf") Dh<-argv$hyletkf.Dh_oi
  }
  if (argv$idiv_instead_of_elev) {
    D<-exp(-0.5*((outer(VecY,VecY,FUN="-")**2.+
                  outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
                 /Dh)**2.)
    diag(D)<-diag(D)+argv$eps2
    InvD<-chol2inv(chol(D))
    rm(D)
    G<-exp(-0.5*((outer(VecY_cv,VecY,FUN="-")**2.+
                  outer(VecX_cv,VecX,FUN="-")**2.)**0.5/1000.
                 /Dh)**2.)
    W<-tcrossprod(G,InvD)
    rm(G,InvD)
    # this is the cross-validation integral data influence ("yidiv")
    elev_for_verif_cv<-rep(1,ncv) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,ncv))
    rm(W)
  } else {
    elev_for_verif_cv<-VecZ_cv
  }
  if (argv$verbose) {
    t11<-Sys.time()
    print(paste("cv-mode prepare output time=",round(t11-t00,1),
                                               attr(t11-t00,"unit")))
  }
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (argv$verbose) print("++ Output")
#------------------------------------------------------------------------------
# Station Points - CVmode
if (argv$cv_mode) {
  # Station Points - CVmode - case of deterministic analysis
  if (is.null(argv$iff_fg.epos)) { 
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
  # Station Points - CVmode - case of case of ensemble analysis
  } else {
    # off_cvvernc_a / off_cvvernc_b
    # write off_cvvernc_a (cv-analysis)
    tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                      format="%Y%m%d%H%M",tz="GMT")
    dim_t<-list(name="time",
                units="H",
                vals=tstamp_nc,
                is_time=T,
                unlim=T,
                times.unit="S",
                times.format="%Y%m%d%H%M",
                timeref="197001010000",
                timeref.format="%Y%m%d%H%M")
    dim_lt<-list(name="leadtime",units="",vals=0)
    dim_loc<-list(name="location",units="",vals=VecS_cv)
    dim_ens<-list(name="ensemble_member",units="",vals=argv$iff_fg.e)
    dim_list<-list(dim_t,dim_lt,dim_loc,dim_ens)
    nr<-length(argv$off_vernc.thresholds)
    nq<-length(argv$off_vernc.quantile)
    nlt<-1
    nt<-1
    nloc<-ncv
    dim_r<-list(name="threshold",units="",vals=argv$off_vernc.thresholds)
    dim_list[[(length(dim_list)+1)]]<-dim_r
    dim_q<-list(name="quantile",units="",vals=argv$off_vernc.quantile)
    dim_list[[(length(dim_list)+1)]]<-dim_q
    # variables
    var_lat<-list(name="lat",
                  units="",
                  dim_names=c("location"),
                  vals=array(VecLat_cv,dim=c(nloc)))
    var_lon<-list(name="lon",
                  units="",
                  dim_names=c("location"),
                  vals=array(VecLon_cv,dim=c(nloc)))
    var_ele<-list(name="altitude",
                  units="",
                  dim_names=c("location"),
                  vals=array(elev_for_verif_cv,dim=c(nloc)))
    obs<-array(data=NA,dim=c(nloc,nlt,nt))
    obs[,nlt,nt]<-yo_cv
    var_obs<-list(name="obs",
                  units="",
                  dim_names=c("location","leadtime","time"),
                  vals=obs)
    rm(obs)
    fcst<-array(data=NA,dim=c(nloc,nlt,nt))
    fcst[,nlt,nt]<-rowMeans(ya_cv,na.rm=T)
    var_fcst<-list(name="fcst",
                   units="",
                   dim_names=c("location","leadtime","time"),
                   vals=fcst)
    rm(fcst)
    xval<-array(data=NA,dim=c(length(argv$iff_fg.e),nloc,nlt,nt))
    xval[,,nlt,nt]<-ya_cv
    var_ens<-list(name="ensemble",
                  units="",
                  dim_names=c("ensemble_member","location","leadtime","time"),
                  vals=xval)
    rm(xval)
    var_list<-list(var_lat,var_lon,var_ele,var_obs,var_fcst,var_ens)
    pit<-array(data=NA,dim=c(nloc,nlt,nt))
    pit[,nlt,nt]<-apply(cbind(ya_cv,yo_cv),
                     MARGIN=1,
                     FUN=function(x){
                       if (is.na(x[length(x)])) return(NA)
                       ix<-which(!is.na(x[1:(length(x)-1)]))
                       if (length(ix)==0) return(NA)
                       ecdf(x[ix])(x[length(x)])})
    var_pit<-list(name="pit",
                  units="",
                  dim_names=c("location","leadtime","time"),
                  vals=pit)
    rm(pit)
    var_list[[(length(var_list)+1)]]<-var_pit
    rm(var_pit)
    cdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    pdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    cdf[1:nr,1:nloc,nlt,nt]<-apply(ya_cv, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        ecdf(x[ix])(argv$off_vernc.thresholds)})
    pdf[1:nr,1:nloc,nlt,nt]<-apply(ya_cv, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        approxfun(density(x[ix]),yleft=0,yright=0)(argv$off_vernc.thresholds)})
    var_cdf<-list(name="cdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=cdf)
    rm(cdf)
    var_list[[(length(var_list)+1)]]<-var_cdf
    var_pdf<-list(name="pdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=pdf)
    rm(pdf)
    var_list[[(length(var_list)+1)]]<-var_pdf
    qx<-array(data=NA,dim=c(nq,nloc,nlt,nt))
    qx[1:nq,1:nloc,nlt,nt]<-apply(ya_cv,MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        as.vector(quantile(x[ix],probs=argv$off_vernc.quantile))})
    var_qx<-list(name="x",
                 units="",
                 dim_names=c("quantile","location","leadtime","time"),
                 vals=qx)
    rm(qx)
    var_list[[(length(var_list)+1)]]<-var_qx
    gatt_longn<-list(attname="long_name",attval=argv$off_vernc.varname,prec="text")
    gatt_stdn<-list(attname="standard_name",attval=argv$off_vernc.stdvarname,prec="text")
    gatt_units<-list(attname="units",attval=argv$off_vernc.varunits,prec="text")
    res<-write_generic_dotnc(nc.file=argv$off_cvvernc_a,
                             dims=dim_list,
                             vars=var_list,
                             glob_attrs=list(gatt_longn,gatt_stdn,gatt_units))
    if (is.null(res)) {
      print(paste("ERROR while writing ",argv$off_cvvernc_a))
    }
    print(paste("written file",argv$off_cvvernc_a))
    # write off_cvvernc_b (background)
    fcst<-array(data=NA,dim=c(nloc,nlt,nt))
    fcst[,nlt,nt]<-rowMeans(yb_cv,na.rm=T)
    var_fcst<-list(name="fcst",
                   units="",
                   dim_names=c("location","leadtime","time"),
                   vals=fcst)
    rm(fcst)
    xval<-array(data=NA,dim=c(length(argv$iff_fg.e),nloc,nlt,nt))
    xval[,,nlt,nt]<-yb_cv
    var_ens<-list(name="ensemble",
                  units="",
                  dim_names=c("ensemble_member","location","leadtime","time"),
                  vals=xval)
    rm(xval)
    var_list<-list(var_lat,var_lon,var_ele,var_obs,var_fcst,var_ens)
    pit<-array(data=NA,dim=c(nloc,nlt,nt))
    pit[,nlt,nt]<-apply(cbind(yb_cv,yo_cv),
                     MARGIN=1,
                     FUN=function(x){
                       if (is.na(x[length(x)])) return(NA)
                       ix<-which(!is.na(x[1:(length(x)-1)]))
                       if (length(ix)==0) return(NA)
                       ecdf(x[ix])(x[length(x)])})
    var_pit<-list(name="pit",
                  units="",
                  dim_names=c("location","leadtime","time"),
                  vals=pit)
    rm(pit)
    var_list[[(length(var_list)+1)]]<-var_pit
    rm(var_pit)
    cdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    pdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    cdf[1:nr,1:nloc,nlt,nt]<-apply(yb_cv, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        ecdf(x[ix])(argv$off_vernc.thresholds)})
    pdf[1:nr,1:nloc,nlt,nt]<-apply(yb_cv, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        approxfun(density(x[ix]),yleft=0,yright=0)(argv$off_vernc.thresholds)})
    var_cdf<-list(name="cdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=cdf)
    rm(cdf)
    var_list[[(length(var_list)+1)]]<-var_cdf
    var_pdf<-list(name="pdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=pdf)
    rm(pdf)
    var_list[[(length(var_list)+1)]]<-var_pdf
    qx<-array(data=NA,dim=c(nq,nloc,nlt,nt))
    qx[1:nq,1:nloc,nlt,nt]<-apply(yb_cv,MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        as.vector(quantile(x[ix],probs=argv$off_vernc.quantile))})
    var_qx<-list(name="x",
                 units="",
                 dim_names=c("quantile","location","leadtime","time"),
                 vals=qx)
    rm(qx)
    var_list[[(length(var_list)+1)]]<-var_qx
    res<-write_generic_dotnc(nc.file=argv$off_cvvernc_b,
                             dims=dim_list,
                             vars=var_list,
                             glob_attrs=list(gatt_longn,gatt_stdn,gatt_units))
    if (is.null(res)) {
      print(paste("ERROR while writing ",argv$off_cvvernc_b))
    }
    print(paste("written file",argv$off_cvvernc_b))
  }
#------------------------------------------------------------------------------
# non-CVmode
} else {
  # Station Points - non-CVmode - case of deterministic analysis
  if (is.null(argv$iff_fg.epos)) { 
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
  # Station Points - non-CVmode - case of ensemble analysis
  } else {
    # off_vernc_a / off_vernc_b
    # write off_vernc_a (analysis)
    tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                      format="%Y%m%d%H%M",tz="GMT")
    dim_t<-list(name="time",
                units="H",
                vals=tstamp_nc,
                is_time=T,
                unlim=T,
                times.unit="S",
                times.format="%Y%m%d%H%M",
                timeref="197001010000",
                timeref.format="%Y%m%d%H%M")
    dim_lt<-list(name="leadtime",units="",vals=0)
    dim_loc<-list(name="location",units="",vals=VecS)
    dim_ens<-list(name="ensemble_member",units="",vals=argv$iff_fg.e)
    dim_list<-list(dim_t,dim_lt,dim_loc,dim_ens)
    nr<-length(argv$off_vernc.thresholds)
    nq<-length(argv$off_vernc.quantile)
    nlt<-1
    nt<-1
    nloc<-n0
    dim_r<-list(name="threshold",units="",vals=argv$off_vernc.thresholds)
    dim_list[[(length(dim_list)+1)]]<-dim_r
    dim_q<-list(name="quantile",units="",vals=argv$off_vernc.quantile)
    dim_list[[(length(dim_list)+1)]]<-dim_q
    # variables
    var_lat<-list(name="lat",
                  units="",
                  dim_names=c("location"),
                  vals=array(VecLat,dim=c(nloc)))
    var_lon<-list(name="lon",
                  units="",
                  dim_names=c("location"),
                  vals=array(VecLon,dim=c(nloc)))
    var_ele<-list(name="altitude",
                  units="",
                  dim_names=c("location"),
                  vals=array(elev_for_verif,dim=c(nloc)))
    obs<-array(data=NA,dim=c(nloc,nlt,nt))
    obs[,nlt,nt]<-yo
    var_obs<-list(name="obs",
                  units="",
                  dim_names=c("location","leadtime","time"),
                  vals=obs)
    rm(obs)
    fcst<-array(data=NA,dim=c(nloc,nlt,nt))
    fcst[,nlt,nt]<-rowMeans(ya,na.rm=T)
    var_fcst<-list(name="fcst",
                   units="",
                   dim_names=c("location","leadtime","time"),
                   vals=fcst)
    rm(fcst)
    xval<-array(data=NA,dim=c(length(argv$iff_fg.e),nloc,nlt,nt))
    xval[,,nlt,nt]<-ya
    var_ens<-list(name="ensemble",
                  units="",
                  dim_names=c("ensemble_member","location","leadtime","time"),
                  vals=xval)
    rm(xval)
    var_list<-list(var_lat,var_lon,var_ele,var_obs,var_fcst,var_ens)
    pit<-array(data=NA,dim=c(nloc,nlt,nt))
    pit[,nlt,nt]<-apply(cbind(ya,yo),
                     MARGIN=1,
                     FUN=function(x){
                       if (is.na(x[length(x)])) return(NA)
                       ix<-which(!is.na(x[1:(length(x)-1)]))
                       if (length(ix)==0) return(NA)
                       ecdf(x[ix])(x[length(x)])})
    var_pit<-list(name="pit",
                  units="",
                  dim_names=c("location","leadtime","time"),
                  vals=pit)
    rm(pit)
    var_list[[(length(var_list)+1)]]<-var_pit
    rm(var_pit)
    cdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    pdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    cdf[1:nr,1:nloc,nlt,nt]<-apply(ya, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        ecdf(x[ix])(argv$off_vernc.thresholds)})
    pdf[1:nr,1:nloc,nlt,nt]<-apply(ya, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        approxfun(density(x[ix]),yleft=0,yright=0)(argv$off_vernc.thresholds)})
    var_cdf<-list(name="cdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=cdf)
    rm(cdf)
    var_list[[(length(var_list)+1)]]<-var_cdf
    var_pdf<-list(name="pdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=pdf)
    rm(pdf)
    var_list[[(length(var_list)+1)]]<-var_pdf
    qx<-array(data=NA,dim=c(nq,nloc,nlt,nt))
    qx[1:nq,1:nloc,nlt,nt]<-apply(ya,MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        as.vector(quantile(x[ix],probs=argv$off_vernc.quantile))})
    var_qx<-list(name="x",
                 units="",
                 dim_names=c("quantile","location","leadtime","time"),
                 vals=qx)
    rm(qx)
    var_list[[(length(var_list)+1)]]<-var_qx
    gatt_longn<-list(attname="long_name",attval=argv$off_vernc.varname,prec="text")
    gatt_stdn<-list(attname="standard_name",attval=argv$off_vernc.stdvarname,prec="text")
    gatt_units<-list(attname="units",attval=argv$off_vernc.varunits,prec="text")
    res<-write_generic_dotnc(nc.file=argv$off_vernc_a,
                             dims=dim_list,
                             vars=var_list,
                             glob_attrs=list(gatt_longn,gatt_stdn,gatt_units))
    if (is.null(res)) {
      print(paste("ERROR while writing ",argv$off_vernc_a))
    }
    print(paste("written file",argv$off_vernc_a))
    # write off_vernc_b (background)
    fcst<-array(data=NA,dim=c(nloc,nlt,nt))
    fcst[,nlt,nt]<-rowMeans(yb,na.rm=T)
    var_fcst<-list(name="fcst",
                   units="",
                   dim_names=c("location","leadtime","time"),
                   vals=fcst)
    rm(fcst)
    xval<-array(data=NA,dim=c(length(argv$iff_fg.e),nloc,nlt,nt))
    xval[,,nlt,nt]<-yb
    var_ens<-list(name="ensemble",
                  units="",
                  dim_names=c("ensemble_member","location","leadtime","time"),
                  vals=xval)
    rm(xval)
    var_list<-list(var_lat,var_lon,var_ele,var_obs,var_fcst,var_ens)
    pit<-array(data=NA,dim=c(nloc,nlt,nt))
    pit[,nlt,nt]<-apply(cbind(yb,yo),
                     MARGIN=1,
                     FUN=function(x){
                       if (is.na(x[length(x)])) return(NA)
                       ix<-which(!is.na(x[1:(length(x)-1)]))
                       if (length(ix)==0) return(NA)
                       ecdf(x[ix])(x[length(x)])})
    var_pit<-list(name="pit",
                  units="",
                  dim_names=c("location","leadtime","time"),
                  vals=pit)
    rm(pit)
    var_list[[(length(var_list)+1)]]<-var_pit
    rm(var_pit)
    cdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    pdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    cdf[1:nr,1:nloc,nlt,nt]<-apply(yb, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        ecdf(x[ix])(argv$off_vernc.thresholds)})
    pdf[1:nr,1:nloc,nlt,nt]<-apply(yb, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        approxfun(density(x[ix]),yleft=0,yright=0)(argv$off_vernc.thresholds)})
    var_cdf<-list(name="cdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=cdf)
    rm(cdf)
    var_list[[(length(var_list)+1)]]<-var_cdf
    var_pdf<-list(name="pdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=pdf)
    rm(pdf)
    var_list[[(length(var_list)+1)]]<-var_pdf
    qx<-array(data=NA,dim=c(nq,nloc,nlt,nt))
    qx[1:nq,1:nloc,nlt,nt]<-apply(yb,MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        as.vector(quantile(x[ix],probs=argv$off_vernc.quantile))})
    var_qx<-list(name="x",
                 units="",
                 dim_names=c("quantile","location","leadtime","time"),
                 vals=qx)
    rm(qx)
    var_list[[(length(var_list)+1)]]<-var_qx
    res<-write_generic_dotnc(nc.file=argv$off_vernc_b,
                             dims=dim_list,
                             vars=var_list,
                             glob_attrs=list(gatt_longn,gatt_stdn,gatt_units))
    if (is.null(res)) {
      print(paste("ERROR while writing ",argv$off_vernc_b))
    }
    print(paste("written file",argv$off_vernc_b))
  }
#------------------------------------------------------------------------------
  # Gridpoints - non-CVmode(only)
  # NOTE: deterministic and ensemble analysis share the same output session
  # to customize the output, use the config argument "off_grd.variables"
  # default output: just the analysis
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
