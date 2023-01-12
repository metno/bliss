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
    if (!argv$twostep_nogrid) {
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
    }
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

#------------------------------------------------------------------------------
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

