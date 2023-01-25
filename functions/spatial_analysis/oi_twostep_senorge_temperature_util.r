#+
vertical_profile_at_centroid_senorge2018 <- function( i) {
#print(i)
  if ( length( ix <- which( envtmp$nn2$nn.idx[i,] != 0)) == 0) 
    return(c( NA, NA, NA, NA, NA, NA))
 
  z_opt <- y_env$super_yo$z[envtmp$nn2$nn.idx[i,ix]]
  v_opt <- y_env$super_yo$value[envtmp$nn2$nn.idx[i,ix]]

  dz <- as.numeric( diff( quantile( z_opt, probs=c(0.05,0.95))))

  if (dz < argv$oi2step.bg_vertprof_dzmin) {
    par  <- mean( v_opt, na.rm=T)
    opt  <- optimize( f=vert_prof_basic_opt, 
                      interval=c(argv$oi2step.bg_vertprof_vmin, argv$oi2step.bg_vertprof_vmax),
                      z=z_opt, v=v_opt, gamma=argv$oi2step.bg_vertprof_gamma)
    tpar <- c( 1, opt$minimum, argv$oi2step.bg_vertprof_gamma, NA, NA, NA)
  } else {
    ui <- matrix( ncol=5, nrow=10, data=0)
    ci <- vector( length=10)
    q  <- as.numeric( quantile( z_opt, probs=c( 0.1, 0.25, 0.75)))
    v_range <- range(v_opt,na.rm=T)
    for (i in 1:5) {
      ui[(2*i-1),i] <- 1
      ui[(2*i),i]   <- (-1)
    }
    ci[1:2] <- c( v_range[1], -v_range[2])
    aux <- sort( c( 0.5*argv$oi2step.bg_vertprof_gamma, 2*argv$oi2step.bg_vertprof_gamma))
    ci[3:4] <- c( aux[1], -aux[2])
    ci[5:6] <- c(      0, -20)
    if ( q[1] == q[2]) {
      aux <- q[2] + mean(z_opt,na.rm=T) / 10
    } else {
      aux <- q[2]
    }
    ci[7:8] <- c( -10, -aux)
    aux <- q[3] - q[2] 
    ci[9:10] <- c( 0, -aux)
    par <- c( mean( v_range),
              argv$oi2step.bg_vertprof_gamma,
              5,
              q[1],
             -ci[10] / 2)
    opt <- constrOptim( theta=par,
                        f=vert_prof_Frei_opt,
                        ui=ui,
                        ci=ci,
                        grad=NULL,
                        z=z_opt, v=v_opt)
    tpar <- c( 2, opt$par[1], opt$par[2], opt$par[3], opt$par[4], opt$par[5])
  }  
#print(tpar)
  return(tpar)

}

#+
blend_vertical_profiles_senorge2018 <- function( i, 
                                                 corr = "soar",
                                                 dh = 10000) {
  if (i==1) cat(".")

  rloc <- corr1d( envtmp$nn2$nn.dist[i,], dh, corr)
  values <- vector( mode="numeric", length=y_env$centroids$n)
  if ( envtmp$n1 > 0) 
    values[envtmp$ix1] <- vert_prof_basic( rep( envtmp$z[i], envtmp$n1), 
                           y_env$centroids$vert_prof[envtmp$ix1,2], 
                           y_env$centroids$vert_prof[envtmp$ix1,3]) 
  if ( envtmp$n2 > 0) 
    values[envtmp$ix2] <- vert_prof_Frei_1( rep( envtmp$z[i], envtmp$n2), 
                           y_env$centroids$vert_prof[envtmp$ix2,2], 
                           y_env$centroids$vert_prof[envtmp$ix2,3],
                           y_env$centroids$vert_prof[envtmp$ix2,4],
                           y_env$centroids$vert_prof[envtmp$ix2,5],
                           y_env$centroids$vert_prof[envtmp$ix2,6])

  sum(rloc * values) / sum(rloc)

}

#+ vertical profile of temperature (linear)
vert_prof_basic <- function( z, t0, gamma) {
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  return( t0 + gamma * z)
}


#+ vertical profile of temperature (Frei, 2014)
vert_prof_Frei<-function( z, t0, gamma, a, h0, h1i) {
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
  t   <- z
  t[] <- NA
  h1  <- h0 + h1i
  if ( length( (z.le.h0 <- which( z <= h0))) > 0) 
    t[z.le.h0] <- t0 + gamma * z[z.le.h0] - a 
  if ( length( (z.ge.h1 <- which( z >= h1))) > 0) 
    t[z.ge.h1] <- t0 + gamma * z[z.ge.h1]
  if ( length( (z.in <- which( z > h0 & z < h1)))>0) 
    t[z.in] <- t0 + gamma * z[z.in] - a/2 * (1 + cos(pi*(z[z.in]-h0)/h1i))
  return(t)
}

#+ vertical profile of temperature (Frei, 2014)
vert_prof_Frei_1<-function( z, t0, gamma, a, h0, h1i) {
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
  t   <- z
  t[] <- NA
  h1  <- h0 + h1i
  if ( length( (z.le.h0 <- which( z <= h0))) > 0) 
    t[z.le.h0] <- t0[z.le.h0] + gamma[z.le.h0] * z[z.le.h0] - a[z.le.h0] 
  if ( length( (z.ge.h1 <- which( z >= h1))) > 0) 
    t[z.ge.h1] <- t0[z.ge.h1] + gamma[z.ge.h1] * z[z.ge.h1]
  if ( length( (z.in <- which( z > h0 & z < h1)))>0) 
    t[z.in] <- t0[z.in] + gamma[z.in] * z[z.in] - a[z.in]/2 * (1 + cos(pi*(z[z.in]-h0[z.in])/h1i[z.in]))
  return(t)
}

#+ cost function used for optimization of tvertprof parameter
vert_prof_Frei_opt<-function(par, z, v) {
  return( log( sqrt( mean( (vert_prof_Frei(z=z,t0=par[1],gamma=par[2],a=par[3],h0=par[4],h1i=par[5]) - v)**2))))
}


#+ cost function used for optimization of tvertprof parameter
vert_prof_basic_opt <- function( par, z, v, gamma) {
#  te <- vert_prof_basic( z=zopt, t0=par[1], gamma=gamma)
  return( log( sqrt( mean( (vert_prof_basic( z=z, t0=par, gamma=gamma) - v)**2))))
}

#+
Theil_Sen_regression <-function( x, y) {
# y = res[1] + res[2] * x
#------------------------------------------------------------------------------
# REF: Wilks (2019) p. 283
  bnum <- outer( y, y, FUN="-")
  bden <- outer( x, x, FUN="-")
  bset <- bnum / bden
  b <- median( bset[row(bset)>col(bset) & is.finite(bset)], na.rm=T)
  if ( !(!is.finite(b) | is.na(b) | is.null(b))) {
    residuals <- y - b * x
    a   <- median( residuals)
    res <- c( a, b)
  } else {
    res <- c( NA, NA)
  }
  res
}

#+
regional_background<-function( argv, y_env) { # sub-region centroids
#                           maxboxl=250000,
#                           nmin=50,
#                           refdist=50000,
#                           dzmin=30,
#                           eps2=0.25,
#                           closeNth=4
#                           ) {
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

