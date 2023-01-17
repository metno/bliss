#+
oi_twostep_senorge_temperature <- function( argv, y_env, env) {
# Spatial scale definitions. 
#  Regional=whole domain
#  Sub-regional (or local scale)=dozens of observations (10-100)
#  small-scale=few observations (1-10)
#  sub-grid scale=not observed
#------------------------------------------------------------------------------
  t0a <- Sys.time()

  cat( "-- OI two-step developed for seNorge --\n")

  y_env$super_yo  <- list()
  y_env$centroids <- list()

  if (argv$twostep_superobbing) {
    # not yet implemented
    y_env$super_yo$x     <- y_env$yo$x
    y_env$super_yo$y     <- y_env$yo$y
    y_env$super_yo$z     <- y_env$yo$z
    y_env$super_yo$value <- y_env$yo$value
    y_env$super_yo$laf   <- y_env$yo$laf
  } else {
    y_env$super_yo$x     <- y_env$yo$x
    y_env$super_yo$y     <- y_env$yo$y
    y_env$super_yo$z     <- y_env$yo$z
    y_env$super_yo$value <- y_env$yo$value
    y_env$super_yo$laf   <- y_env$yo$laf
  }
  y_env$super_yo$n <- length(y_env$yo$value)

  # Regional backgrounds
  cat("+----------------------------+\n")
  cat("Regional backgroundg field\n")
  # Regional background field is the average of local vertical temperature profiles, 
  # each of them computed around a so-called centroid (a point at the center of a
  # sub-region)

  # define grid used to compute local backgrounds 
  # (grid nodes are centroid candidates)
  y_env$centroids$r <- raster( extent(env$rmaster),
                               ncol=argv$oi2step.bg_centroids_nrnc[2],
                               nrow=argv$oi2step.bg_centroids_nrnc[1],
                               crs=crs(env$rmaster))
  y_env$centroids$res_mean <- round( mean( res(y_env$centroids$r)))
  xy    <- xyFromCell( y_env$centroids$r, 1:ncell(y_env$centroids$r))
  # xr and yr are the candidate centroids
  xr    <- xy[,1]
  yr    <- xy[,2]
  rm(xy)
  y_env$centroids$r[] <- 1:ncell(y_env$centroids$r)
  r_notmasked <- extract( env$rmaster, cbind( xr, yr), 
                          buffer=argv$oi2step.bg_centroids_buffer, 
                          na.rm=T, fun=mean)
  # selection of centroids
  #  centroids are nodes of r where -within a predefined distance (oi2step.bg_obsbufferlength)-
  #  these two conditions hold true
  #  (a) there is at least one grid point of rmaster that is not masked
  #  (b) there are at least oi2step.bg_obsnmin4centroid observations
#  irx <- which( ( (1:ncell(r)) %in% as.vector( na.omit( unique( 
#                   extract( r, cbind(y_env$super_yo$x,y_env$super_yo$y)))))) & 
#                ( !is.na( extract( env$rmaster, cbind(xr,yr), 
#                          buffer=argv$oi2step.bg_obsbufferlength, 
#                          na.rm=T, fun=mean))))
#  irx <- which( ( !is.na( extract( env$rmaster, cbind(xr,yr), 
#                          buffer=argv$oi2step.bg_obsbufferlength, 
#                          na.rm=T, fun=mean))))
  nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
                     query = cbind( xr, yr), 
                     k = argv$oi2step.bg_centroids_nobsmin, 
                     searchtype = "radius", 
                     radius = argv$oi2step.bg_centroids_buffer)
  y_env$centroids$i <- integer(0)
  for (i in 1:length(xr)) {
    if ( (length( which(nn2$nn.idx[i,]!=0)) < argv$oi2step.bg_centroids_nobsmin) | 
         is.na(r_notmasked[i])) next
    y_env$centroids$i <- c(y_env$centroids$i, i)
  } 
  y_env$centroids$n <- length(y_env$centroids$i)

  cat(paste("the master grid has been divided in",
              argv$oi2step.bg_centroids_nrnc[1],"x",argv$oi2step.bg_centroids_nrnc[2],"boxes\n"))
  cat(paste("sub-regional area extensions (length x (m),length y (m))=",
        round(res(y_env$centroids$r)[1]),round(res(y_env$centroids$r)[2])),"\n")
  cat(paste("reference (horizontal) length scale to weight the sub-regional backgrounds (m)=",round(mean(res(y_env$centroids$r))),"\n"))
  cat(paste("# sub-regional centroids", y_env$centroids$n, "\n"))

  # save results in a global variable
  y_env$centroids$x <- xr[y_env$centroids$i]
  y_env$centroids$y <- yr[y_env$centroids$i]

  # clean memory
  rm(nn2, r_notmasked, xr, yr)

  #
  envtmp$nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
                            query = cbind( y_env$centroids$x, y_env$centroids$y), 
                            k = y_env$super_yo$n, 
                            searchtype = "radius", 
                            radius = argv$oi2step.bg_centroids_buffer)

  # run OI gridpoint by gridpoint
  if (!is.na(argv$cores)) {
    res <- t( mcmapply( vertical_profile_at_centroid_senorge2018,
                        1:y_env$centroids$n,
                        mc.cores=argv$cores,
                        SIMPLIFY=T))
  # no-multicores
  } else {
    res <- t( mapply( vertical_profile_at_centroid_senorge2018,
                      1:y_env$centroids$n,
                      SIMPLIFY=T))
  }
save(file="tmp.rdata",envtmp,y_env,res,env,argv)
q()
  # count the number of observations in each box
  rnobs <- rasterize( cbind( y_env$super_yo$x, y_env$super_yo$y), y_env$centroids$r, 
                      y_env$super_yo$value,fun=function(x,...)length(x))
  nr <- getValues(rnobs)
  # create the 4D array for the function call ...
  # ... via apply centroids=(xr[irx],yr[irx])
#  ixyn <- cbind( (1:ncell(r))[irx], xr[irx], yr[irx], nr[irx])
save(file="tmp.rdata",y_env,env,nn2,xr,yr,r_notmasked)
  # Bye-bye
  t1a <- Sys.time()
  cat( paste( "OI two-step developed for seNorge, total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))
q()
  # Regional backgrounds
  cat("+----------------------------+")
  print("two-step interpolation: regional bg field")
  # Spatial scale definitions. 
  #  Regional=whole domain
  #  Sub-regional (or local scale)=dozens of observations (10-100)
  #  small-scale=few observations (1-10)
  #  sub-grid scale=not observed
  # define grid used to compute local backgrounds 
  #  (grid nodes are candidates for the centroids)
  r<-raster( extent(xmn,xmx,ymn,ymx),
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
    if (!(argv$cv_mode|argv$cv_mode_random)) {
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
    if (argv$twostep_nogrid) {
      if (!any(yb==na)) {
        b_ok<-T
        print(paste("maxboxl=",(i*argv$maxboxl),"m"))
        break
      }
    } else {
      if (!(any(yb==na) | ((!(argv$cv_mode|argv$cv_mode_random)) & (any(xb==na))))) {
        b_ok<-T
        print(paste("maxboxl=",(i*argv$maxboxl),"m"))
        break
      }
    }
  }
  if (!b_ok) {
    if (any(yb==na)) 
      boom("ERROR: found NAs in background values at station locations (try to increase \"maxboxl\"")
    if (!(argv$cv_mode|argv$cv_mode_random)) {
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
    if (!(argv$cv_mode|argv$cv_mode_random)) {
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
  if (!(argv$cv_mode|argv$cv_mode_random)) {
    if (!argv$twostep_nogrid) {
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
      for (i in 1:length(argv$off_x.variables)) {
        if (!exists("r.list")) r.list<-list()
        if (argv$off_x.variables[i]=="analysis") {
          ra<-rmaster
          ra[]<-NA
          ra[mask]<-xout[1,]
          r.list[[i]]<-matrix(data=getValues(ra),
                              ncol=length(y),
                              nrow=length(x))
          xa  <- getValues(ra)
          ng  <- length(getValues(ra))
          aix <- 1:length(getValues(ra))
        } else if (argv$off_x.variables[i]=="background") {
          r<-rmaster
          r[]<-NA
          r[mask]<-xb
          r.list[[i]]<-matrix(data=getValues(r),
                              ncol=length(y),
                              nrow=length(x))
          xb  <- getValues(r)
          ng  <- length(getValues(r))
          aix <- 1:length(getValues(r))
          rm(r)
        } else if (argv$off_x.variables[i]=="idi") {
          r<-rmaster
          r[]<-NA
          r[mask]<-xout[2,]
          r.list[[i]]<-matrix(data=getValues(r),
                              ncol=length(y),
                              nrow=length(x))
          xidi <- getValues(r)
          ng  <- length(getValues(r))
          aix  <- 1:length(getValues(r))
          rm(r)
        } else if (argv$off_x.variables[i]=="dh") {
          r<-rmaster
          r[]<-NA
          r[mask]<-xdh_oi
          r.list[[i]]<-matrix(data=getValues(r),
                              ncol=length(y),
                              nrow=length(x))
          xdh <- getValues(r)
          ng  <- length(getValues(r))
          aix <- 1:length(getValues(r))
          rm(r)
        }
      }
      rm(xout)
    }
    if (!exists("t000")) t000<-Sys.time()
    # station points
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
}
