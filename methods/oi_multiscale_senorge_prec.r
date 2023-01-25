#+
oi_multiscale_senorge_prec <- function( argv, y_env, fg_env, env) {
#------------------------------------------------------------------------------
  t0a <- Sys.time()

  cat( "-- OI multiscale developed for seNorge --\n")
  # If there is at least one observation of precipitation (not dry everywhere)
  if (y_env$yo$nwet > 0) {

    # We need to have just one background field
    if ( ( nfg <- length( fg_env$fg)) != 1) return( FALSE)
    rrf <- subset( fg_env$fg[[1]]$r_main, subset=1)
    yrf <- extract( rrf, cbind(y_env$yo$x, y_env$yo$y))

    # Master grid average grid spacing
    res_master <- mean( res(env$rmaster))

    # Distances between all stations
    envtmp$dist2 <- outer(y_env$yo$x,y_env$yo$x,FUN="-")**2.+
                    outer(y_env$yo$y,y_env$yo$y,FUN="-")**2.

    # IDI
    if ( !is.na(argv$dh_idi)) {
      cat("IDI")
      D <- corr2d( par=argv$dh_idi, label=argv$corrfun, values_as_globalVar=T)
      diag(D) <- diag(D) + 0.1
      envtmp$SRinv_di <- crossprod( chol2inv(chol(D)), array(data=rep(1,y_env$yo$n),dim=c(y_env$yo$n,1)))
      envtmp$Eb  <- array( data=0, dim=c( env$ngrid, 1))
      envtmp$x <- env$xgrid[env$mask]
      envtmp$y <- env$ygrid[env$mask]
      envtmp$m_dim <- env$ngrid
      envtmp$obs_x <- y_env$yo$x
      envtmp$obs_y <- y_env$yo$y
      # run OI gridpoint by gridpoint
      if (!is.na(argv$cores)) {
        res <- t( mcmapply( enoi_basicFaster_gridpoint_by_gridpoint,
                            1:envtmp$m_dim,
                            mc.cores=argv$cores,
                            SIMPLIFY=T,
                            MoreArgs = list( corr=argv$corrfun,
                                             dh=argv$dh_idi,
                                             idi=F)))
      # no-multicores
      } else {
        res <- t( mapply( enoi_basicFaster_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corrfun, 
                                           dh=argv$dh_idi,
                                           idi=F)))
      }
      xidi <- res[,1]
      rm(res)
      cat("OK\n")
    } else {
      xidi <- rep( 0, env$ngrid)
    }

    # Define sequence of spatial scales (coarser to finer)
    if (y_env$yo$n < 100) {
      kseq <- rev( unique( c(2,3,4,seq(5,y_env$yo$n,by=5), y_env$yo$n)))
    } else {
      kseq <- rev( unique( c(2,3,4,seq(5,100,by=5),seq(100,y_env$yo$n,by=200), y_env$yo$n)))
    }
    # distances between each station and all the others (increasing order) 
    nn2 <- nn2( cbind( y_env$yo$x,y_env$yo$y), 
                       query = cbind(y_env$yo$x,y_env$yo$y), 
                       k = y_env$yo$n, searchtype = "radius", 
                       radius = 100000000)
    vecd_tmp <- unique( sort( c( (res_master*round(colMeans(nn2$nn.dists)[kseq]/res_master)), (2:100)*res_master), decreasing=T))
    vecf_tmp <- pmin( res_master*min(c(env$nx,env$ny))/3, pmax( 1*res_master, round(vecd_tmp/2)))
    # vecd: sequence of spatial scales
    vecd <- unique( sort( c( vecd_tmp[which(!duplicated(vecf_tmp,fromLast=T) & vecd_tmp >= (2*res_master))],
                             vecd_tmp[which(!duplicated(vecf_tmp,fromLast=F) & vecd_tmp >= (2*res_master))]),
                             decreasing=T))
    nl   <- length(vecd)
    # vecf: sequence of aggregation factor to save computational time
    # aggregation factor is half the horizontal decorrelation length
    vecf <- round( vecd/(2*1000), 0)
    rm( vecf_tmp, vecd_tmp)
    # vece: observation to background ratios
    vece <- rep( argv$eps2, length(vecf))
    cat("+---------------------------------------------------------------+\n")
    cat("vecd vecf vece\n")
    cat(" m    -    -\n")
    for (l in 1:nl) cat(paste(vecd[l],vecf[l],vece[l],"\n"))
    cat("+---------------------------------------------------------------+\n")

    # Define relative anomalies 
    if (argv$use_relativeAnomalies) {
      if (any(yrf==0)) {
        cat("Warning: rescaling factor is equal to 0 for some grid point, this may cause weird patterns in the analysis")
        yrf[which(yrf==0)] <- 1
      }
      yo <- y_env$yo$value / yrf
    } else {
      yo <- y_env$yo$value
    }

    # Transformation 
    zero <- 0
    if ( argv$transf == "Box-Cox") {
      yo   <- boxcox( yo, argv$transf.boxcox_lambda)
      zero <- boxcox(  0, argv$transf.boxcox_lambda)
    } else if (argv$transf!="none") {
      boom("transformation not defined")
    }

    # Multi-scale OI loop
    envtmp$obs_x <- y_env$yo$x
    envtmp$obs_y <- y_env$yo$y
    for (l in 1:nl) {
      t0b <- Sys.time()
      cat( paste0( "scale # ",formatC(l,width=3,flag="0"),
                   " of ",nl,
                   " (",formatC(vecd[l],width=8,flag="0",format="d"),"m ",
                   "fact=",formatC(vecf[l],width=3,flag="0"),")"))
      # prepare grid for the l-th iteration
      if (l==nl | vecf[l]==1) {
        r <- env$rmaster
      } else {
        r <- aggregate(env$rmaster, fact=vecf[l], expand=T, na.rm=T)
      }
      xy.l    <- xyFromCell( r, 1:ncell(r))
      mask.l  <- which( !is.na( getValues(r)))
      envtmp$x <- xy.l[mask.l,1]
      envtmp$y <- xy.l[mask.l,2]
      envtmp$m_dim <- length(envtmp$x)
      envtmp$Eb  <- array( data=NA, dim=c( envtmp$m_dim, 1))
      rm(xy.l)
      # prepare background for the l-th iteration
      if (l==1) {
        # first one, background is the regional average of observations
        rb <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=r, field=yo, fun=mean, na.rm=T) 
      } else {
        rb <- resample( ra, r, method="bilinear")
        if ("scale" %in% argv$off_x.variables) 
          xl_tmp <- getValues( resample( rl, r, method="ngb"))[mask.l]
        # fill in NAs in the scale field
        if (any(is.na(xl_tmp))) {
          ix_to   <- which(  is.na(xl_tmp))
          ix_from <- which( !is.na(xl_tmp))
          xl_tmp[ix_to] <- ngb_gridpoint_by_gridpoint( x_from=envtmp$x[ix_from], 
                                                       y_from=envtmp$y[ix_from], 
                                                       val_from=xl_tmp[ix_from], 
                                                       x_to=envtmp$x[ix_to], 
                                                       y_to=envtmp$y[ix_to])
        }
        if (any(is.na(xl_tmp))) cat("Warning: NAs found is scale field\n") 
      }
      envtmp$Eb[,1] <- getValues(rb)[mask.l]
      # fill in NAs in the background
      if ( any( is.na(envtmp$Eb[,1]))) {
        ix_to   <- which(  is.na(envtmp$Eb[,1]))
        ix_from <- which( !is.na(envtmp$Eb[,1]))
        envtmp$Eb[ix_to,1] <- ngb_gridpoint_by_gridpoint( x_from=envtmp$x[ix_from], 
                                                          y_from=envtmp$y[ix_from], 
                                                          val_from=envtmp$Eb[ix_from,1], 
                                                          x_to=envtmp$x[ix_to], 
                                                          y_to=envtmp$y[ix_to])
      }
      rb[mask.l] <- envtmp$Eb[,1]
      yb <- extract( rb, cbind(y_env$yo$x,y_env$yo$y), method="bilinear")
      rm(rb)
      if (any(is.na(envtmp$Eb[,1]))) cat("Warning: xb is NA\n")
      if (any(is.na(yb))) cat("Warning: yb is NA\n")
      # compute correlations among observation locations
      D <- corr2d( par=vecd[l], label=argv$corrfun, values_as_globalVar=T)
      # prepare for OI faster
      diag(D) <- diag(D) + vece[l] * y_env$yo$ovarc
      envtmp$SRinv <- chol2inv(chol(D))
      envtmp$di <- array(data=(yo-yb),dim=c(y_env$yo$n,1))
      envtmp$SRinv_di <- crossprod( envtmp$SRinv, envtmp$di)
      # run OI gridpoint by gridpoint
      if (!is.na(argv$cores)) {
        res <- t( mcmapply( enoi_basicFaster_gridpoint_by_gridpoint,
                            1:envtmp$m_dim,
                            mc.cores=argv$cores,
                            SIMPLIFY=T,
                            MoreArgs = list( corr=argv$corrfun,
                                             dh=vecd[l],
                                             safecheck=T,
                                             idi=F)))
      # no-multicores
      } else {
        res <- t( mapply( enoi_basicFaster_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corrfun, 
                                           dh=vecd[l],
                                           safecheck=T,
                                           idi=F)))
      }
      xa.l <- res[,1]
      rm(res)
      ra <- r
      ra[] <- NA
      ra[mask.l] <- xa.l
      if ("scale" %in% argv$off_x.variables) {
        rl <- ra
        rl[mask.l] <- vecd[l]
        if (l>1) {
#          ixl <- which( abs(xa.l-envtmp$Eb[,1]) < 0.005)
          cond <- (envtmp$Eb[,1]==0 & xa.l==0) | 
                  ((abs(xa.l-envtmp$Eb[,1])/abs(envtmp$Eb[,1]))<0.001)
          ixl <- which(cond)
          if (length(ixl)>0) rl[mask.l[ixl]] <- xl_tmp[ixl]
        }
      }
      t1b <- Sys.time()
      cat( paste( " time", round(t1b-t0b,1), attr(t1b-t0b,"unit"), "\n"))
    } # end of multi-scale OI

    # back to precipitation values (from relative anomalies)
    if ( argv$transf == "Box-Cox") {
      xa <- tboxcox( xa.l, argv$transf.boxcox_lambda)
    } else if ( argv$transf != "none") {
      xa <- xa.l
    }
    xa_rel <- xa
    if (argv$use_relativeAnomalies) {
      xa <- xa * getValues(rrf)[mask.l]
    }

    # refine analysis
    cat("refine prec/no-prec borders and remove wet regions with no obs in them\n")
    for (frac_rr in c(100,50,10)) {
      if (frac_rr==100) xa[which(!is.na(xa) & xa<(y_env$rain/frac_rr))]<-0
      ra[mask.l]<-xa
      xa_aux<-getValues(ra)
      xa_aux[which(!is.na(xa_aux) & xa_aux<(y_env$rain/frac_rr))]<-NA
      ra[]<-xa_aux
      rclump<-clump(ra)
      oclump<-extract(rclump,cbind(y_env$yo$x[y_env$yo$ixwet],y_env$yo$y[y_env$yo$ixwet]))
      fr<-freq(rclump)
      # remove clumps of YESprec cells less than (4x4)km^2 or not including wet obs
      ix<-which(!is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=16) | !(fr[,1] %in% oclump)) )
      xa[which(getValues(rclump)[mask.l] %in% fr[ix,1])]<-0
      rm(xa_aux,rclump,oclump,fr,ix)
      ra[mask.l]<-xa
    }

  # If dry everywhere
  } else {
    ra <- env$rmaster
    ra[]     <- NA
    ra[mask] <- 0
    rl <- ra
  } # END IF dry everywhere

  # Initialization
  env$Xa <- array( data=NA, dim=c( env$ngrid, 1))
  env$Xa_rel <- array( data=NA, dim=c( env$ngrid, 1))
  env$Xidi <- array( data=NA, dim=c( env$ngrid, 1))
  env$Xscale <- array( data=NA, dim=c( env$ngrid, 1))
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, 1))
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, 1))

  # Output  
  env$Xa[,1] <- getValues(ra)[mask.l]
  env$Xa_rel[,1] <- xa_rel
  env$Xscale[,1] <- getValues(rl)[mask.l]
  env$Xidi[,1] <- xidi
  y_env$yo$value_a[,1] <- extract( ra, cbind( y_env$yo$x, y_env$yo$y), method="bilinear")
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$value_a[,1] <- extract( ra, cbind( y_env$yov$x, y_env$yov$y), method="simple")

  # Bye-bye
  t1a <- Sys.time()
  cat( paste( "OI multiscale developed for seNorge, total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))

}
