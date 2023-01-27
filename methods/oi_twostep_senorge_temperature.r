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
  y_env$super_yb  <- list()
  y_env$super_ya  <- list()
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

  do_cv <- FALSE
  if (env$cv_mode | env$cv_mode_random) do_cv <- TRUE
  do_xa <- FALSE
  if ( !is.na( argv$off_x)) do_xa <- TRUE
  do_ya <- argv$twostep_superobbing

  cat(paste("Perform analysis over observation points is always TRUE\n"))
  cat(paste("(but if superobbing is requested then we need to perform analysis over actual observation points. "))
  cat(paste("Superobbing is",do_ya,")\n"))
  cat(paste("Perform cv-analysis over selected observation points is",do_cv,"\n"))
  cat(paste("Perform analysis over grid points is",do_xa,"\n"))

  # Regional backgrounds
  cat("+--------------------------------------------------------------+\n")
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

  # Alternative way to find the centroids via clustering (instead of using a coarser grid)
#  obs_clusters <- kmeans( cbind(y_env$super_yo$x,y_env$super_yo$y), centers=argv$oi2step.bg_centroids_nclusters)
#  xr <- obs_clusters$centers[,1]
#  yr <- obs_clusters$centers[,2]
#  nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
#                     query = cbind( xr, yr), 
#                     k = argv$oi2step.bg_centroids_nobsmin, 
#                     searchtype = "radius", 
#                     radius = argv$oi2step.bg_centroids_buffer)
#  y_env$centroids$i <- integer(0)
#  for (i in 1:length(xr)) {
##    if (length( which(nn2$nn.idx[i,]!=0)) < argv$oi2step.bg_centroids_nobsmin) next
#    if ( any( nn2$nn.idx[i,] == 0)) next
#    y_env$centroids$i <- c( y_env$centroids$i, i)
#  } 
#  y_env$centroids$n <- length(y_env$centroids$i)
#
#  cat(paste("the master grid has been divided in",
#              argv$oi2step.bg_centroids_nclusters,"clusters\n"))
#  cat(paste("sub-regional area extensions (length x (m),length y (m))=",
#        round(res(y_env$centroids$r)[1]),round(res(y_env$centroids$r)[2])),"\n")
#  cat(paste("reference (horizontal) length scale to weight the sub-regional backgrounds (m)=",round(mean(res(y_env$centroids$r))),"\n"))
#  cat(paste("# sub-regional centroids", y_env$centroids$n, "\n"))

  # save results in a global variable
  y_env$centroids$x <- xr[y_env$centroids$i]
  y_env$centroids$y <- yr[y_env$centroids$i]

  # clean memory
  rm(nn2, r_notmasked, xr, yr)

  # estimate vertical profiles at centroids (sub-regional)
  envtmp$nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
                            query = cbind( y_env$centroids$x, y_env$centroids$y), 
                            k = y_env$super_yo$n, 
                            searchtype = "radius", 
                            radius = argv$oi2step.bg_centroids_buffer)

  # multicores 
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
  y_env$centroids$vert_prof <- res

  # NOTE: move to argparser
  argv$oi2step.bg_blending_deltam <- 10000
  argv$oi2step.bg_blending_corr <- "soar"
  argv$oi2step.bg_blending_dh <- 100000

  # blend at observation points the sub-regional profiles into a regional one

  cat("At observation points, blend sub-regional profiles into a regional one\n")
  # helpers
  envtmp$n1 <- length( envtmp$ix1 <- which( y_env$centroids$vert_prof[,1] == 1))
  envtmp$n2 <- length( envtmp$ix2 <- which( y_env$centroids$vert_prof[,1] == 2))

  # initialization
  y_env$super_yb$value <- array( data=NA, dim=c(y_env$super_yo$n,1))
  envtmp$x <- y_env$super_yo$x
  envtmp$y <- y_env$super_yo$y
  envtmp$z <- y_env$super_yo$z
  envtmp$m_dim <- y_env$super_yo$n 
  envtmp$nn2 <- nn2( cbind( y_env$centroids$x, y_env$centroids$y), 
                            query = cbind( envtmp$x, envtmp$y), 
                            k = y_env$centroids$n, 
                            searchtype = "radius", 
                            radius = 99999999)
  # blending based on IDI
  if (!is.na(argv$cores)) {
    res <- t( mcmapply( blend_vertical_profiles_senorge2018,
                        1:envtmp$m_dim,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                         dh=argv$oi2step.bg_blending_dh)))
  # no-multicores
  } else {
    res <- t( mapply( blend_vertical_profiles_senorge2018,
                      1:envtmp$m_dim,
                      SIMPLIFY=T,
                      MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                       dh=argv$oi2step.bg_blending_dh)))
  }
  # save results in background data structure
  y_env$super_yb$value[,1] <- res

  cat(paste("RMS(yo-yb)=",round(sqrt(mean((y_env$super_yo$value-y_env$super_yb$value[,1])**2)),2),"degC\n"))
  
  # blend at observation points (if superobbing is TRUE) the sub-regional profiles into a regional one
  if (do_ya) {
    cat("At actual observation points (not at superobbing points), blend sub-regional profiles into a regional one\n")

    # initialization
    y_env$yo$value_b <- array( data=NA, dim=c(y_env$yo$n,1))
    envtmp$x <- y_env$yo$x
    envtmp$y <- y_env$yo$y
    envtmp$z <- y_env$yo$z
    envtmp$m_dim <- y_env$yo$n 
    envtmp$nn2 <- nn2( cbind( y_env$centroids$x, y_env$centroids$y), 
                              query = cbind( envtmp$x, envtmp$y), 
                              k = y_env$centroids$n, 
                              searchtype = "radius", 
                              radius = 99999999)
    # blending based on IDI
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( blend_vertical_profiles_senorge2018,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                           dh=argv$oi2step.bg_blending_dh)))
    # no-multicores
    } else {
      res <- t( mapply( blend_vertical_profiles_senorge2018,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                         dh=argv$oi2step.bg_blending_dh)))
    }
    # save results in background data structure
    y_env$yo$value_b[,1] <- res

    cat(paste("RMS(yo-yb)=",round(sqrt(mean((y_env$yo$value-y_env$yo$value_b[,1])**2)),2),"degC\n"))
  } # end if do_cv

  # blend at cv-points the sub-regional profiles into a regional one
  if (do_cv) {
    cat("At observation points, blend sub-regional profiles into a regional one\n")

    # initialization
    y_env$yov$value_b <- array( data=NA, dim=c(y_env$yov$n,1))
    envtmp$x <- y_env$yov$x
    envtmp$y <- y_env$yov$y
    envtmp$z <- y_env$yov$z
    envtmp$m_dim <- y_env$yov$n 
    envtmp$nn2 <- nn2( cbind( y_env$centroids$x, y_env$centroids$y), 
                              query = cbind( envtmp$x, envtmp$y), 
                              k = y_env$centroids$n, 
                              searchtype = "radius", 
                              radius = 99999999)
    # blending based on IDI
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( blend_vertical_profiles_senorge2018,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                           dh=argv$oi2step.bg_blending_dh)))
    # no-multicores
    } else {
      res <- t( mapply( blend_vertical_profiles_senorge2018,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                         dh=argv$oi2step.bg_blending_dh)))
    }
    # save results in background data structure
    y_env$yov$value_b[,1] <- res

    cat(paste("RMS(yo-yb)=",round(sqrt(mean((y_env$yov$value-y_env$yov$value_b[,1])**2)),2),"degC\n"))
  } # end if do_cv

  # blend at gridpoints the sub-regional profiles into a regional one
  if (do_xa) {
  
    cat("At grid points, blend sub-regional profiles into a regional one\n")

    # initialization
    m <- 0
    env$Xb <- array( data=NA, dim=c(env$ngrid,1))
    while (m <= env$ngrid) {
      m1 <- m + 1; m2 <- min( c( m + argv$oi2step.bg_blending_deltam, env$ngrid))
      envtmp$x <- env$xgrid[env$mask][m1:m2]
      envtmp$y <- env$ygrid[env$mask][m1:m2]
      envtmp$z <- getValues(env$rdem)[env$mask][m1:m2]
      envtmp$m_dim <- m2 - m1 + 1
      #
      envtmp$nn2 <- nn2( cbind( y_env$centroids$x, y_env$centroids$y), 
                                query = cbind( envtmp$x, envtmp$y), 
                                k = y_env$centroids$n, 
                                searchtype = "radius", 
                                radius = 99999999)
      # blending based on IDI
      if (!is.na(argv$cores)) {
        res <- t( mcmapply( blend_vertical_profiles_senorge2018,
                            1:envtmp$m_dim,
                            mc.cores=argv$cores,
                            SIMPLIFY=T,
                            MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                             dh=argv$oi2step.bg_blending_dh)))
      # no-multicores
      } else {
        res <- t( mapply( blend_vertical_profiles_senorge2018,
                          1:envtmp$m_dim,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$oi2step.bg_blending_corr,
                                           dh=argv$oi2step.bg_blending_dh)))
      }
      # save results in background data structure
      env$Xb[m1:m2,1] <- res
      # next bunch of gridpoints
      m <- m + argv$oi2step.bg_blending_deltam
    } # end loop over gridpoints
  
    cat( paste( "range(xb) (min max, degC)=",
                round(min(env$Xb[,1],na.rm=T),2),
                round(max(env$Xb[,1],na.rm=T),2),"\n"))
    cat(paste("is any of xb set to NA?",any(is.na(env$Xb[,1])),"\n"))
  } # end if do_xa

  # NOTE: adaptive Dh must be implemented
  envtmp$dh_obs <- rep( 60000, y_env$yo$n)
  if (do_ya) envtmp$dh_yaobs <- rep( 60000, y_env$yo$n) 
  if (do_cv) envtmp$dh_cvobs <- rep( 60000, y_env$yov$n)
  if (do_xa) envtmp$dh_grid  <- rep( 60000, env$ngrid)

  # NOTE: move to argparser
  argv$oi2step.loop_deltam <- 10000
  argv$oi2step.analysis_dz <- c(800,600,400,200)
  argv$oi2step.analysis_eps2 <- c(0.5,0.5,0.5,0.25)
  argv$oi2step.analysis_lafmin <- 0 
  argv$oi2step.analysis_nclose <- 50
  argv$oi2step.analysis_corr <- c("soar","gaussian","linear")

  # Analysis 
  cat("+--------------------------------------------------------------+\n")
  cat("Analysis\n")

  # initalization
  nzloop <- length(argv$oi2step.analysis_dz)
  # observation points
  y_env$super_ya$value <- array( data=NA, dim=c(y_env$super_yo$n,1))
  y_env$super_ya$idi_loop <- array( data=NA, dim=c(y_env$super_yo$n,nzloop))
  y_env$super_yb$value_loop <- array( data=NA, dim=c(y_env$super_yo$n,nzloop))
  y_env$super_yb$value_loop[,1] <- y_env$super_yb$value[,1]
  # observation points (if different from superobbing)
  if (do_ya) {
    y_env$yo$value_a <- array( data=NA, dim=c(y_env$super_yo$n,1))
    y_env$yo$idi_loop <- array( data=NA, dim=c(y_env$super_yo$n,nzloop))
    y_env$yo$value_b_loop <- array( data=NA, dim=c(y_env$yo$n,nzloop))
    y_env$yo$value_b_loop[,1] <- y_env$yo$value_b[,1]
  }
  # cv points
  if (do_cv) {
    y_env$yov$value_a <- array( data=NA, dim=c(y_env$super_yo$n,1))
    y_env$yov$idi_loop <- array( data=NA, dim=c(y_env$super_yo$n,nzloop))
    y_env$yov$value_b_loop <- array( data=NA, dim=c(y_env$yov$n,nzloop))
    y_env$yov$value_b_loop[,1] <- y_env$yov$value_b[,1]
  }
  # grid points
  if (do_xa) {
    env$Xa <- array( data=NA, dim=c(env$ngrid,1))
    env$Xidi <- array( data=NA, dim=c(env$ngrid,1))
    envtmp$Eb_loop <- array( data=NA, dim=c(env$ngrid,nzloop))
    envtmp$Eb_loop[,1] <- env$Xb[,1] 
    envtmp$Xidi_loop <- array( data=NA, dim=c(env$ngrid,nzloop))
  }

  # Matrix of distances between all stations
  envtmp$dist2     <- outer(y_env$super_yo$x,y_env$super_yo$x,FUN="-")**2.+
                      outer(y_env$super_yo$y,y_env$super_yo$y,FUN="-")**2.
  # Matrix of elevation differences between all stations
  envtmp$dist2_z   <- outer(y_env$super_yo$z,y_env$super_yo$z,FUN="-")**2.
  # Matrix of distances between all stations
  envtmp$dist2_laf <- outer(y_env$super_yo$laf,y_env$super_yo$laf,FUN="-")**2.

  envtmp$obs_z   <- y_env$super_yo$z
  envtmp$obs_laf <- y_env$super_yo$laf
 
  # Loop over vertical decorrelation lengths
  for (zloop in 1:nzloop) {

    cat(paste("set vertical decorrelation length=",argv$oi2step.analysis_dz[zloop],"\n"))

    # analysis at observation points
    envtmp$x <- y_env$super_yo$x
    envtmp$y <- y_env$super_yo$y 
    envtmp$z <- y_env$super_yo$z
    envtmp$laf <- y_env$super_yo$laf
    envtmp$m_dim <- y_env$super_yo$n
    envtmp$par <- cbind( envtmp$dh_obs, 
                         rep(argv$oi2step.analysis_dz[zloop],y_env$super_yo$n),
                         rep(argv$oi2step.analysis_lafmin,y_env$super_yo$n))
    envtmp$eps2 <- rep(argv$oi2step.analysis_eps2[zloop],envtmp$m_dim)
    envtmp$nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
                              query = cbind( envtmp$x, envtmp$y), 
                              k = argv$oi2step.analysis_nclose, 
                              searchtype = "radius", 
                              radius = 99999999)
    envtmp$Eb <- array( data=y_env$super_yb$value_loop[,zloop], dim=c(envtmp$m_dim,1))
    envtmp$D  <- array( data=(y_env$super_yo$value - y_env$super_yb$value_loop[,zloop]), dim=c(envtmp$m_dim,1))
    # analysis based on IDI
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
    # no-multicores
    } else {
      res <- t( mapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
    }
    y_env$super_ya$value[,1]        <- res[,1]
    y_env$super_ya$idi_loop[,zloop] <- res[,2]
    cat(paste("RMS(yo-ya)=",round(sqrt(mean((y_env$super_yo$value-y_env$super_ya$value[,1])**2)),2),"degC\n"))

    # analysis at observation points (if different from superobbing)
    if (do_ya) {
      envtmp$x <- y_env$yo$x
      envtmp$y <- y_env$yo$y 
      envtmp$z <- y_env$yo$z
      envtmp$laf <- y_env$yo$laf
      envtmp$m_dim <- y_env$yo$n
      envtmp$par <- cbind( envtmp$dh_yaobs, 
                           rep(argv$oi2step.analysis_dz[zloop],y_env$yo$n),
                           rep(argv$oi2step.analysis_lafmin,y_env$yo$n))
      envtmp$eps2 <- rep(argv$oi2step.analysis_eps2[zloop],envtmp$m_dim)
      envtmp$nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
                                query = cbind( envtmp$x, envtmp$y), 
                                k = argv$oi2step.analysis_nclose, 
                                searchtype = "radius", 
                                radius = 99999999)
      envtmp$Eb <- array( data=y_env$yo$value_b_loop[,zloop], dim=c(envtmp$m_dim,1))
      # analysis based on IDI
      if (!is.na(argv$cores)) {
        res <- t( mcmapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                            1:envtmp$m_dim,
                            mc.cores=argv$cores,
                            SIMPLIFY=T,
                            MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
      # no-multicores
      } else {
        res <- t( mapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
      }
      y_env$yo$value_a[,1]      <- res[,1]
      y_env$yo$idi_loop[,zloop] <- res[,2]
      cat(paste("not-at-superobbing-sites RMS(yo-ya)=",round(sqrt(mean((y_env$super_yo$value-y_env$super_ya$value[,1])**2)),2),"degC\n"))
    } # end if do_ya

    # analysis at cv-observation points
    if (do_cv) {
      envtmp$x <- y_env$yov$x
      envtmp$y <- y_env$yov$y 
      envtmp$z <- y_env$yov$z
      envtmp$laf <- y_env$yov$laf
      envtmp$m_dim <- y_env$yov$n
      envtmp$par <- cbind( envtmp$dh_cvobs, 
                           rep(argv$oi2step.analysis_dz[zloop],y_env$yov$n),
                           rep(argv$oi2step.analysis_lafmin,y_env$yov$n))
      envtmp$eps2 <- rep(argv$oi2step.analysis_eps2[zloop],envtmp$m_dim)
      envtmp$nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
                                query = cbind( envtmp$x, envtmp$y), 
                                k = argv$oi2step.analysis_nclose, 
                                searchtype = "radius", 
                                radius = 99999999)
      envtmp$Eb <- array( data=y_env$yov$value_b_loop[,zloop], dim=c(envtmp$m_dim,1))
      # analysis based on IDI
      if (!is.na(argv$cores)) {
        res <- t( mcmapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                            1:envtmp$m_dim,
                            mc.cores=argv$cores,
                            SIMPLIFY=T,
                            MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
      # no-multicores
      } else {
        res <- t( mapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
      }
      y_env$yov$value_a[,1]      <- res[,1]
      y_env$yov$idi_loop[,zloop] <- res[,2]
      cat(paste("RMS(yo-yav)=",round(sqrt(mean((y_env$super_yo$value-y_env$super_ya$value[,1])**2)),2),"degC\n"))
    } # end if do_cv

    # analysis at grid points
    if (do_xa) {
      m <- 0
      count_loop <- 1
      while (m <= env$ngrid) {
        cat( paste( count_loop, "/", ceiling(env$ngrid/argv$oi2step.loop_deltam)))
        m1 <- m + 1; m2 <- min( c( m + argv$oi2step.loop_deltam, env$ngrid))
        envtmp$x <- env$xgrid[env$mask][m1:m2]
        envtmp$y <- env$ygrid[env$mask][m1:m2]
        envtmp$z <- getValues(env$rdem)[env$mask][m1:m2]
        envtmp$laf <- getValues(env$rlaf)[env$mask][m1:m2]
        envtmp$m_dim <- m2 - m1 + 1
        envtmp$par <- cbind( envtmp$dh_grid[m1:m2], 
                             rep(argv$oi2step.analysis_dz[zloop],envtmp$m_dim),
                             rep(argv$oi2step.analysis_lafmin,envtmp$m_dim))
        envtmp$eps2 <- rep(argv$oi2step.analysis_eps2[zloop],envtmp$m_dim)
        envtmp$nn2 <- nn2( cbind( y_env$super_yo$x, y_env$super_yo$y), 
                                  query = cbind( envtmp$x, envtmp$y), 
                                  k = argv$oi2step.analysis_nclose, 
                                  searchtype = "radius", 
                                  radius = 99999999)
        envtmp$Eb <- array( data=envtmp$Eb_loop[m1:m2,zloop], dim=c(envtmp$m_dim,1))
        # analysis based on IDI
        if (!is.na(argv$cores)) {
          res <- t( mcmapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                              1:envtmp$m_dim,
                              mc.cores=argv$cores,
                              SIMPLIFY=T,
                              MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
        # no-multicores
        } else {
          res <- t( mapply( oi_senorge_temperature_gridpoint_by_gridpoint,
                            1:envtmp$m_dim,
                            SIMPLIFY=T,
                            MoreArgs = list( corr=argv$oi2step.analysis_corr, idi=T)))
        }
        # save results in background data structure
        env$Xa[m1:m2,1]               <- res[,1]
        envtmp$Xidi_loop[m1:m2,zloop] <- res[,2]
        # next bunch of gridpoints
        m <- m + argv$oi2step.loop_deltam
        count_loop <- count_loop + 1
      } # end loop over gridpoints
      cat( paste( "range(xa-xb) (min max, degC)=",
                  round(min(env$Xa[,1]-envtmp$Eb_loop[,zloop],na.rm=T),2),
                  round(max(env$Xa[,1]-envtmp$Eb_loop[,zloop],na.rm=T),2),"\n"))
      cat(paste("is any of xa set to NA?",any(is.na(env$Xa)),"\n"))

    } # end if do_xa

    # prepare for next iteration
    if (zloop < nzloop) {
      y_env$super_yb$value_loop[,zloop+1] <- y_env$super_ya$value[,1] 
      if (do_ya) y_env$yo$value_b_loop[,zloop+1] <- y_env$yo$value_a[,1]
      if (do_cv) y_env$yov$value_b_loop[,zloop+1] <- y_env$yov$value_a[,1]
      if (do_xa) envtmp$Eb_loop[,zloop+1] <- env$Xa[,1]
    }
  } # end loop over elevation decorellation parameters

  save(file="tmp.rdata",envtmp,y_env,res,env,argv,fg_env)
  # save special results for Output
  if (do_xa) {

    # Initialization
    env$Xidi <- array( data=NA, dim=c( env$ngrid, 1))
    env$Xscale <- array( data=NA, dim=c( env$ngrid, 1))

    # Output  
    env$Xscale[,1] <- envtmp$dh_grid
    env$Xidi[,1] <- rowSums(envtmp$Xidi_loop)

    # save background in the output data structures
    env$k_dim <- 1
    fg_env$ixs <- 1
    fg_env$ixf <- 1
    fg_env$ixe <- 1
    fg_env$fg[[1]]$r_main <- env$rmaster
    fg_env$fg[[1]]$r_main[] <- NA
    fg_env$fg[[1]]$r_main[env$mask] <- env$Xb[,1]

  }

  # Bye-bye
  t1a <- Sys.time()
  cat( paste( "OI two-step developed for seNorge, total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))

}
