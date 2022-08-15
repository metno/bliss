#+ Change-of-Resolution Ensemble Kalman Smoother
corenks <- function( argv, y_env, fg_env, env, dir_plot=NA) {
#
#------------------------------------------------------------------------------

argv$corenks_ididense <- 0.8
argv$corenks_pmax <- 30
argv$corenks_corrfun <- "toar"
y_env$rain<-0.1

  t0a <- Sys.time()

  cat( "-- corenks --\n")

  # number of background ensemble members
  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

  # set domain
  xmn <- as.numeric(argv$grid_master.x1) - as.numeric(argv$grid_master.resx)/2
  xmx <- as.numeric(argv$grid_master.xn) + as.numeric(argv$grid_master.resx)/2
  ymn <- as.numeric(argv$grid_master.y1) - as.numeric(argv$grid_master.resy)/2
  ymx <- as.numeric(argv$grid_master.yn) + as.numeric(argv$grid_master.resy)/2

  nx <- ncol( env$rmaster)
  ny <- nrow( env$rmaster)

  # definition of the tree-structure, coarser spatial level has a resolution of 2**jmax cells
  jmax <- floor( log( max(nx,ny), 2))

  # -~- Observations -~-
  # multi-resolution tree-structured data model
  mrtree <- list()
  # multi-resolution observation, tree-structured data model
  mrobs <- list()
  # loop over the spatial levels in the tree-structured data model
  for (j in 1:jmax) {
    mrtree$raster[[j]] <- list()
    # j=1 -> finest level has the same grid as the master grid (nx,ny)
    if ( j == 1) {
      # initialization from merged observations(they covers the wider region where we have observations)
      mrtree$raster[[j]]$r <- env$rmaster
      mrtree$raster[[j]]$r[] <- NA
      robs <- mrtree$raster[[j]]$r
      robs[] <- env$mergeobs$value
      ridi <- mrtree$raster[[j]]$r
      ridi[] <- env$mergeobs$idi
      ridi[is.na(ridi)] <- 0
    # j-th level has a grid of approximately (nx/2**j,ny/2**j)
    } else {
      raux <- mrtree$raster[[j-1]]$r
      raux[] <- mrobs$val_all[[j-1]]
      robs <- aggregate( raux, fact=2, fun="mean", expand=T, na.rm=T)
      raux[] <- mrobs$idi[[j-1]]
      ridi <- aggregate( raux, fact=2, fun="mean", expand=T, na.rm=T)
      mrtree$raster[[j]]$r <- ridi
      mrtree$raster[[j]]$r[] <- NA
    }
    # select only gridpoints where observations are ok (IDI is larger than a threshold)
    ix <- which( getValues(ridi) >= argv$corenks_ididense & 
                 !is.na( getValues(robs)))
    # stop if no observations found 
    if ( length(ix) == 0) { 
      jstop <- j-1
      break 
    # else save observations in the tree structure
    } else {
      mrtree$m_dim[[j]] <- ncell(mrtree$raster[[j]]$r)
      mrtree$res[[j]] <- res(mrtree$raster[[j]]$r)
      mrtree$mean_res[[j]]  <- mean(mrtree$res[[j]])
      xy <- xyFromCell( mrtree$raster[[j]]$r, 1:mrtree$m_dim[[j]])
      mrtree$x[[j]] <- xy[,1]
      mrtree$y[[j]] <- xy[,2]
      # idi for all gridpoints
      mrobs$idi[[j]] <- getValues(ridi)
      # observed values only for selected gridpoints
      mrobs$ix[[j]] <- ix
      mrobs$d_dim[[j]] <- length(ix)
      mrobs$val[[j]] <- getValues(robs)[ix]
      mrobs$x[[j]] <- xy[ix,1]
      mrobs$y[[j]] <- xy[ix,2]
      mrobs$val_all[[j]] <- getValues(robs)
    }
  } # end loop over spatal levels
  # safe-check, exit when no ok observations found 
  if (jstop == 0) return(NULL)

  # -~- Fine-to-Coarse CorEnKS sweep -~-
  # tree-structured data model (we identify this sweep with filtering step in EnKS)
  mrenkf <- list()
  # loop over spatial levels, from fine (j=1) to coarse (j=jstop, i.e. coarsest level with observations we can trust)
  for (j in 1:jstop) {
    mrenkf$data[[j]] <- list()
    # loop over ensemble members
    for (e in 1:env$k_dim) {
      # finest level, get the background from the ensemble on the master grid
      if ( j == 1) {
        i <- fg_env$ixs[e]
        r <- subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i])
      # jth-level, perform the saptial aggregation
      } else {
        s <- mrtree$raster[[j-1]]$r
        s[] <- mrenkf$data[[j-1]]$Ea[,e]
        r <- aggregate( s, fact=2, fun="mean", expand=T, na.rm=T)
      }
      # initializations (only once per level)
      if ( e == 1) {
        mrenkf$data[[j]]$Eb  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrenkf$data[[j]]$Xb  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrenkf$data[[j]]$Y   <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
        mrenkf$data[[j]]$HEb <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
      }
      # ensemble members on the grid
      mrenkf$data[[j]]$Eb[,e] <- getValues(r)
      # ensemble members at observation locations
      mrenkf$data[[j]]$HEb[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
    } # END loop over ensemble members
    # ensemble mean on the grid
    mrenkf$data[[j]]$xb <- rowMeans(mrenkf$data[[j]]$Eb)
    # safe-check and selection of gridpoints having all values not NAs
    if ( (mrenkf$m_dim[[j]] <- length( ix <- which( !is.na( mrenkf$data[[j]]$xb)))) == 0) return(NULL)
    # ensemble mean at observation loactions
    mrenkf$data[[j]]$Hxb <- rowMeans(mrenkf$data[[j]]$HEb)
    # safe-check and selection of observations having all background values not NAs
    if ( (mrenkf$d_dim[[j]] <- length( iy <- which( !is.na( mrenkf$data[[j]]$Hxb)))) == 0) return(NULL)
    # ensemble anomalies 
    for (e in 1:env$k_dim) { 
      mrenkf$data[[j]]$Xb[,e] <- 1/sqrt(env$k_dim-1) * (mrenkf$data[[j]]$Eb[,e] - mrenkf$data[[j]]$xb)
      r[] <- mrenkf$data[[j]]$Xb[,e]
      mrenkf$data[[j]]$Y[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
    }
    # init structure used for EnKF(EnOI)
    envtmp$x <- mrtree$x[[j]][ix]
    envtmp$y <- mrtree$y[[j]][ix]
    envtmp$obs_x <- mrobs$x[[j]][iy]
    envtmp$obs_y <- mrobs$y[[j]][iy]
    envtmp$m_dim <- mrenkf$m_dim[[j]]
    envtmp$k_dim <- env$k_dim
    envtmp$obs_val <- mrobs$val[[j]][iy]
    envtmp$Eb <- mrenkf$data[[j]]$Eb[ix,]
    envtmp$HEb <- mrenkf$data[[j]]$HEb[iy,]
    envtmp$eps2 <- rep( 0.1, envtmp$m_dim)
    envtmp$D <- envtmp$obs_val - envtmp$HEb
    # helper to get the neighbours
    envtmp$nn2 <- nn2( cbind(envtmp$obs_x,envtmp$obs_y), 
                       query = cbind(envtmp$x,envtmp$y), 
                       k = min( c(argv$corenks_pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( enoi_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corenks_corrfun,
                                           dh=mrtree$mean_res[[j]])))
    # no-multicores
    } else {
      res <- t( mapply( enoi_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corenks_corrfun, 
                                         dh=mrtree$mean_res[[j]])))
    }
    mrenkf$data[[j]]$Ea <- res[,1:env$k_dim]
    if (!is.na(y_env$rain)) mrenkf$data[[j]]$Ea[mrenkf$data[[j]]$Ea<y_env$rain] <- 0
  } # END loop over spatial levels
  
  # -~- Coarse-to-fine CorEnKS sweep -~-
  # tree-structured data model (we identify this sweep with smoothing step in EnKS)
  mrenks <- list()
  for (j in (jstop-1):1) {
    mrenks$data[[j]] <- list()
    for (e in 1:env$k_dim) {
      s <- mrtree$raster[[j+1]]$r
      r <- mrtree$raster[[j]]$r
      if ( j == (jstop-1)) {
        s[] <- mrenkf$data[[jstop]]$Ea[,e]
      } else {
        s[] <- mrenks$data[[j+1]]$Ea[,e]
      }
      r <- resample( s, r, method="bilinear", na.rm=T)
      # initializations (only once per level)
      if ( e == 1) {
        mrenks$data[[j]]$Eb  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrenks$data[[j]]$Xb  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrenks$data[[j]]$Y   <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
        mrenks$data[[j]]$HEb <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
      }
      # ensemble members on the grid
      mrenks$data[[j]]$Eb[,e] <- getValues(r)
      # ensemble members at observation locations
      mrenks$data[[j]]$HEb[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
    } # END loop over ensemble members
    # ensemble mean on the grid
    mrenks$data[[j]]$xb <- rowMeans(mrenks$data[[j]]$Eb)
    # safe-check and selection of gridpoints having all values not NAs
    if ( (mrenks$m_dim[[j]] <- length( ix <- which( !is.na( mrenks$data[[j]]$xb)))) == 0) return(NULL)
    # ensemble mean at observation loactions
    mrenks$data[[j]]$Hxb <- rowMeans(mrenks$data[[j]]$HEb)
    # safe-check and selection of observations having all background values not NAs
    if ( (mrenks$d_dim[[j]] <- length( iy <- which( !is.na( mrenks$data[[j]]$Hxb)))) == 0) return(NULL)
    # ensemble anomalies 
    for (e in 1:env$k_dim) { 
      mrenks$data[[j]]$Xb[,e] <- 1/sqrt(env$k_dim-1) * (mrenks$data[[j]]$Eb[,e] - mrenks$data[[j]]$xb)
      r[] <- mrenks$data[[j]]$Xb[,e]
      mrenks$data[[j]]$Y[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
    }
    # init structure used for EnKF(EnOI)
    envtmp$x <- mrtree$x[[j]][ix]
    envtmp$y <- mrtree$y[[j]][ix]
    envtmp$obs_x <- mrobs$x[[j]][iy]
    envtmp$obs_y <- mrobs$y[[j]][iy]
    envtmp$m_dim <- mrenks$m_dim[[j]]
    envtmp$k_dim <- env$k_dim
    envtmp$obs_val <- mrobs$val[[j]][iy]
    envtmp$Eb <- mrenks$data[[j]]$Eb[ix,]
    envtmp$HEb <- mrenks$data[[j]]$HEb[iy,]
    envtmp$eps2 <- rep( 0.1, envtmp$m_dim)
    envtmp$D <- envtmp$obs_val - envtmp$HEb
    # helper to get the neighbours
    envtmp$nn2 <- nn2( cbind(envtmp$obs_x,envtmp$obs_y), 
                       query = cbind(envtmp$x,envtmp$y), 
                       k = min( c(argv$corenks_pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( enoi_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corenks_corrfun,
                                           dh=mrtree$mean_res[[j]])))
    # no-multicores
    } else {
      res <- t( mapply( enoi_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corenks_corrfun, 
                                         dh=mrtree$mean_res[[j]])))
    }
    mrenks$data[[j]]$Ea <- res[,1:env$k_dim]
    if (!is.na(y_env$rain)) mrenks$data[[j]]$Ea[mrenks$data[[j]]$Ea<y_env$rain] <- 0
  }
save(file="tmp.rdata",mrobs,mrtree,mrenkf,mrenks)
#e<-6; j<-6; r<-mrtree$raster[[j]]$r; r[] <- mrenks$data[[j]]$Ea[,e]; image(r,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
q()
  ixb <- which( getValues(rfidi) >= argv$wise_preproc_idisparse)
  ixa <- which( getValues(rfidi) >= argv$wise_preproc_ididense)
  rfobsb <- rdyad
  rfobsb[ixb] <- rfobs[ixb]

  lXb_dyad <- list()
  for (e in 1:env$k_dim) {
    i <- fg_env$ixs[e]
    # interpolate onto the dyadic grid
    # background at grid  points
    rfxb <- subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i])
save(file="tmp.rdata",rfxb,robs,robsidi)
q()

    rfxb[rfxb<y_env$rain] <- 0
#    Xb_dyad_original[,e]  <- getValues(rfxb)
#    Xb_dyad[,e]  <- getValues(rfxb)
    lXb_dyad[[e]] <- list()
    lXb_dyad[[e]]$mr[[1]] <- getValues(rfxb)
    for (j in 1:nn) {
      cat(paste("e j",e,j,"\n"))
      dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=j, boundary=env$boundary)
      lXb_dyad[[e]]$mr[[j+1]] <- dwt[[3*j+1]]
    }
  }

  # set dyadic domain
  xmn <- as.numeric(argv$grid_master.x1) - as.numeric(argv$grid_master.resx)/2
  xmx <- as.numeric(argv$grid_master.xn) + as.numeric(argv$grid_master.resx)/2
  ymn <- as.numeric(argv$grid_master.y1) - as.numeric(argv$grid_master.resy)/2
  ymx <- as.numeric(argv$grid_master.yn) + as.numeric(argv$grid_master.resy)/2

  nx <- ncol( env$rmaster)
  ny <- nrow( env$rmaster)

  n <- ceiling( log( max(nx,ny), 2))
  rdyad <- raster( extent( xmn, xmx, ymn, ymx), 
                   ncol=2**n, nrow=2**n, crs=argv$grid_master.proj4)
  rdyad[] <- 0

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))

  # ---~--------------
  # observations on the dyadic grid

argv$wise_preproc_ididense <- 0.8
argv$wise_preproc_idisparse <- 0.1

  r <-  env$rmaster
  
  # rfobs covers the wider region where we have observations
  r[] <- env$mergeobs$value
  rfobs <- resample( r, rdyad, method="bilinear", na.rm=T)
  rfobs[is.na(rfobs)] <- 0
  if (!is.na(y_env$rain)) rfobs[rfobs<y_env$rain] <- 0

  # rfidi covers wider areas than rfobs (=1 obs-dense; =0 obs-void)
  r[] <- env$mergeobs$idi
  rfidi <- resample( r, rdyad, method="bilinear", na.rm=T)
  rfidi[rfidi<0] <- 0
  rfidi[rfidi>1] <- 1
  rfidi[is.na(rfidi)] <- 0

  # helper
  ixb <- which( getValues(rfidi) >= argv$wise_preproc_idisparse)
  ixa <- which( getValues(rfidi) >= argv$wise_preproc_ididense)
  rfobsb <- rdyad
  rfobsb[ixb] <- rfobs[ixb]

  # idi region a
  rfidia <- rdyad
  rfidia[ixa] <- 1
  # idi region b
  rfidib <- rdyad
  rfidib[ixb] <- 1
  # =1 obs-void; =0 obs_dense
  rfidic <- 1 - rfidib

#  # rfobsall covers wider areas than rfobs
#  rfobsb[rfidi<argv$wise_preproc_idisparse] <- 0
#  rfobsall <- resample( r, rdyad, method="bilinear", na.rm=T)
#  rfobsall[rfobsall<y_env$rain] <- 0
#
#  # helper
#  c_xy <- which( !is.na( getValues(rfobs)))
#  # save for future use
#  env$rfobs <- rfobs
#  env$rfidi <- rfidi
#  env$rfobsall <- rfobsall
#
#  # paddle with 0s
#  rfobs[is.na(rfobs)] <- 0
#  rfidi[is.na(rfidi)] <- 0
#  rfobsall[is.na(rfobsall)] <- 0
#
#  vidi <- getValues(rfidi)
#  vobsall <- getValues(rfobsall)
#
#  # =1 obs-void; =0 obs_dense
#  rfnoidi <- 1 - rfidi
#
#  # wet gridpoints in observed regions
#  rfidiwet <- rfidi
#  rfidiwet[] <- 0
#  rfidiwet[which(vobsall >= y_env$rain & vidi >= 0.1)] <- 1
#
#  # dry gridpoints in observed regions
#  rfididry <- rfobs
#  rfididry[] <- 0
#  rfididry[which(vobsall < y_env$rain & vidi >= 0.1)] <- 1
  
  # ---~--------------
  # define constants
  # number of grid points on the dyadic grid
  env$m_dim <- length( getValues(rdyad))
  sqrt_m_dim <- sqrt( env$m_dim)
  # max number of wavelet coefficients is the same as the number of grid points 
  env$n_dim <- env$m_dim
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))
  cat( paste( "number of observations (sparse) >", length(ixb), "\n"))
  cat( paste( "number of observations (dense) >",  length(ixa), "\n"))

  # ---~--------------
  # -- Observations --

  cat("Transform observations and compute error covariance matrix ")
 
  # multires representations
  mrobs <- list()
  mrobsb <- list()
  mrnoidi <- list()
  mridi <- list()
  mridia <- list()
  mridib <- list()
  mridic <- list()
  mridiwet <- list()
  mrididry <- list()
  mrnorain<- list()
  mrw<- list()
  mrwobs<- list()
  mrwbkg<- list()
  mrwwet<- list()
  mrwdry<- list()
  mrwa<- list()
  mrwb<- list()
  mrwc<- list()
  # aux
  rf1<-rdyad
  rf1[]<-1
  lambda <- 1
  env$wf<-"haar"
  
  # unlist the result and compute squared energies
  nn <- 0
  for (l in 1:n) {
    aux1 <- dwt.2d( as.matrix(rf1), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]
    mridib[[l]] <- dwt.2d( as.matrix(rfidib), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridib[[l]])[1] > 2)
      mridib[[l]] <- gauss2dsmooth( x=mridib[[l]], lambda=lambda, nx=dim(mridib[[l]])[1], ny=dim(mridib[[l]])[2])
    if ( any( mridib[[l]] > 1)) mridib[[l]][mridib[[l]]>1] <- 1
    if ( any( mridib[[l]] < 0)) mridib[[l]][mridib[[l]]<0] <- 0

    mridic[[l]] <- dwt.2d( as.matrix(rfidic), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridic[[l]])[1] > 2)
      mridic[[l]] <- gauss2dsmooth( x=mridic[[l]], lambda=lambda, nx=dim(mridic[[l]])[1], ny=dim(mridic[[l]])[2])
    if ( any( mridic[[l]] > 1)) mridic[[l]][ mridic[[l]]>1] <- 1
    if ( any( mridic[[l]] < 0)) mridic[[l]][ mridic[[l]]<0] <- 0

    if ( any( mridib[[l]] > mridic[[l]])) nn <- l
  }


  lXb_dyad <- list()
  for (e in 1:env$k_dim) {
    i <- fg_env$ixs[e]
    # interpolate onto the dyadic grid
    # background at grid  points
    rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
    rfxb[rfxb<y_env$rain] <- 0
#    Xb_dyad_original[,e]  <- getValues(rfxb)
#    Xb_dyad[,e]  <- getValues(rfxb)
    lXb_dyad[[e]] <- list()
    lXb_dyad[[e]]$mr[[1]] <- getValues(rfxb)
    for (j in 1:nn) {
      cat(paste("e j",e,j,"\n"))
      dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=j, boundary=env$boundary)
      lXb_dyad[[e]]$mr[[j+1]] <- dwt[[3*j+1]]
    }
  }

save(file="tmp.rdata",rfxb,env,lXb_dyad)
q()


#  lambda <- c( seq( 2**n/2, 1, length=nn), 1)
#  lambda <- c( seq( 2**nn/2, 1, length=nn), 1)
  
argv$wise_preproc_cellsmooth <- 10
  lambda <- c( (nn+1)/((nn+1)-1) * (argv$wise_preproc_cellsmooth-1) * 1/1:(nn+1) + ((nn+1)-argv$wise_preproc_cellsmooth)/((nn+1)-1), 1)
  print(lambda)

  for (l in 1:(nn+1)) {
    aux1 <- dwt.2d( as.matrix(rf1), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]
    mrobsb[[l]] <- dwt.2d( as.matrix(rfobsb), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]
##    if ( any( mrobsb[[l]] < mrnorain[[l]])) mrobsb[[l]][mrobsb[[l]]<mrnorain[[l]]] <- 0
    if ( any( mrobsb[[l]] < 0)) mrobsb[[l]][mrobsb[[l]]<0] <- 0

#    mridiwet[[l]] <- dwt.2d( as.matrix(rfidiwet), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
#    if ( dim(mridiwet[[l]])[1] > 2)
#      mridiwet[[l]] <- gauss2dsmooth(x=mridiwet[[l]],lambda=lambda[l],nx=dim(mridiwet[[l]])[1],ny=dim(mridiwet[[l]])[2])
#    if ( any( mridiwet[[l]] > 1)) mridiwet[[l]][mridiwet[[l]]>1] <- 1
#    if ( any( mridiwet[[l]] < 0)) mridiwet[[l]][mridiwet[[l]]<0] <- 0
#    mrididry[[l]] <- dwt.2d( as.matrix(rfididry), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
#    if ( dim(mrididry[[l]])[1] > 2)
#      mrididry[[l]] <- gauss2dsmooth(x=mrididry[[l]],lambda=lambda[l],nx=dim(mrididry[[l]])[1],ny=dim(mrididry[[l]])[2])
#    if ( any( mrididry[[l]] > 1)) mrididry[[l]][mrididry[[l]]>1] <- 1
#    if ( any( mrididry[[l]] < 0)) mrididry[[l]][mrididry[[l]]<0] <- 0
##    res <- set_weights(cbind(as.vector(mridiwet[[l]]),as.vector(mrididry[[l]])))
##    mrwwet[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
##    mrwdry[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    res <- set_weights(cbind(as.vector(mrididry[[l]]),as.vector(mridiwet[[l]])))
#    mrwwet[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    mrwdry[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    if ( any( mrwdry[[l]] > mrwwet[[l]])) mrobs[[l]][mrwdry[[l]] > mrwwet[[l]]] <- 0

    mridia[[l]] <- dwt.2d( as.matrix(rfidia), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridia[[l]])[1] > 2)
      mridia[[l]] <- gauss2dsmooth( x=mridia[[l]], lambda=lambda[l+1], nx=dim(mridia[[l]])[1], ny=dim(mridia[[l]])[2])
    if ( any( mridia[[l]] > 1)) mridia[[l]][mridia[[l]]>1] <- 1
    if ( any( mridia[[l]] < 0)) mridia[[l]][mridia[[l]]<0] <- 0

    mridib[[l]] <- dwt.2d( as.matrix(rfidib), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridib[[l]])[1] > 2)
      mridib[[l]] <- gauss2dsmooth( x=mridib[[l]], lambda=lambda[l+1], nx=dim(mridib[[l]])[1], ny=dim(mridib[[l]])[2])
    if ( any( mridib[[l]] > 1)) mridib[[l]][mridib[[l]]>1] <- 1
    if ( any( mridib[[l]] < 0)) mridib[[l]][mridib[[l]]<0] <- 0

    mridic[[l]] <- dwt.2d( as.matrix(rfidic), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridic[[l]])[1] > 2)
      mridic[[l]] <- gauss2dsmooth( x=mridic[[l]], lambda=lambda[l+1], nx=dim(mridic[[l]])[1], ny=dim(mridic[[l]])[2])
    if ( any( mridic[[l]] > 1)) mridic[[l]][ mridic[[l]]>1] <- 1
    if ( any( mridic[[l]] < 0)) mridic[[l]][ mridic[[l]]<0] <- 0

    res <- set_weights( cbind( as.vector(mridic[[l]]), as.vector(mridia[[l]]), as.vector(mridib[[l]])))
#    res <- set_weights( cbind( as.vector(mridic[[l]]), as.vector(mridia[[l]]), as.vector(mridib[[l]])))
#    mrwobs[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    mrwbkg[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
    mrwa[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
    mrwb[[l]] <- array(data=as.matrix(res[,3]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
    mrwc[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    res <- set_weights(cbind(as.vector(mrnoidi[[l]]),as.vector(mridi[[l]])))
#    mrwobs[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    mrwbkg[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
  }

  mridia0 <- gauss2dsmooth( x=as.matrix(rfidia), lambda=lambda[1], nx=sqrt_m_dim, ny=sqrt_m_dim)
  if ( any( mridia0 > 1)) mridia0[mridia0>1] <- 1
  if ( any( mridia0 < 0)) mridia0[mridia0<0] <- 0
  mridib0 <- gauss2dsmooth( x=as.matrix(rfidib), lambda=lambda[1], nx=sqrt_m_dim, ny=sqrt_m_dim)
  if ( any( mridib0 > 1)) mridib0[mridib0>1] <- 1
  if ( any( mridib0 < 0)) mridib0[mridib0<0] <- 0
  mridic0 <- gauss2dsmooth( x=as.matrix(rfidic), lambda=lambda[1], nx=sqrt_m_dim, ny=sqrt_m_dim)
  if ( any( mridic0 > 1)) mridic0[mridic0>1] <- 1
  if ( any( mridic0 < 0)) mridic0[mridic0<0] <- 0
  res <- set_weights( cbind( as.vector(mridic0), as.vector(mridia0), as.vector(mridib0)))
  mrwa0 <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim,sqrt_m_dim))
  mrwb0 <- array(data=as.matrix(res[,3]),dim=c(sqrt_m_dim,sqrt_m_dim))
  mrwc0 <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim,sqrt_m_dim))

 
save(file="tmp.rdata",lambda,mridi,mrnoidi,mridiwet,mrididry,mrobs,mrnorain,rdyad,mrw,mrwobs,mrwbkg,mrwwet,mrwdry,env,mrwa,mrwb,mrwc,mridia,mridib,mridic,mrobsb,mridia0,mridib0,mridic0,mrwa0,mrwb0,mrwc0)
print(nn)
#q()
#i<-1;r<-aggregate(rdyad,fact=2**i);r[]<-mrobs[[i]];image(r,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#i<-1;r<-aggregate(rdyad,fact=2**i);r[]<-mrwobs[[i]];contour(r,levels=c(0,0.4,0.8),add=T)


  if (nn == 0) return(NULL)
  env$nn_dim <- max( ijFromLev( n, (nn+1), F))
  cat("\n")

  #--------------------------------------------------------    
  # Pre-processing
  cat("Pre-processing on the transformed space \n")

  Xb_dyad <- array( data=NA, dim=c( env$m_dim, env$k_dim))
  Xb_dyad_original <- array( data=NA, dim=c( env$m_dim, env$k_dim))

  #--------------------------------------------------------    
  # MAIN LOOP
  for (loop in 1:max_it) {

    t0 <- Sys.time()
    cat( paste(" loop", loop, " "))

    #--------------------------------------------------------    
    # set the background

    Ub     <- array( data=NA, dim=c( env$nn_dim, env$k_dim))
    Vb     <- array( data=NA, dim=c( env$nn_dim, env$k_dim))
    vo_Vb  <- array( data=NA, dim=c( env$nn_dim, env$k_dim))

    if (plot & loop == 1) fxb <- vector()
    
    for (e in 1:env$k_dim) {

      # first iteration - background from ensemble 
      if ( loop == 1) {
        i <- fg_env$ixs[e]
        # interpolate onto the dyadic grid
        # background at grid  points
        rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
        rfxb[rfxb<y_env$rain] <- 0
        Xb_dyad_original[,e]  <- getValues(rfxb)
        Xb_dyad[,e]  <- getValues(rfxb)
        Xbpp_dyad  <- Xb_dyad
        Xbpp_dyad[] <- NA
      # second iteration onwards - background from previous iteration 
      } else {
        # -- background at grid  points --
        rfxb[] <- Xbpp_dyad[,e]
        rfxb[rfxb<y_env$rain] <- 0
        Xb_dyad[,e]  <- getValues(rfxb)
      }

      # background at observation points
      rfyb <- rdyad; rfyb[ixa] <- rfxb[ixa]
      # innovation 
      rfin <- rdyad; rfin[ixa] <- rfobs[ixa] - rfxb[ixa]
      # transformation operator
      dwtxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtyb <- dwt.2d( as.matrix(rfyb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtin <- dwt.2d( as.matrix(rfin), wf=env$wf, J=(nn+1), boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in 1:(nn+1)) {
        ij <- ijFromLev( n, l, F)
        Ub[ij[1,1]:ij[3,2],e] <- c( dwtxb[[3*l-2]], dwtxb[[3*l-1]], dwtxb[[3*l]])
        Vb[ij[1,1]:ij[3,2],e] <- c( dwtyb[[3*l-2]], dwtyb[[3*l-1]], dwtyb[[3*l]])
        vo_Vb[ij[1,1]:ij[3,2],e] <- c( dwtin[[3*l-2]], dwtin[[3*l-1]], dwtin[[3*l]])
      }
    } # end loop over ensemble members

    # set the background: end
    #--------------------------------------------------------    

    #--------------------------------------------------------    
    # Background error covariance matrices
    Ab <- Ub - rowMeans(Ub)
    Db <- Vb - rowMeans(Vb)
    Ao_b <- vo_Vb - rowMeans(vo_Vb)
    var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) 
#    var_d1_yy <- 1/(env$k_dim-1) * rowSums( vo_Vb *vo_Vb) 
    var_d1_yy <- 1/(env$k_dim-1) * rowSums( Ao_b * Ao_b) 
    # scaling coefficient used in pre-processing
    coeff <- var_b1_xy / var_d1_yy
    coeff[!is.finite(coeff)] <- 0
    rm( Ab, Db)

#    #--------------------------------------------------------    
#    # sort of root mean squared error of the decomposed innovation
#    env$costf[loop] <- mean( sqrt( var_d1_yy))
#    cat( paste("rmse", round(env$costf[loop],5)), "\n")

    #--------------------------------------------------------    
    # Rescale the different components
env$mrxa<-list()    
    Ua  <- array( data=NA, dim=c(     env$nn_dim, env$k_dim))
    for (e in 1:env$k_dim) {
      cat(".")
      Ua[,e] <- Ub[,e] + coeff * vo_Vb[,e]
      rfxb[] <- Xb_dyad[,e]
      r <- rdyad
# revise the weighting function as in set_weights
      for (l in (nn+1):1) {
        dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=l, boundary=env$boundary)
        for (i in 1:(length(dwt)-1)) dwt[[i]][] <- 0
        if (l < (nn+1)) dwt[[length(dwt)]][] <- mrxa
        ij <- ijFromLev( n, l, F)
        dwt[[3*l-2]][] <- Ua[ij[1,1]:ij[1,2],e]
        dwt[[3*l-1]][] <- Ua[ij[2,1]:ij[2,2],e]
        dwt[[3*l]][]   <- Ua[ij[3,1]:ij[3,2],e]
        xbadj <- idwt.2d( dwt) 
        if (!is.na(y_env$rain)) xbadj[xbadj<y_env$rain] <- 0 
        r[] <- array(data=as.matrix(xbadj),dim=c(sqrt_m_dim,sqrt_m_dim))
        if (l > 1) {
          mrxbadj <- dwt.2d( as.matrix(r), wf=env$wf, J=(l-1), boundary=env$boundary)[[3*(l-1)+1]]
          mrxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(l-1), boundary=env$boundary)[[3*(l-1)+1]]
##          mrxbadj <- mrwbkg[[l-1]] * mrxb + mrwobs[[l-1]] * mrxbadj
###          if (!is.na(y_env$rain)) { mrxb[mrxb<y_env$rain] <- 0; mrxbadj[mrxbadj<y_env$rain] <- 0 }
##          mrxa <- mrwbkg[[l-1]] * mrxbadj + mrwobs[[l-1]] * mrobs
          mrxa <- mrwa[[l-1]] * mrobsb[[l-1]] +  mrwb[[l-1]] * mrxbadj + mrwc[[l-1]] * mrxb
env$mrxa[[l-1]]<-mrxa
        } else {
#          xbadj <- getValues(r)
##          xbadj <- getValues(rfidiwet)/(getValues(rfidiwet)+getValues(rfididry)) * rfobs +  getValues(rfididry)/(getValues(rfidiwet)+getValues(rfididry)) * xbadj
##          xbadj[getValues(rfididry)>getValues(rfidiwet)] <- 0
###          xbadj <- getValues(rfidiwet) * xbadj
##          if (!is.na(y_env$rain)) { xbadj[xbadj<y_env$rain] <- 0 }
##          Xbpp_dyad[,e] <- getValues(rfnoidi)/(getValues(rfnoidi)+getValues(rfidi)) * Xb_dyad[,e] + getValues(rfidi)/(getValues(rfnoidi)+getValues(rfidi)) * xbadj
          Xbpp_dyad[,e] <-  mrwa0 * getValues(rfobsb) + mrwb0 * getValues(r) + mrwc0 * Xb_dyad[,e]
          if (!is.na(y_env$rain)) Xbpp_dyad[,e][Xbpp_dyad[,e]<y_env$rain] <- 0 
        }
      }  # end loop over scales
    }  # end loop over ensembles
save(file="tmp1.rdata",env,Xb_dyad,Xbpp_dyad,y_env,Xb_dyad_original,rfobsb,r)

    #
    #--------------------------------------------------------    
    # strict condition for convergence
    r[] <- rowMeans(Xb_dyad)
    vxb <- extract( r, cbind( y_env$yo$x, y_env$yo$y))
    r[] <- rowMeans(Xbpp_dyad)
    vxbpp <- extract( r, cbind( y_env$yo$x, y_env$yo$y))

    ets_thr <- y_env$rain
    if (any(y_env$yo$value>y_env$rain)) ets_thr <- median(y_env$yo$value[y_env$yo$value>y_env$rain])

    ets_b <- wise_ets( vxb, y_env$yo$value, ets_thr)
    ets_bpp <- wise_ets( vxbpp, y_env$yo$value, ets_thr)
    env$costf[loop] <- ets_bpp - ets_b
    cat(paste("\n","ets before after Delta% =",round(ets_b,3),round(ets_bpp,3),
                                             round(100*(ets_b-ets_bpp)/ets_b,1),"\n"))
    # break out of the main loop 
    if ( env$costf[loop] < 0) {
      cat( paste( "break out of the loop and keep the adjustments of loop", (loop-1), "\n"))
      Xbpp_dyad <- Xb_dyad
      break
    }
    t1 <- Sys.time()
#    if (loop==20) {
#      save(file="tmp.rdata",coeff,env,mrobs,mrnoidi,mridi,mridiwet,mrididry,mrxa,n,nn,rdyad,Xb_dyad,Ua,Ub,vo_Vb,y_env,mrxb,mrxbadj,Xb_dyad_original,rfobs,rfnoidi,rfidi,rfidiwet,rfididry)
#i<-1; s<-rdyad; s[]<-Xbpp_dyad[,i]; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#i<-1; s<-rdyad; s[]<-Xb_dyad_original[,i]; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#points(y_env$yo$x,y_env$yo$y,cex=0.25,pch=21,bg="beige",col="beige")
# i<-1; s[]<-Xb_dyad_original[,i]; yb<-extract(s,cbind(y_env$yo$x,y_env$yo$y)); plot(y_env$yo$value,yb)
# i<-1; s[]<-Xbpp_dyad[,i]; ya<-extract(s,cbind(y_env$yo$x,y_env$yo$y)); plot(y_env$yo$value,ya)
#s<-env$rfobs; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
# 
#      q()
#    }
  } # end main loop

  # save the adjusted background fields in "env"
  env$Xbpp_dyad <- Xbpp_dyad

} # END FUNCTION

wise_ets <- function(pred,ref,thr) {
  flag <- !is.na(pred) & !is.na(ref)
  hits <- length( which( flag & ref >= thr & pred >= thr))
  fals <- length( which( flag & ref <  thr & pred >= thr))
  miss <- length( which( flag & ref >= thr & pred <  thr))
  corn <- length( which( flag & ref <  thr & pred <  thr))
  hits_random <- (hits+miss) * (hits+fals) / (hits+fals+miss+corn)
  den <- (hits+fals+miss-hits_random)
  if (den == 0) return(NA)
  return( (hits-hits_random)/den )
}

#+ set IDI-based weights for the blending of multiple fields
set_weights <- function( idis, 
                         small=0.0000001) {
#------------------------------------------------------------------------------
# We have n fields to blend toghether in one field. Different fields are the
# results of spatial analysis with the same statistical interpolation (SI)
# method, applied with different settings.
# SI produces a field where each value is an average over a spatial support, 
# and the sizes of those spatial supports vary across the domain.
# Some settings are such that -on average- the spatial support is
# smaller than for others. 
# For each point in each field, an indicator -named IDI (0-1)- is available 
# that summarizes the amount of information derived from nearby observations
# which is used in the SI:
#  - 0= no info from nearby obs= spatial support large, compared to SI setup
#  - 1= a lot of info from nearby obs= spatial support comparable to SI setup
#
# The blending is done such that more weith is given to the field with smaller
# spatial support (higher resolution) and weight has been given to the other 
# (coarser) fields only if they provides additional information.
# The relative content of information in a field wrt another field is 
# computed through a function inspired by the Shannon's measure of information
# content (Tarantola, sec. 1.2.5).
# The fields are ordered from the one with smaller spatial supports to the one
# with larger spatial supports. Then, relative content of information in one
# field wrt its immediate neighbour with smaller spatial supports is computed.
# The first field, the one with smaller spatial support, is assigned the 
# completmentary of the relative content of information in the second field.
#
# The weights are the normalized relative contents of information.
#
# - Inputs -
# idis, matrix (npoints x nscales) IDI, finer scale [,1] to coarser scale [,n]
# small, constant
# - Output -
# w, matrix (npoints x nscales) weights for the blending
#------------------------------------------------------------------------------
  npoints <- dim(idis)[1]
  nscales <- dim(idis)[2]
  if ( nscales < 2) return(NULL)
  # initialization
  idis[idis<=0] <- small
  idis[idis>=1] <- 1-small
  #
  w <- array( data=NA, dim=dim(idis))
  # relative content of information in idi(i) wrt idi(i-1), coarser wrt finer
  for (i in 2:nscales) 
    w[,i] <- idis[,i] * log( idis[,i] / idis[,(i-1)])
  # adjust small things
  w[w<small] <- 0
  w[w>(1-small)] <- 1
  # weight for the finer resolution
  w[,1] <- 1 - w[,2]
  # normalize weights to 1
  w <- w / rowSums( w)
  w
}
#w <- set_weights( cbind( getValues(i200), getValues(i400), 
#                         getValues(i600), getValues(i800)))
#r <- w[,1] * r200 + w[,2] * r400 + w[,3] * r600 + w[,4] * r800


