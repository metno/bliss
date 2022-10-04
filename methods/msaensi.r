#+ Multiscale Alignment Ensemble Statistical Interpolation
msaensi <- function( argv, y_env, fg_env, env, use_fg_env=T) {
# MA Multiscale Alignment Method for Spatial analysis with Displacement Errors
#
# Multi-resolution tree-structured data models
# mrtree 
# mrobs 
# mrbkg
#  data[[j]]$Eor background ensemble (original)
#  data[[j]]$E background ensemble (updated)
#  data[[j]]$HE ensemble at observation locations (updated)

# argv$msaensi_ididense
# argv$analysis_range
# argv$pmax
# argv$corrfun

#
# Inspired to:
# Ying, Y. (2019). A Multiscale Alignment Method for Ensemble Filtering with 
#  Displacement Errors, Monthly Weather Review, 147(12), 4553-4565.
#------------------------------------------------------------------------------
library(smoothie)
  t0a <- Sys.time()

  cat( "-- MSA-EnSI --\n")

  # set domain
  nx <- ncol( env$rmaster)
  ny <- nrow( env$rmaster)

  # definition of the tree-structure, coarser spatial level has a resolution of 2**jmax cells
  jmax <- floor( log( max(nx,ny), 2))

  # --- Multi-resolution tree-structured data models --- 
  mrtree <- list()
  # multi-resolution observation, tree-structured data model
  mrobs <- list()

  # loop over the spatial levels in the tree-structured data model
  for (j in 1:jmax) {
    mrtree$raster[[j]] <- list()
    # j=1 -> finest level has the same grid as the master grid (nx,ny)
    if ( j == 1) {
      # initialization from merged observations(they covers the wider region where we have observations)
      mrtree$raster[[j]]$r   <- env$rmaster
      mrtree$raster[[j]]$r[] <- NA
      robs   <- mrtree$raster[[j]]$r
      robs[] <- env$mergeobs$value
      ridi   <- mrtree$raster[[j]]$r
      ridi[] <- env$mergeobs$idi
      ridi[is.na(ridi)] <- 0
      ridi[ridi>1]      <- 1
    # j-th level has a grid of approximately (nx/2**j,ny/2**j)
    } else {
      raux    <- mrtree$raster[[j-1]]$r
      raux[]  <- mrobs$val_all[[j-1]]
      robs    <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      raux[]  <- mrobs$idi[[j-1]]
      ridi    <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      mrtree$raster[[j]]$r   <- ridi
      mrtree$raster[[j]]$r[] <- NA
    }

    # select only gridpoints where observations are ok (IDI is larger than a threshold)
    ix <- which( getValues(ridi) >= argv$msaensi_ididense & !is.na( getValues(robs)))

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
      rm(xy)
    }
  } # end loop over spatial levels

  # safe-check, exit when no ok observations found 
  if (jstop == 0) return(NULL)

  # multi-resolution background, tree-structured data model
  mrbkg <- list()

  # Can this be made faster?
env$k_dim <- 3
  for (j in jstop:1) {
    mrbkg$data[[j]] <- list()
    mrbkg$data[[j]]$E <- array( data=NA, dim=c( mrtree$m_dim[[1]], env$k_dim))
    mrbkg$data[[j]]$Eor <- array( data=NA, dim=c( mrtree$m_dim[[1]], env$k_dim))
    mrbkg$data[[j]]$HE <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
    # loop over ensemble members
    for (e in 1:env$k_dim) {
      i <- fg_env$ixs[e]
      r <- subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i])
      Ee <- gauss2dsmooth( x=matrix(getValues(r),ncol=ny,nrow=nx), lambda=2**j, nx=nx, ny=ny)
      r[] <- t(Ee)
      mrbkg$data[[j]]$HE[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
      mrbkg$data[[j]]$E[,e] <- getValues(r)
      mrbkg$data[[j]]$Eor[,e] <- getValues(r)
    }
    if (!is.na(y_env$rain)) mrbkg$data[[j]]$E[mrbkg$data[[j]]$E<y_env$rain] <- 0
    if (!is.na(y_env$rain)) mrbkg$data[[j]]$HE[mrbkg$data[[j]]$HE<y_env$rain] <- 0
  }

  # --- MSA-EnSI ---

  # Initializations
  envtmp$x <- mrtree$x[[1]]
  envtmp$y <- mrtree$y[[1]]
  envtmp$m_dim <- mrtree$m_dim[[1]]
  ra <- mrtree$raster[[1]]$r
  rb <- mrtree$raster[[1]]$r

  # Loop over spatial scales
  for (j in jstop:1) {
    envtmp$obs_x <- mrobs$x[[j]]
    envtmp$obs_y <- mrobs$y[[j]]
    envtmp$k_dim <- env$k_dim
    envtmp$obs_val <- mrobs$val[[j]]
    envtmp$Eb <- mrbkg$data[[j]]$E
    envtmp$HE <- mrbkg$data[[j]]$HE
    envtmp$D <- envtmp$obs_val - envtmp$HE
    envtmp$eps2 <- rep( 1, envtmp$m_dim) # 1? is that ok?
    envtmp$nn2 <- nn2( cbind(mrobs$x[[j]],mrobs$y[[j]]), 
                       query = cbind(mrtree$x[[1]],mrtree$y[[1]]), 
                       k = min( c(argv$pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
# see corens_up_gridpoint_by_gridpoint
      res <- t( mcmapply( enoi_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corrfun,
                                           dh=mrtree$mean_res[[j]],
                                           idi=F)))
    # no-multicores
    } else {
      res <- t( mapply( enoi_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corrfun, 
                                         dh=mrtree$mean_res[[j]],
                                         idi=F)))
    }
    Ea <- res[,1:env$k_dim]
    if (!is.na(y_env$rain)) Ea[Ea<y_env$rain] <- 0
    if (any(is.na(Ea))) Ea[is.na(Ea)] <- 0

    # Align smaller scales
    if (j>1) {
      # Loop over ensembles
      for (e in 1:env$k_dim) {
        ra[] <- Ea[,e]
        rb[] <- mrbkg$data[[j]]$E[,e]
        of <- optical_flow_HS( rb, ra, nlevel=j, niter=100, w1=0.25, w2=0.03)
        # Loop over scales
        for (jj in (j-1):1) {
          rb[] <- mrbkg$data[[jj]]$E[,e]
          rbmod <- warp( rb, -of$u, -of$v)
          if ( any( is.na( getValues(rbmod)))) rbmod[is.na(rbmod)] <- 0
#          mrbkg$data[[jj]]$E[,e] <- getValues(rbmod)
#          Ee <- gauss2dsmooth( x=matrix(getValues(rbmod),ncol=ny,nrow=nx), lambda=2, nx=nx, ny=ny)
#          r[] <- t(Ee)
          r<-rbmod
          mrbkg$data[[jj]]$HE[,e] <- extract( r, cbind( mrobs$x[[jj]], mrobs$y[[jj]]))
          mrbkg$data[[jj]]$E[,e] <- getValues(r)
          if (!is.na(y_env$rain)) mrbkg$data[[jj]]$E[mrbkg$data[[jj]]$E<y_env$rain] <- 0
          if (!is.na(y_env$rain)) mrbkg$data[[jj]]$HE[mrbkg$data[[jj]]$HE<y_env$rain] <- 0
        } # END - Loop over scales
      } # END - Loop over ensembles
    } # END - Align smaller scales
    save(file=paste0("tmp_",j,".rdata"),mrobs,mrtree,mrbkg,env,Ea)
    print(paste0("tmp_",j,".rdata"))
  } # END - Loop over spatial scales

  save(file="tmp.rdata",mrobs,mrtree,mrbkg,env)
#save(file="tmp.rdata",mrdw,mrup,mrtree,mrobs,argv,env)
#all<-T; j<-1; r<-mrtree$raster[[j]]$r; r[]<-NA; if (all) { r[] <- mrobs$val_all[[j]] } else { r[mrobs$ix[[j]]] <- mrobs$val[[j]] }; image(r,breaks=c(0,seq(0.1,45,length=10)),col=c("gray",rev(rainbow(9))))
#a<-T; e<-12; j<-1; r<-mrtree$raster[[j]]$r; if (a) { r[]<-mrdw$data[[j]]$Ea[,e] } else { r[]<-mrup$data[[j]]$E[,e] }; image(r,breaks=c(0,seq(0.1,45,length=10)),col=c("gray",rev(rainbow(9))))
#q()

  # save gridded analysis data for output 
  env$Xa <- Ea
  env$Xb <- mrbkg$data[[1]]$Eor
  env$Xidi <- mrobs$idi[[1]]
  # Safe checks
  if (!is.na(argv$analysis_range[1])) 
    env$Xa[env$Xa<argv$analysis_range[1]] <- argv$analysis_range[1]
  if (!is.na(argv$analysis_range[2])) 
    env$Xa[env$Xa>argv$analysis_range[2]] <- argv$analysis_range[2]

  # save IDI at observation locations
  r <- env$rmaster
  r[] <- env$Xidi 
  y_env$yo$mergedobs_idi  <- extract( r, cbind( y_env$yo$x, y_env$yo$y), method="simple")
  if (env$cv_mode | env$cv_mode_random)
    y_env$yov$mergedobs_idi <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
  
  # save analysis at observation locations
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$k_dim))
  if (env$cv_mode | env$cv_mode_random) {
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$k_dim))
  }
  # loop over ensemble members
  for (e in 1:env$k_dim) {
    # extract point values (crossvalidation)
    r[] <- env$Xa[,e]
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a[,e] <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y), method="simple")
  } # END loop over ensemble members

} # END FUNCTION
