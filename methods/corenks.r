#+ Change-of-Resolution Ensemble Kalman Smoother
corenks <- function( argv, y_env, fg_env, env) {
#
#------------------------------------------------------------------------------

  t0a <- Sys.time()

  cat( "-- corenks --\n")

  # number of background ensemble members
#  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

  # set domain
#  xmn <- as.numeric(argv$grid_master.x1) - as.numeric(argv$grid_master.resx)/2
#  xmx <- as.numeric(argv$grid_master.xn) + as.numeric(argv$grid_master.resx)/2
#  ymn <- as.numeric(argv$grid_master.y1) - as.numeric(argv$grid_master.resy)/2
#  ymx <- as.numeric(argv$grid_master.yn) + as.numeric(argv$grid_master.resy)/2

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
      ridi[ridi>1] <- 1
      rerrvar <- mrtree$raster[[j]]$r
      rerrvar[] <- env$mergeobs$mergedobs_errvar
    # j-th level has a grid of approximately (nx/2**j,ny/2**j)
    } else {
      raux <- mrtree$raster[[j-1]]$r
      raux[] <- mrobs$val_all[[j-1]]
#      raux[mrobs$ix[[j-1]]] <- mrobs$val[[j-1]]
      robs <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      rerrvar <- aggregate( raux, fact=2, fun=sd, expand=T, na.rm=T)**2/4
      raux[] <- mrobs$idi[[j-1]]
      ridi <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
#      raux[] <- mrobs$errvar_all[[j-1]]
#      raux[mrobs$ix[[j-1]]] <- mrobs$errvar[[j-1]]
#      rerrvar <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
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
      mrobs$errvar[[j]] <- getValues(rerrvar)[ix]
      mrobs$x[[j]] <- xy[ix,1]
      mrobs$y[[j]] <- xy[ix,2]
      mrobs$val_all[[j]] <- getValues(robs)
      mrobs$errvar_all[[j]] <- getValues(rerrvar)
#      print(paste("j d_dim",j,mrobs$d_dim[[j]]))
    }
  } # end loop over spatal levels
  # safe-check, exit when no ok observations found 
  if (jstop == 0) return(NULL)

  # -~- Fine-to-Coarse CorEnKS sweep -~-
  cat("Fine-to-Coarse CorEnKS sweep \n")
  # tree-structured data model (we identify this sweep with filtering step in EnKS)
  mrenkf <- list()
  # loop over spatial levels, from fine (j=1) to coarse (j=jstop, i.e. coarsest level with observations we can trust)
  for (j in 1:jstop) {
    cat( paste("Spatial level",j))
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
        r <- aggregate( s, fact=2, fun=mean, expand=T, na.rm=T)
      }
      # initializations (only once per level)
      if ( e == 1) {
        mrenkf$data[[j]]$E  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrenkf$data[[j]]$Xb  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrenkf$data[[j]]$Z   <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrenkf$data[[j]]$Y   <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
        mrenkf$data[[j]]$HE <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
      }
      # background ensemble members on the grid
      if (!is.na(y_env$rain)) r[r<y_env$rain] <- 0
      mrenkf$data[[j]]$E[,e] <- getValues(r)
      # background ensemble members at observation locations
      mrenkf$data[[j]]$HE[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
      if (!is.na(y_env$rain)) mrenkf$data[[j]]$HE[,e][mrenkf$data[[j]]$HE[,e]<y_env$rain] <- 0
    } # END loop over ensemble members
    # background ensemble mean on the grid
    mrenkf$data[[j]]$xb <- rowMeans(mrenkf$data[[j]]$E)
    # safe-check and selection of gridpoints having all values not NAs
    if ( (mrenkf$m_dim[[j]] <- length( ix <- which( !is.na( mrenkf$data[[j]]$xb)))) == 0) return(NULL)
    # background ensemble mean at observation loactions
#    mrenkf$data[[j]]$Hxb <- rowMeans(mrenkf$data[[j]]$HE)
    # safe-check and selection of observations having all background values not NAs
#    if ( (mrenkf$d_dim[[j]] <- length( iy <- which( !is.na( mrenkf$data[[j]]$Hxb)))) == 0) return(NULL)
    # background ensemble anomalies 
    mrenkf$data[[j]]$Esd <- apply( mrenkf$data[[j]]$E, FUN=function(x){sd(x)}, MAR=1)
    for (e in 1:env$k_dim) { 
# covariances
      mrenkf$data[[j]]$Xb[,e] <- 1/sqrt(env$k_dim-1) * (mrenkf$data[[j]]$E[,e] - mrenkf$data[[j]]$xb)
# correlations
      mrenkf$data[[j]]$Z[,e] <- 1/sqrt(env$k_dim-1) * (mrenkf$data[[j]]$E[,e] - mrenkf$data[[j]]$xb) / mrenkf$data[[j]]$Esd
      mrenkf$data[[j]]$Z[,e][!is.finite(mrenkf$data[[j]]$Z[,e])] <- 1/sqrt(env$k_dim) 
      r[] <- mrenkf$data[[j]]$Z[,e]
      mrenkf$data[[j]]$Y[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]), method="simple")
    }
    # init structure used for EnKF(EnOI)
    envtmp$x <- mrtree$x[[j]]
    envtmp$y <- mrtree$y[[j]]
    envtmp$obs_x <- mrobs$x[[j]]
    envtmp$obs_y <- mrobs$y[[j]]
    envtmp$m_dim <- mrtree$m_dim[[j]]
    envtmp$k_dim <- env$k_dim
    envtmp$obs_val <- mrobs$val[[j]]
    envtmp$E <- mrenkf$data[[j]]$E
    envtmp$HE <- mrenkf$data[[j]]$HE
    if (length(argv$corenks_eps2_range) == 2) {
      eps2_guess <- mean( mrobs$errvar[[j]], na.rm=T) / mean( rowMeans( mrenkf$data[[j]]$Xb**2))
      eps2 <- max( c( argv$corenks_eps2_range[1], min( c( eps2_guess, argv$corenks_eps2_range[2]))))
    } else {
      eps2 <- argv$corenks_eps2_range[1]
    }
    envtmp$eps2 <- rep( eps2, envtmp$m_dim)
    envtmp$Z <- mrenkf$data[[j]]$Z
    envtmp$Y <- mrenkf$data[[j]]$Y
#    print(paste("eps2 guess def:",round(eps2_guess,2),round(eps2,2)))
    envtmp$D <- envtmp$obs_val - envtmp$HE
    # helper to get the neighbours
    envtmp$nn2 <- nn2( cbind(mrobs$x[[j]],mrobs$y[[j]]), 
                       query = cbind(mrtree$x[[j]],mrtree$y[[j]]), 
                       k = min( c(argv$corenks_pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( enks_up_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corenks_corrfun,
                                           idi=T,
                                           dh=mrtree$mean_res[[j]])))
    # no-multicores
    } else {
      res <- t( mapply( enks_up_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corenks_corrfun, 
                                         idi=T,
                                         dh=mrtree$mean_res[[j]])))
    }
    mrenkf$data[[j]]$Ea  <- res[,1:env$k_dim]
    mrenkf$data[[j]]$idi <- res[,(env$k_dim+1)]
    if (!is.na(y_env$rain)) mrenkf$data[[j]]$Ea[mrenkf$data[[j]]$Ea<y_env$rain] <- 0
    # background ensemble mean on the grid
    mrenkf$data[[j]]$xa <- rowMeans(mrenkf$data[[j]]$Ea)
    mrenkf$data[[j]]$Xa <- 1/sqrt(env$k_dim-1) * (mrenkf$data[[j]]$Ea - mrenkf$data[[j]]$xa)
    cat("\n")
  } # END loop over spatial levels Fine-to-Coarse CorEnKS sweep

  # -~- Coarse-to-fine CorEnKS sweep -~-
  cat("Coarse-to-fine CorEnKS sweep \n")
  # tree-structured data model (we identify this sweep with smoothing step in EnKS)
  mrenks <- list()
  mrenks$data[[jstop]] <- list()
  mrenks$data[[jstop]]$Ea <- mrenkf$data[[jstop]]$Ea
  # loop over spatial levels, from coarser (j=jstop-1) to finer (j=1)
  for (j in (jstop-1):1) {
    cat( paste("Spatial level",j))
    mrenks$data[[j]] <- list()
    envtmp$x <- mrtree$x[[j]]
    envtmp$y <- mrtree$y[[j]]
    envtmp$obs_x <- mrtree$x[[j+1]]
    envtmp$obs_y <- mrtree$y[[j+1]]
    envtmp$m_dim <- mrtree$m_dim[[j]]
    envtmp$Ea <- mrenkf$data[[j]]$Ea
    envtmp$Xa_j  <- mrenkf$data[[j]]$Xa
    envtmp$Xa_j1 <- mrenkf$data[[j+1]]$Xa
    envtmp$E  <- mrenkf$data[[j+1]]$E
    eps2 <- argv$corenks_eps2_range[1]
    envtmp$eps2 <- rep( eps2, envtmp$m_dim)
    envtmp$D <- mrenks$data[[j+1]]$Ea - mrenkf$data[[j+1]]$E
    # helper to get the neighbours
    envtmp$nn2 <- nn2( cbind(envtmp$obs_x,envtmp$obs_y), 
                       query = cbind(envtmp$x,envtmp$y), 
                       k = min( c(argv$corenks_pmax,mrtree$m_dim[[j+1]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( enks_down_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corenks_corrfun,
                                           dh=mrtree$mean_res[[j]])))
    # no-multicores
    } else {
      res <- t( mapply( enks_down_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corenks_corrfun, 
                                         dh=mrtree$mean_res[[j]])))
    }
    mrenks$data[[j]]$Ea <- res[,1:env$k_dim]
    if (!is.na(y_env$rain)) 
      mrenks$data[[j]]$Ea[mrenks$data[[j]]$Ea<y_env$rain] <- 0
    cat("\n")
  } # END loop over spatial levels Coarse-to-fine CorEnKS sweep
save(file="tmp.rdata",mrenks,mrenkf,mrtree,mrobs,argv,env)
q()
  # save gridded analysis data for output 
  env$Xa <- mrenks$data[[1]]$Ea
  env$Xb <- mrenkf$data[[1]]$E
  env$Xidi <- mrobs$idi[[1]]
  # Safe checks
  if (!is.na(argv$corenks_range[1])) 
    env$Xa[env$Xa<argv$corenks_range[1]] <- argv$corenks_range[1]
  if (!is.na(argv$corenks_range[2])) 
    env$Xa[env$Xa>argv$corenks_range[2]] <- argv$corenks_range[2]

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
