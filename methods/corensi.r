#+ Change-of-Resolution Ensemble Rauch-Tung-Striebel Smoother
corensi <- function( argv, y_env, fg_env, env, use_fg_env=T) {
#
#------------------------------------------------------------------------------

  t0a <- Sys.time()

  cat( "-- corensi --\n")

  # set domain
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
    # j-th level has a grid of approximately (nx/2**j,ny/2**j)
    } else {
      raux <- mrtree$raster[[j-1]]$r
      raux[] <- mrobs$val_all[[j-1]]
#      raux[mrobs$ix[[j-1]]] <- mrobs$val[[j-1]]
      robs <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      raux[] <- mrobs$idi[[j-1]]
      ridi <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      mrtree$raster[[j]]$r <- ridi
      mrtree$raster[[j]]$r[] <- NA
    }
    # select only gridpoints where observations are ok (IDI is larger than a threshold)
    ix <- which( getValues(ridi) >= argv$corensi_ididense & !is.na( getValues(robs)))

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
#      print(paste("j d_dim",j,mrobs$d_dim[[j]]))
    }
  } # end loop over spatal levels
  # safe-check, exit when no ok observations found 
  if (jstop == 0) return(NULL)

  # -~- Fine-to-Coarse corensi sweep -~-
  cat("Fine-to-Coarse corensi sweep \n")
  # tree-structured data model (we identify this sweep with filtering step in EnKS)
  mrup <- list()
  # loop over spatial levels, from fine (j=1) to coarse (j=jstop, i.e. coarsest level with observations we can trust)
  for (j in 1:jstop) {
    cat( paste("Spatial level",j))
    mrup$data[[j]] <- list()
    # loop over ensemble members
    for (e in 1:env$k_dim) {
      # finest level, get the background from the ensemble on the master grid
      if ( j == 1) {
        if (use_fg_env) {
          i <- fg_env$ixs[e]
          r <- subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i])
        } else {
          r <- env$rmaster
          r[] <- env$Xb[,e]
        }
      # jth-level, perform the saptial aggregation
      } else {
        s <- mrtree$raster[[j-1]]$r
        s[] <- mrup$data[[j-1]]$Ea[,e]
        r <- aggregate( s, fact=2, fun=mean, expand=T, na.rm=T)
      }
      # initializations (only once per level)
      if ( e == 1) {
        mrup$data[[j]]$E  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrup$data[[j]]$X  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrup$data[[j]]$Z   <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
        mrup$data[[j]]$Y   <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
        mrup$data[[j]]$HE <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
      }
      # background ensemble members on the grid
      if (!is.na(y_env$rain)) r[r<y_env$rain] <- 0
      mrup$data[[j]]$E[,e] <- getValues(r)
      # background ensemble members at observation locations
      mrup$data[[j]]$HE[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
      if (!is.na(y_env$rain)) mrup$data[[j]]$HE[,e][mrup$data[[j]]$HE[,e]<y_env$rain] <- 0
    } # END loop over ensemble members
    # background ensemble mean on the grid
    mrup$data[[j]]$x <- rowMeans(mrup$data[[j]]$E)
    # safe-check and selection of gridpoints having all values not NAs
    if ( (mrup$m_dim[[j]] <- length( ix <- which( !is.na( mrup$data[[j]]$x)))) == 0) return(NULL)
    # background ensemble mean at observation loactions
#    mrup$data[[j]]$Hx <- rowMeans(mrup$data[[j]]$HE)
    # safe-check and selection of observations having all background values not NAs
#    if ( (mrup$d_dim[[j]] <- length( iy <- which( !is.na( mrup$data[[j]]$Hx)))) == 0) return(NULL)
    # background ensemble anomalies 
    mrup$data[[j]]$Esd <- apply( mrup$data[[j]]$E, FUN=function(x){sd(x)}, MAR=1)
    for (e in 1:env$k_dim) { 
      # covariances
      mrup$data[[j]]$X[,e] <- 1/sqrt(env$k_dim-1) * (mrup$data[[j]]$E[,e] - mrup$data[[j]]$x)
      # correlations
      mrup$data[[j]]$Z[,e] <- 1/sqrt(env$k_dim-1) * (mrup$data[[j]]$E[,e] - mrup$data[[j]]$x) / mrup$data[[j]]$Esd
      mrup$data[[j]]$Z[,e][!is.finite(mrup$data[[j]]$Z[,e])] <- 1/sqrt(env$k_dim) 
      r[] <- mrup$data[[j]]$Z[,e]
      mrup$data[[j]]$Y[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]), method="simple")
    }
    # init structure used for EnKF(EnOI)
    envtmp$x <- mrtree$x[[j]]
    envtmp$y <- mrtree$y[[j]]
    envtmp$obs_x <- mrobs$x[[j]]
    envtmp$obs_y <- mrobs$y[[j]]
    envtmp$m_dim <- mrtree$m_dim[[j]]
    envtmp$k_dim <- env$k_dim
    envtmp$obs_val <- mrobs$val[[j]]
    envtmp$E <- mrup$data[[j]]$E
    envtmp$HE <- mrup$data[[j]]$HE
    eps2 <- argv$corensi_eps2
    envtmp$eps2 <- rep( eps2, envtmp$m_dim)
    envtmp$Z <- mrup$data[[j]]$Z
    envtmp$Y <- mrup$data[[j]]$Y
#    print(paste("eps2 guess def:",round(eps2_guess,2),round(eps2,2)))
    envtmp$D <- envtmp$obs_val - envtmp$HE
#    sig2b <- mean( apply( envtmp$HE, FUN=function(x){sd(x)}, MAR=1)**2)
#    sig2o_sig2b <- mean( (envtmp$obs_val - rowMeans(envtmp$HE))**2)
#    print( paste("j eps2",j,round(sig2o_sig2b,3),round(sig2b,3),round(sig2o_sig2b/sig2b-1,3)))
    # helper to get the neighbours
    envtmp$nn2 <- nn2( cbind(mrobs$x[[j]],mrobs$y[[j]]), 
                       query = cbind(mrtree$x[[j]],mrtree$y[[j]]), 
                       k = min( c(argv$pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( corensi_up_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corrfun,
                                           dh=mrtree$mean_res[[j]],
                                           alpha=argv$corensi_alpha,
                                           k_dim_corr=argv$corensi_k_dim_corr,
                                           idi=T)))
    # no-multicores
    } else {
      res <- t( mapply( corensi_up_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corrfun, 
                                         dh=mrtree$mean_res[[j]],
                                         alpha=argv$corensi_alpha,
                                         k_dim_corr=argv$corensi_k_dim_corr,
                                         idi=T)))
    }
    mrup$data[[j]]$Ea  <- res[,1:env$k_dim]
    mrup$data[[j]]$idi <- res[,(env$k_dim+1)]
    if (!is.na(y_env$rain)) mrup$data[[j]]$Ea[mrup$data[[j]]$Ea<y_env$rain] <- 0
    # background ensemble mean on the grid
    mrup$data[[j]]$xa <- rowMeans(mrup$data[[j]]$Ea)
    mrup$data[[j]]$Xa <- 1/sqrt(env$k_dim-1) * (mrup$data[[j]]$Ea - mrup$data[[j]]$xa)
    cat("\n")
  } # END loop over spatial levels Fine-to-Coarse corensi sweep

  # -~- Coarse-to-fine corensi sweep -~-
  cat("Coarse-to-fine corensi sweep \n")
  # tree-structured data model (we identify this sweep with smoothing step in EnKS)
  mrdw <- list()
  mrdw$data[[jstop]] <- list()
  mrdw$data[[jstop]]$Ea <- mrup$data[[jstop]]$Ea
  # loop over spatial levels, from coarser (j=jstop-1) to finer (j=1)
  for (j in (jstop-1):1) {
    cat( paste("Spatial level",j))
    mrdw$data[[j]] <- list()
    envtmp$x <- mrtree$x[[j]]
    envtmp$y <- mrtree$y[[j]]
    envtmp$obs_x <- mrtree$x[[j+1]]
    envtmp$obs_y <- mrtree$y[[j+1]]
    envtmp$m_dim <- mrtree$m_dim[[j]]
    envtmp$Ea <- mrup$data[[j]]$Ea
    envtmp$Xa_j  <- mrup$data[[j]]$Xa
    envtmp$X_j1 <- mrup$data[[j+1]]$X
    envtmp$E  <- mrup$data[[j+1]]$E
    eps2 <- argv$corensi_eps2
    envtmp$eps2 <- rep( eps2, envtmp$m_dim)
    envtmp$D <- mrdw$data[[j+1]]$Ea - mrup$data[[j+1]]$E
    print( paste("j qq2",j,round( 
      mean( apply( mrdw$data[[j+1]]$Ea, FUN=function(x){sd(x)}, MAR=1)**2),3)))
    # helper to get the neighbours
    envtmp$nn2 <- nn2( cbind(envtmp$obs_x,envtmp$obs_y), 
                       query = cbind(envtmp$x,envtmp$y), 
                       k = min( c(argv$pmax,mrtree$m_dim[[j+1]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( corensi_down_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T))
    # no-multicores
    } else {
      res <- t( mapply( corensi_down_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T))
    }
    mrdw$data[[j]]$Ea <- res[,1:env$k_dim]
    if (!is.na(y_env$rain)) 
      mrdw$data[[j]]$Ea[mrdw$data[[j]]$Ea<y_env$rain] <- 0
    cat("\n")
  } # END loop over spatial levels Coarse-to-fine corensi sweep
#save(file="tmp.rdata",mrdw,mrup,mrtree,mrobs,argv,env)
#all<-T; j<-1; r<-mrtree$raster[[j]]$r; r[]<-NA; if (all) { r[] <- mrobs$val_all[[j]] } else { r[mrobs$ix[[j]]] <- mrobs$val[[j]] }; image(r,breaks=c(0,seq(0.1,45,length=10)),col=c("gray",rev(rainbow(9))))
#a<-T; e<-12; j<-1; r<-mrtree$raster[[j]]$r; if (a) { r[]<-mrdw$data[[j]]$Ea[,e] } else { r[]<-mrup$data[[j]]$E[,e] }; image(r,breaks=c(0,seq(0.1,45,length=10)),col=c("gray",rev(rainbow(9))))
#q()
  # save gridded analysis data for output 
  env$Xa <- mrdw$data[[1]]$Ea
  env$Xb <- mrup$data[[1]]$E
  env$Xidi <- mrobs$idi[[1]]
  # Safe checks
  if (!is.na(argv$range[1])) 
    env$Xa[env$Xa<argv$range[1]] <- argv$range[1]
  if (!is.na(argv$range[2])) 
    env$Xa[env$Xa>argv$range[2]] <- argv$range[2]

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

id <- function(j,e,up,inno,sm_anres) {
  str<-paste( "/ lev=",j,"ens=",e)
  if (up & inno) {
    main<-paste("up, innovation",str)
    s<-mrtree$raster[[j]]$r; s[]<-NA
    s[mrobs$ix[[j]]]<-mrobs$val[[j]]-mrup$data[[j]]$HE[,e]
  } else if (up & !inno) {
    main<-paste("up, analysis residual",str)
    s<-mrtree$raster[[j]]$r; s[]<-NA
    s[]<-mrup$data[[j]]$Ea[,e]-mrup$data[[j]]$E[,e]
  } else if (!up & inno) {
    main<-paste("dw, innovation (Esm(j+1)-E(j+1))",str)
    s<-mrtree$raster[[j+1]]$r; s[]<-NA
    s[] <- mrdw$data[[j+1]]$Ea[,e] - mrup$data[[j+1]]$E[,e]
  } else if (!up & !inno & sm_anres) {
    main<-paste("dw, analysis residual",str)
    s<-mrtree$raster[[j]]$r; s[]<-NA
    s[]<-mrdw$data[[j]]$Ea[,e]-mrup$data[[j]]$Ea[,e]
  } else if (!up & !inno & !sm_anres) {
    main<-paste("dw, sm_analysis minus fg",str)
    s<-mrtree$raster[[j]]$r; s[]<-NA
    s[]<-mrdw$data[[j]]$Ea[,e]-mrup$data[[j]]$E[,e]
  }
  range <- range(getValues(s),na.rm=T)
  print(range)
  image(s,breaks=c(range[1]-0.1,range[1]/2,-0.1,0.1,range[2]/2,range[2]+0.1),col=c("blue","cornflowerblue","gray","pink","red"),main=main)
  if (!inno | (inno&!up)) {
    ix<-which(mrobs$val[[j]]>=0.1)
    points(mrobs$x[[j]][ix],mrobs$y[[j]][ix],pch=21,bg="black",col="black",cex=0.2)
    ix<-which(mrobs$val[[j]]<=0.1)
    points(mrobs$x[[j]][ix],mrobs$y[[j]][ix],pch=21,bg="white",col="white",cex=0.2)
  }
}

