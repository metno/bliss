#+
msaensi <- function( argv, y_env, fg_env, env) {
#------------------------------------------------------------------------------
  library(waveslim)
  library(smoothie)

  t0a <- Sys.time()

  cat( "-- MSA-EnSI --\n")

  # number of background ensemble members
  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

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
  rfxb <- rdyad
  rfyb <- rdyad
  rfin <- rdyad

  nx <- ncol( rdyad)
  ny <- nrow( rdyad)

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))


  # definition of the tree-structure, coarser spatial level has a resolution of 2**jmax cells
  jmax <- floor( log( max(nx,ny), 2))

  # --- Multi-resolution tree-structured data models --- 
  mrtree <- list()
  # multi-resolution observation, tree-structured data model
  mrobs <- list()

  # loop over the spatial levels in the tree-structured data model
  for (j in 1:jmax) {
    print( paste( "j jmax", j, jmax))
    mrtree$raster[[j]] <- list()
    # j=1 -> finest level has the same grid as the master grid (nx,ny)
    if ( j == 1) {
      # initialization from merged observations(they covers the wider region where we have observations)
      mrtree$raster[[j]]$r   <- rdyad 
      mrtree$raster[[j]]$r[] <- NA
      raux <- env$rmaster
      raux[] <- env$mergeobs$value
      robs   <- resample( raux, rdyad, method="bilinear")
      raux[] <- env$mergeobs$idi
      ridi   <- resample( raux, rdyad, method="bilinear")
      ridi[is.na(ridi)] <- 0
      ridi[ridi>1]      <- 1
    # j-th level has a grid of approximately (nx/[2**(j-1)],ny/[2**(j-1)]; approx because the coarser domain is expanded to cover the finer one)
    } else {
      raux    <- mrtree$raster[[j-1]]$r
      raux[]  <- mrobs$val_all[[j-1]]
      robs    <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      raux[]  <- mrobs$idi[[j-1]]
      ridi    <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      rm(raux)
      mrtree$raster[[j]]$r   <- ridi
      mrtree$raster[[j]]$r[] <- NA
    }

    # select only gridpoints where observations are ok (IDI is larger than a threshold)
    ix <- which( getValues(ridi) >= argv$msa_ididense & !is.na( getValues(robs)))

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
      rm(xy,ix)
    }
  } # end loop over spatial levels to get the multi-resolution observations
  if ( exists( "robs" )) rm(robs)
  if ( exists( "ridi" )) rm(ridi)

  # safe-check, exit when no ok observations found 
  if (jstop == 0) return(NULL)

  mrbkg <- list()

  dwt <- list() 
  # Loop over spatial scales
  for (j in jstop:2) {
    cat( paste( "alligning spatial level", j))
    t0b <- Sys.time()
    jw <- j-1

    mrbkg$data[[j]] <- list()
    mrbkg$data[[j]]$E <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
    mrbkg$data[[j]]$HE <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
    mrbkg$data[[j]]$X  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
    mrbkg$data[[j]]$Z   <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
    mrbkg$data[[j]]$Y   <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
    if (j == jstop) { 
      mrbkg$data[[1]] <- list()
      mrbkg$data[[1]]$Eor <- array( data=NA, dim=c( mrtree$m_dim[[1]], env$k_dim))
    }

    # loop over ensemble members - prepare the background 
    for (e in 1:env$k_dim) {
      # first iteration, background is from the original ensemble
      if ( j == jstop ) {
        r <- resample( subset( fg_env$fg[[fg_env$ixf[fg_env$ixs[e]]]]$r_main, subset=fg_env$ixe[fg_env$ixs[e]]), rdyad, method="bilinear")
        if (!is.na(y_env$rain)) r[r<y_env$rain] <- 0
        mrbkg$data[[1]]$Eor[,e] <- getValues(r) 
      # 2nd,3rd,... iterations, background is from the previous iteration
      } else {
        aux <- idwt.2d( dwt[[e]])
        r <- s <- mrtree$raster[[1]]$r 
        r[] <- array(data=as.matrix(aux),dim=c(sqrt(mrtree$m_dim[1]),sqrt(mrtree$m_dim[1])))
        s[] <- mrbkg$data[[1]]$Eor[,e]
        r <- clean_and_smooth_the_field( r, y_env$rain, s, 
                                         mrobs$val[[1]], mrobs$x[[1]], mrobs$y[[1]], 
                                         lambda=2, area_small_clumps=argv$area_small_clumps)
      }
      #discrete wavelet transformation
      dwt[[e]] <- dwt.2d( as.matrix(r), wf=argv$msaensi_wf, J=jw, boundary=argv$msaensi_boundary)
      # father wavelet as the reconstructed field at the j-th spatial level
      r   <- mrtree$raster[[j]]$r
      r[] <- dwt[[e]][[3*jw+1]] / 2**jw
      # store the background in the multiresolution tree
      mrbkg$data[[j]]$E[,e]  <- getValues(r)
      mrbkg$data[[j]]$HE[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
    } # end loop over ensemble members
    mrbkg$data[[j]]$x <- rowMeans(mrbkg$data[[j]]$E)
    # background ensemble anomalies 
    mrbkg$data[[j]]$Esd <- apply( mrbkg$data[[j]]$E, FUN=function(x){sd(x)}, MAR=1)
    for (e in 1:env$k_dim) { 
      # covariances
      mrbkg$data[[j]]$X[,e] <- 1/sqrt(env$k_dim-1) * (mrbkg$data[[j]]$E[,e] - mrbkg$data[[j]]$x)
      # correlations
      mrbkg$data[[j]]$Z[,e] <- 1/sqrt(env$k_dim-1) * (mrbkg$data[[j]]$E[,e] - mrbkg$data[[j]]$x) / mrbkg$data[[j]]$Esd
      mrbkg$data[[j]]$Z[,e][!is.finite(mrbkg$data[[j]]$Z[,e])] <- 1/sqrt(env$k_dim) 
      r   <- mrtree$raster[[j]]$r
      r[] <- mrbkg$data[[j]]$Z[,e]
      mrbkg$data[[j]]$Y[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]), method="simple")
    }

    # Spatial analysis with multiresolution background and observations
    envtmp$x <- mrtree$x[[j]]
    envtmp$y <- mrtree$y[[j]]
    envtmp$m_dim <- mrtree$m_dim[[j]]
    ra <- mrtree$raster[[j]]$r
    rb <- mrtree$raster[[j]]$r
    envtmp$obs_x <- mrobs$x[[j]]
    envtmp$obs_y <- mrobs$y[[j]]
    envtmp$k_dim <- env$k_dim
    envtmp$obs_val <- mrobs$val[[j]]
    envtmp$Eb <- mrbkg$data[[j]]$E
    envtmp$HE <- mrbkg$data[[j]]$HE
    envtmp$Y <- mrbkg$data[[j]]$Y
    envtmp$Z <- mrbkg$data[[j]]$Z
    envtmp$D <- envtmp$obs_val - envtmp$HE
    envtmp$eps2 <- rep( argv$msa_eps2, envtmp$m_dim) 
    envtmp$nn2 <- nn2( cbind(mrobs$x[[j]],mrobs$y[[j]]), 
                       query = cbind(mrtree$x[[j]],mrtree$y[[j]]), 
                       k = min( c(argv$pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( enoi_Evensen2003_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corrfun,
                                           dh=mrtree$mean_res[[j]],
                                           dh_loc=(2*mrtree$mean_res[[j]]),
                                           alpha=argv$alpha,
                                           k_dim_corr=envtmp$k_dim,
                                           idi=F)))
    # no-multicores
    } else {
      res <- t( mapply( enoi_Evensen2003_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corrfun, 
                                         dh=mrtree$mean_res[[j]],
                                         dh_loc=(2*mrtree$mean_res[[j]]),
                                         alpha=argv$alpha,
                                         k_dim_corr=envtmp$k_dim,
                                         idi=F)))
    }
    Ea <- res[,1:env$k_dim]
    if (any(is.na(Ea))) Ea[is.na(Ea)] <- 0

    # Loop over ensembles
    for (e in 1:env$k_dim) {
      ra[] <- Ea[,e]
      rb[] <- mrbkg$data[[j]]$E[,e]
      of_nlevel <- min( c( floor(log2(nrow(rb))), floor(log2(ncol(rb)))))
      if ( floor(log2(nrow(rb))) == of_nlevel & floor(log2(ncol(rb))) == of_nlevel) 
        of_nlevel <- of_nlevel - 1
      of <- optical_flow_HS( rb, ra, nlevel=of_nlevel, niter=100, w1=100, w2=0, tol=0.00001)
#        print( paste( "j e of_par",j,e,of_par$par[1],of_par$par[2],of_par$par[3],
#                      round(range(getValues(of$u))[1],0),round(range(getValues(of$u))[2],0),
#                      round(range(getValues(of$v))[1],0),round(range(getValues(of$v))[2],0)))
      # Align smaller scales
      of_u <- of$u
      of_v <- of$v
      # do not align when idi is 0
      # note: idi is calculated by mergeobs
      of_u[mrobs$idi[[j]]==0] <- 0
      of_v[mrobs$idi[[j]]==0] <- 0
      # Loop over scales
      for (jj in j:2) {
        jjw <- jj-1
#        of_u <- resample( of_u, mrtree$raster[[jj]]$r, method="bilinear")
#        of_v <- resample( of_v, mrtree$raster[[jj]]$r, method="bilinear")
        of_u <- resample( of_u, mrtree$raster[[jj]]$r, method="ngb")
        of_v <- resample( of_v, mrtree$raster[[jj]]$r, method="ngb")
        rb_finer <- mrtree$raster[[jj]]$r
        for (ww in 1:3) {
          rb_finer[] <- dwt[[e]][[3*(jjw-1)+ww]]
#          rbmod_finer <- warp( rb_finer, -of_u, -of_v, method="bilinear")
          rbmod_finer <- warp( rb_finer, -of_u, -of_v, method="simple")
          if ( any( is.na(getValues(rbmod_finer)))) 
            rbmod_finer[is.na(rbmod_finer)] <- 0
          dwt[[e]][[3*(jjw-1)+ww]][] <- as.matrix(rbmod_finer)
        }
      }  # END - Loop over scales
      r <- mrtree$raster[[j]]$r  
      r[] <- Ea[,e] * 2**jw
      dwt[[e]][[3*jw+1]][] <- as.matrix(r)
#      rm( rb_finer, of_u, of_v, rbmod_finer)
#e<-1;Eaj <- idwt.2d( dwt[[e]]);r <- mrtree$raster[[1]]$r;r[] <- array(data=as.matrix(Eaj),dim=c(sqrt(mrtree$m_dim[1]),sqrt(mrtree$m_dim[1])));r[r<0.1]<-0;image(r,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#s<-r;s[]<-mrbkg$data[[1]]$Eor[,e]; image(s,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#t<-r;t[]<-mrobs$idi[[1]]; image(t,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#png(file="test.png",width=1200,height=1200)
#image(r,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4)))) 

#save(file="tmp.rdata",r,env,dwt,mrbkg,envtmp,j,jw,mrtree,mrobs,argv,y_env)
#dev.off()
#q()
    } # END - Loop over ensembles
    t1b <- Sys.time()
    cat( paste( "total time", round(t1b-t0b,1), attr(t1b-t0b,"unit"), "\n"))
  } # END - Loop over spatial scales

  #initialization
  env$Xa <- array( data=NA, dim=c( ncell(env$rmaster), env$k_dim))
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$k_dim))
  if (env$cv_mode | env$cv_mode_random) {
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$k_dim))
  }

  # loop over ensemble members
  for (e in 1:env$k_dim) {
    aux <- idwt.2d( dwt[[e]])
    r <- s <- mrtree$raster[[1]]$r 
    r[] <- array(data=as.matrix(aux),dim=c(sqrt(mrtree$m_dim[1]),sqrt(mrtree$m_dim[1])))
    s[] <- mrbkg$data[[1]]$Eor[,e]
    r <- clean_and_smooth_the_field( r, y_env$rain, s, 
                                     mrobs$val[[1]], mrobs$x[[1]], mrobs$y[[1]], 
                                     lambda=2, area_small_clumps=argv$area_small_clumps)
    t <- resample( r, env$rmaster, method="ngb")
    env$Xa[,e] <- getValues(t)
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a[,e] <- extract( t, cbind( y_env$yov$x, y_env$yov$y), method="simple")
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( t, cbind(  y_env$yo$x, y_env$yo$y), method="simple")
  } # END loop over ensemble members

} # END FUNCTION
