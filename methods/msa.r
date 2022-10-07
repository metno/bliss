#+ Multiscale Alignment of gridded fields (pre-proc step of Ensemble Statistical Interpolation)
msa <- function( argv, y_env, fg_env, env, use_fg_env=T) {
# MSA Multiscale Alignment Method for Spatial analysis with Displacement Errors
#
# Multi-resolution tree-structured data models
# mrtree 
# mrobs 
# mrbkg
#  data[[j]]$Eor background ensemble (original)
#  data[[j]]$E background ensemble (updated)
#  data[[j]]$HE ensemble at observation locations (updated)

# argv$msa_ididense
# argv$msa_eps2
# argv$analysis_range
# argv$pmax
# argv$corrfun

#
# Inspired to:
# Ying, Y. (2019). A Multiscale Alignment Method for Ensemble Filtering with 
#  Displacement Errors, Monthly Weather Review, 147(12), 4553-4565.
#------------------------------------------------------------------------------
#  library(smoothie)

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

  # multi-resolution background, tree-structured data model
  mrbkg <- list()

  for (j in 1:jstop) {
    mrbkg$data[[j]] <- list()
    mrbkg$data[[j]]$E <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
    mrbkg$data[[j]]$Eor <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
    mrbkg$data[[j]]$HE <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
    # loop over ensemble members
    for (e in 1:env$k_dim) {
      if ( j == 1) {
        r <- subset( fg_env$fg[[fg_env$ixf[fg_env$ixs[e]]]]$r_main, subset=fg_env$ixe[fg_env$ixs[e]])
      } else {
        raux    <- mrtree$raster[[j-1]]$r
        raux[]  <- mrbkg$data[[j-1]]$E[,e]
        r       <- aggregate( raux, fact=2, fun=mean, expand=T, na.rm=T)
      }
      mrbkg$data[[j]]$HE[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
      mrbkg$data[[j]]$E[,e] <- getValues(r)
      mrbkg$data[[j]]$Eor[,e] <- getValues(r)
    }
    if (!is.na(y_env$rain)) mrbkg$data[[j]]$E[mrbkg$data[[j]]$E<y_env$rain] <- 0
    if (!is.na(y_env$rain)) mrbkg$data[[j]]$HE[mrbkg$data[[j]]$HE<y_env$rain] <- 0
  } # end loop over spatial scales to get the multi-resolution background

  # Initializations
  envtmp$x <- mrtree$x[[1]]
  envtmp$y <- mrtree$y[[1]]
  envtmp$m_dim <- mrtree$m_dim[[1]]
  ra <- mrtree$raster[[1]]$r
  rb <- mrtree$raster[[1]]$r

  # Loop over spatial scales
  for (j in jstop:1) {
    cat( paste( "alligning spatial level", j))
    t0b <- Sys.time()

    # Spatial anal
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
    envtmp$D <- envtmp$obs_val - envtmp$HE
    envtmp$eps2 <- rep( argv$msa_eps2, envtmp$m_dim) 
    envtmp$nn2 <- nn2( cbind(mrobs$x[[j]],mrobs$y[[j]]), 
                       query = cbind(mrtree$x[[j]],mrtree$y[[j]]), 
                       k = min( c(argv$pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( enoi_basic_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corrfun,
                                           dh=mrtree$mean_res[[j]],
                                           idi=F)))
    # no-multicores
    } else {
      res <- t( mapply( enoi_basic_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corrfun, 
                                         dh=mrtree$mean_res[[j]],
                                         idi=F)))
    }
    Ea <- res[,1:env$k_dim]
    if (!is.na(y_env$rain)) Ea[Ea<y_env$rain] <- 0
    if (any(is.na(Ea))) Ea[is.na(Ea)] <- 0

    # Loop over ensembles
    for (e in 1:env$k_dim) {
      ra[] <- Ea[,e]
      rb[] <- mrbkg$data[[j]]$E[,e]
      of_nlevel <- min( c( floor(log2(nrow(rb))), floor(log2(ncol(rb)))))
      of <- optical_flow_HS( rb, ra, nlevel=of_nlevel, niter=100, w1=100, w2=0, tol=0.00001)
#        print( paste( "j e of_par",j,e,of_par$par[1],of_par$par[2],of_par$par[3],
#                      round(range(getValues(of$u))[1],0),round(range(getValues(of$u))[2],0),
#                      round(range(getValues(of$v))[1],0),round(range(getValues(of$v))[2],0)))
      # Align smaller scales
      if (j>1) {
        of_u <- of$u
        of_v <- of$v
        # Loop over scales
        for (jj in (j-1):1) {
          rb_finer <- mrtree$raster[[jj]]$r
          rb_finer[] <- mrbkg$data[[jj]]$E[,e]
          of_u <- resample( of_u, mrtree$raster[[jj]]$r, method="bilinear")
          of_v <- resample( of_v, mrtree$raster[[jj]]$r, method="bilinear")
          rbmod_finer <- warp( rb_finer, -of_u, -of_v, method="bilinear")
          rbmod_finer[is.na(rbmod_finer)] <- 0
          if (!is.na(y_env$rain)) rbmod_finer[rbmod_finer<y_env$rain] <- 0
          mrbkg$data[[jj]]$E[,e] <- getValues( rbmod_finer)
          mrbkg$data[[jj]]$HE[,e] <- extract( rbmod_finer, cbind( mrobs$x[[jj]], mrobs$y[[jj]]))
        }  # END - Loop over scales 
        rm( rb_finer, of_u, of_v, rbmod_finer)
      # Align the smallest scale
      } else {
        rbmod <- warp( rb, -of$u, -of$v, method="bilinear")
        rbmod[is.na(rbmod)] <- 0
        if (!is.na(y_env$rain)) rbmod[rbmod<y_env$rain] <- 0
        mrbkg$data[[1]]$E[,e] <- getValues( rbmod)
        mrbkg$data[[1]]$HE[,e] <- extract( rbmod, cbind( mrobs$x[[1]], mrobs$y[[1]]))
        rm(rbmod)
      } # END IF - Align smaller scales 
    } # END - Loop over ensembles
    t1b <- Sys.time()
    cat( paste( "total time", round(t1b-t0b,1), attr(t1b-t0b,"unit"), "\n"))
  } # END - Loop over spatial scales
  if ( exists( "ra" ))     rm(ra)
  if ( exists( "rb" ))     rm(rb)
  if ( exists( "mrobs" ))  rm(mrobs)
  if ( exists( "mrtree" ))  rm(mrtree)
  if ( exists( "Ea" ))     rm(Ea)

  # save gridded analysis data for output 
  # initialization
  env$Xa <- array( data=NA, dim=c( ncell(env$rmaster), env$k_dim))
  r <- env$rmaster
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$k_dim))
  if (env$cv_mode | env$cv_mode_random) {
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$k_dim))
  }
  # loop over ensemble members
  for (e in 1:env$k_dim) {
    env$Xa[,e] <- mrbkg$data[[1]]$E[,e]
    # extract point values (crossvalidation)
    r[] <- env$Xa[,e]
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a[,e] <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y), method="simple")
  } # END loop over ensemble members

} # END FUNCTION

#match_rasters <- function(a,nlevel) { 
##  of <- optical_flow_HS( rb, ra, nlevel=nlevel, niter=a[1], w1=a[2], w2=a[3])
##print(a)
##  of <- optical_flow_HS( rb, ra, nlevel=nlevel, niter=a[1], w1=a[2], w2=0,tol=0.00001)
#  of <- optical_flow_HS( rb, ra, nlevel=nlevel, niter=100, w1=a[1], w2=0,tol=0.00001)
#  rbmod <- warp( rb, -of$u, -of$v)
#  if ( any( is.na( getValues(rbmod)))) rbmod[is.na(rbmod)] <- 0
#  costf <- sqrt(mean( getValues(ra-rbmod)**2))
#  costf
#}

###        if (e == 1 & j==jstop) of_par <- optim( c(100,0.25,0.03), match_rasters, method="Nelder-Mead", control=list( maxit=100, parscale=c(1,.1,.1)), nlevel=of_nlevel)
###        of <- optical_flow_HS( rb, ra, nlevel=of_nlevel, niter=of_par$par[1], w1=of_par$par[2], w2=of_par$par[3])
#        if (e == 1) of_par <- optim( 100, match_rasters, method="BFGS", control=list( maxit=100, parscale=c(10000)), nlevel=of_nlevel)
#        of <- optical_flow_HS( rb, ra, nlevel=of_nlevel, niter=100, w1=of_par$par[1], w2=0, tol=0.00001)

#image(rb,breaks=c(0,0.1,1,2,4,8,16,32),col=c("gray",rev(rainbow(6))));plot_arrows(of$u, of$v, fact=1, length=0.03)
#dev.new()
#image(ra,breaks=c(0,0.1,1,2,4,8,16,32),col=c("gray",rev(rainbow(6))))
#dev.new()
#image(rbmod,breaks=c(0,0.1,1,2,4,8,16,32),col=c("gray",rev(rainbow(6))))
#or<-F;e<-1;j<-1;r<-mrtree$raster[[j]]$r;if(or){r[]<-mrbkg$data[[j]]$Eor[,e]}else{r[]<-mrbkg$data[[j]]$E[,e]};image(r,breaks=c(0,0.1,1,2,4,8,16,32),col=c("gray",rev(rainbow(6))));ix0<-which(mrobs$val[[j]]<0.1);points(mrobs$x[[j]][ix0],mrobs$y[[j]][ix0],pch=21,col="gold",bg="gold",cex=0.2);ix1<-which(mrobs$val[[j]]>=0.1);points(mrobs$x[[j]][ix1],mrobs$y[[j]][ix1],pch=21,col="black",bg="black",cex=0.2);
#e<-3;j<-1;r<-mrtree$raster[[j]]$r;r[]<-Ea[,e];image(r,breaks=c(0,0.1,1,2,4,8,16,32),col=c("gray",rev(rainbow(6))))
#j<-2;r<-mrtree$raster[[j]]$r;r[mrobs$ix[[j]]]<-mrobs$val[[j]];image(r,breaks=c(0,0.1,1,2,4,8,16,32),col=c("gray",rev(rainbow(6))))
