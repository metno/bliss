#+ Change-of-Resolution Ensemble Rauch-Tung-Striebel Smoother
msa <- function( argv, y_env, fg_env, env, use_fg_env=T) {
#
#------------------------------------------------------------------------------
library(smoothie)
  t0a <- Sys.time()

  cat( "-- msa --\n")

  # set domain
  nx <- ncol( env$rmaster)
  ny <- nrow( env$rmaster)

  # definition of the tree-structure, coarser spatial level has a resolution of 2**jmax cells
  jmax <- floor( log( max(nx,ny), 2))

  # --- Multi-resolution tree-structured data models --- 
  mrtree <- list()
  # multi-resolution observation, tree-structured data model
  mrobs <- list()
  # background
  mrbkg <- list()
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
    ix <- which( getValues(ridi) >= argv$corens_ididense & 
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

env$k_dim <- 3
  for (j in jstop:1) {
    mrbkg$data[[j]] <- list()
    mrbkg$data[[j]]$E <- array( data=NA, dim=c( mrtree$m_dim[[1]], env$k_dim))
    mrbkg$data[[j]]$Eor <- array( data=NA, dim=c( mrtree$m_dim[[1]], env$k_dim))
    mrbkg$data[[j]]$HE <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
    # loop over ensemble members
    for (e in 1:env$k_dim) {
print(e)
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

  envtmp$x <- mrtree$x[[1]]
  envtmp$y <- mrtree$y[[1]]
  envtmp$m_dim <- mrtree$m_dim[[1]]
  ra <- mrtree$raster[[1]]$r
  rb <- mrtree$raster[[1]]$r
  for (j in jstop:1) {
    envtmp$obs_x <- mrobs$x[[j]]
    envtmp$obs_y <- mrobs$y[[j]]
    envtmp$k_dim <- env$k_dim
    envtmp$obs_val <- mrobs$val[[j]]
    envtmp$Eb <- mrbkg$data[[j]]$E
    envtmp$HE <- mrbkg$data[[j]]$HE
    envtmp$D <- envtmp$obs_val - envtmp$HE
    envtmp$eps2 <- rep( 1, envtmp$m_dim)
    envtmp$nn2 <- nn2( cbind(mrobs$x[[j]],mrobs$y[[j]]), 
                       query = cbind(mrtree$x[[1]],mrtree$y[[1]]), 
                       k = min( c(argv$corens_pmax,mrobs$d_dim[[j]])), 
                       searchtype = "radius", 
                       radius = (7*mrtree$mean_res[[j]]))
    # run EnKF/EnOI gridpoint by gridpoint
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( enoi_gridpoint_by_gridpoint,
                          1:envtmp$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          MoreArgs = list( corr=argv$corens_corrfun,
                                           dh=mrtree$mean_res[[j]],
                                           idi=F)))
    # no-multicores
    } else {
      res <- t( mapply( enoi_gridpoint_by_gridpoint,
                        1:envtmp$m_dim,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corens_corrfun, 
                                         dh=mrtree$mean_res[[j]],
                                         idi=F)))
    }
    Ea <- res[,1:env$k_dim]
    if (!is.na(y_env$rain)) Ea[Ea<y_env$rain] <- 0
    if (any(is.na(Ea))) Ea[is.na(Ea)] <- 0
    if (j>1) {
      for (e in 1:env$k_dim) {
        ra[] <- Ea[,e]
        rb[] <- mrbkg$data[[j]]$E[,e]
        of <- optical_flow_HS( rb, ra, nlevel=8, niter=20)
        for (jj in (j-1):1) {
          rb[] <- mrbkg$data[[jj]]$E[,e]
          rbmod <- warp( rb, -of$u, -of$v)
          if ( any( is.na( getValues(rbmod)))) rbmod[is.na(rbmod)] <- 0
          mrbkg$data[[jj]]$E[,e] <- getValues(rbmod)
        }
      }
    }
    save(file=paste0("tmp_",j,".rdata"),mrobs,mrtree,mrbkg,env,Ea)
    print(paste0("tmp_",j,".rdata"))
  }

  save(file="tmp.rdata",mrobs,mrtree,mrbkg,env)

}
