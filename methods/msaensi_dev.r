#+
msaensi_dev <- function( argv, y_env, fg_env, env) {
#------------------------------------------------------------------------------
  library(waveslim)
  library(smoothie)

  t0a <- Sys.time()

  cat( "-- MSA-EnSI --\n")

  # number of background ensemble members
#  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

  # observations
  r <- env$rmaster
  
  m_dim <- ncell(r)
  xy <- xyFromCell( r, 1:ncell(r))
  envtmp$m_dim <- m_dim
  envtmp$x <- xy[,1]
  envtmp$y <- xy[,2]
#argv$mergeobs_eps2 <- 0.1
#  envtmp$eps2 <- rep( argv$mergeobs_eps2, m_dim)
#  envtmp$obs_x <- y_env$yo$x
#  envtmp$obs_y <- y_env$yo$y
  cat( paste( "         number of grid points, m dim >", m_dim, "\n"))
  cat( paste( "number of in-situ observations, p dim >", y_env$yo$n, "\n"))

  # rasterize observations
  robs_on_master <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=env$rmaster, field=y_env$yo$value, fun=mean, na.rm=T )
#  ix <- which( !is.na(getValues(raux)))
#  envtmp$obs_value <- getValues(raux)[ix]
#  envtmp$obs_x <- xy[ix,1]
#  envtmp$obs_y <- xy[ix,2]
#  envtmp$obs_n <- length(ix)
#  envtmp$eps2 <- rep( argv$mergeobs_eps2, m_dim)
#  cat( paste( "number of rasterized observations, p dim >", envtmp$obs_n, "\n"))
  
  # Loop over ensembles
  envtmp$Eb  <- array( data=NA, dim=c( ncell(env$rmaster), env$k_dim))
  ra <- env$rmaster; ra[] <- NA
  rb <- env$rmaster; rb[] <- NA

  # structure functions
  y_env$yo$value[y_env$yo$value>=0.1] <- 1
  y_env$yo$value[y_env$yo$value< 0.1] <- 0
  dh_vec <- seq(5000,500000,by=10000)
  rmse <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
  for (e in 1:env$k_dim) {
    rb <- subset( fg_env$fg[[1]]$r_main, subset=e)
    yb <- extract( rb, cbind(y_env$yo$x,y_env$yo$y))
    yb[yb>=0.1] <- 1
    yb[yb< 0.1] <- 0
#    dist <- sqrt( outer( y_env$yo$x,y_env$yo$x,FUN="-")**2 + outer( y_env$yo$y,y_env$yo$y,FUN="-")**2)
#save(file="tmp.rdata",yb,y_env,dist)
    cat( paste("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"))
    cat( paste("ensemble",e,"\n"))
    for (j in 1:length(dh_vec)) {
      cat( paste("===========================================================\n"))
      dh <- dh_vec[j]
      fact <- ceiling( (dh / 3) / mean(res(env$rmaster)))
      cat( paste("dh pmax fact",dh,fact,"\n"))
      if (fact == 1) {
        raux_agg <- env$rmaster 
      } else {
        raux_agg <- aggregate(env$rmaster, fact=fact, expand=T)
      }
      rauxo <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                         field=y_env$yo$value, fun=mean, na.rm=T )
      ix <- which( !is.na(getValues(rauxo)))
      vo <- getValues(rauxo)[ix]
      rauxb <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                          field=yb, fun=mean, na.rm=T )
      ix <- which( !is.na(getValues(rauxb)))
      vb <- getValues(rauxb)[ix]
      rmse[j,e] <- sqrt(mean((vo-vb)**2))
    }
  }
save(file="tmp.rdata",rmse,dh_vec)
 plot(dh_vec,rmse[,1],ylim=c(0,3),col="white"); for (i in 1:30) lines(dh_vec,rmse[,i])
 plot(dh_vec,rmse[,1],ylim=c(0.1,0.5),col="white"); for (i in 1:30) lines(dh_vec,rmse[,i])
points(dh_vec,rowMeans(rmse),pch=21,bg="red")
# try with spatial correlaiton instead of rmse (see Casati's euqation)
q()
#  for (e in 1:env$k_dim) {
  for (e in 1:env$k_dim) {
    first <- T
    eps2 <- 0.1
    dh_vec <- c(100000,50000,10000,5000)
    pmax_vec <- c(20,20,20,20)
#    dh_vec <- c(10000,5000)
#    pmax_vec <- c(20,20)
    cat( paste("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"))
    cat( paste("ensemble",e,"\n"))
    for (j in 1:length(dh_vec)) {
      cat( paste("===========================================================\n"))
      dh <- dh_vec[j]
      pmax <- pmax_vec[j]
      fact <- ceiling( (dh / 3) / mean(res(env$rmaster)))
      cat( paste("dh pmax fact",dh,pmax,fact,"\n"))
      if (fact == 1) {
        raux_agg <- env$rmaster 
      } else {
        raux_agg <- aggregate(env$rmaster, fact=fact, expand=T)
      }
      raux <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
#                         field=y_env$yo$value, fun=mean, na.rm=T )
                         field=y_env$yo$value, fun=max, na.rm=T )
      ix <- which( !is.na(getValues(raux)))
      xy <- xyFromCell( raux, 1:ncell(raux))
      envtmp$obs_value <- getValues(raux)[ix]
#      envtmp$obs_x <- xy[ix,1]
#      envtmp$obs_y <- xy[ix,2]
      raux_x <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                           field=y_env$yo$x, fun=mean, na.rm=T )
      envtmp$obs_x <- getValues(raux_x)[ix]
      raux_y <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                           field=y_env$yo$y, fun=mean, na.rm=T )
      envtmp$obs_y <- getValues(raux_y)[ix]
      envtmp$obs_n <- length(ix)
      cat( paste( "number of rasterized observations, p dim >", envtmp$obs_n, "\n"))
      if (first) {
        rb <- subset( fg_env$fg[[1]]$r_main, subset=e)
        first <- F
      } else {
        rb <- rbmod
        rb[is.na(rb)] <- 0
        rb[rb<0] <- 0
      }
      envtmp$xb <- getValues(rb)
      envtmp$yb <- extract( rb, cbind(envtmp$obs_x, envtmp$obs_y))
      envtmp$d  <- envtmp$obs_value - envtmp$yb
      m_step <- 100000
      xa <- getValues(ra); xa[] <- NA
      for (i in 1:ceiling(m_dim/m_step)) {
        m1 <- (i-1) * m_step + 1
        m2 <- min( c( i * m_step, m_dim))
        envtmp$m_dim <- m2 - m1 + 1
        envtmp$eps2 <- rep( eps2, envtmp$m_dim)

        xy <- xyFromCell( r, m1:m2)
        envtmp$x <- xy[,1]
        envtmp$y <- xy[,2]
        envtmp$xb <- getValues(rb)[m1:m2]

        cat( paste("-----------------------------------------------------------\n"))
        cat( paste("m1 m2 range / m_dim",m1,m2,envtmp$m_dim,"/",m_dim,"\n"))
        t00 <- Sys.time()
        envtmp$nn2 <- nn2( cbind( envtmp$obs_x, envtmp$obs_y), 
                           query = cbind( envtmp$x, envtmp$y), 
                           k = pmax, searchtype = "radius", 
                           radius = (7*dh))
        t11 <- Sys.time()
        cat( paste( "nn2 total time", round(t11-t00,1), attr(t11-t00,"unit"), "\n"))
        t00 <- Sys.time()
        # run oi gridpoint by gridpoint (idi the first time only)
        if (!is.na(argv$cores)) {
          res <- t( mcmapply( oi_basic_gridpoint_by_gridpoint,
                              1:envtmp$m_dim,
                              mc.cores=argv$cores,
                              SIMPLIFY=T,
                              MoreArgs=list( corr=argv$mergeobs_corrfun,
                                             dh=dh,
                                             idi=F,
                                             uncertainty=F)))
        # no-multicores
        } else {
          res <- t( mapply( oi_basic_gridpoint_by_gridpoint,
                            1:envtmp$m_dim,
                            SIMPLIFY=T,
                            MoreArgs=list( corr=argv$mergeobs_corrfun,
                                           dh=dh,
                                           idi=F,
                                           uncertainty=F)))
        }
        t11 <- Sys.time()
        cat( paste( "oi total time", round(t11-t00,1), attr(t11-t00,"unit"), "\n"))
        res[,1][res[,1]<0] <- 0
        xa[m1:m2] <- res[,1]
      }
      ra[] <- xa
      ra_norm <- ra; ra_norm[] <- xa / max(xa)
      rb_norm <- rb; rb_norm[] <- getValues(rb) / max(getValues(rb))

      t00 <- Sys.time()
      of_nlevel <- min( c( floor(log2(nrow(env$rmaster))), floor(log2(ncol(env$rmaster)))))
      if ( floor(log2(nrow(env$rmaster))) == of_nlevel & 
           floor(log2(ncol(env$rmaster))) == of_nlevel) 
        of_nlevel <- of_nlevel - 1
      of_nlevel <- min( 1+c( ceiling(max(c(3,log2(fact)))), of_nlevel))
#      of <- optical_flow_HS( rb_norm, ra_norm, nlevel=of_nlevel, 
      of <- optical_flow_HS( rb, ra, nlevel=of_nlevel, 
                             niter=100, w1=100, w2=0, tol=0.00001)
      t11 <- Sys.time()
      cat( paste( "optflow of_nlevel", of_nlevel, "\n"))
      cat( paste( "optflow total time", round(t11-t00,1), attr(t11-t00,"unit"), "\n"))
      u <- env$rmaster
      v <- env$rmaster
#    u[] <- t( matrix( data=as.vector(of$u[,,1]), 
#              ncol=nrow(env$rmaster), nrow=ncol(env$rmaster)))
#    v[] <- t( matrix( data=as.vector(of$v[,,1]), 
#              ncol=nrow(env$rmaster), nrow=ncol(env$rmaster)))
#    rbmod <- warp( rb, -u, -v, method="simple")
      rbmod <- warp( rb, -of$u, -of$v, method="simple")
save(file=paste0("tmp_e",formatC(e,width=2,flag="0"),"_j",j,".rdata"),argv,raux,ra,rb,ra_norm,rb_norm,res,y_env,xa,envtmp,env,of,u,v,rbmod,robs_on_master)
    }
  }
q()
source("~/projects/bliss/functions/optflow/optflow_util.r")
###
save(file=paste0("of_",j,"_",formatC(e,width=2,flag="0"),".rdata"),of,ra,rb)

e<-"03"; load(paste0("tmp_e",e,"_j1.rdata")); rb1<-rb; load(paste0("tmp_e",e,"_j4.rdata"));image(rb1,breaks=c(0,0.1,0.25,0.5,1,2,4,8,16,32,64),col=c("gray",rev(rainbow(9))))
image(rbmod,breaks=c(0,0.1,0.25,0.5,1,2,4,8,16,32,64),col=c("gray",rev(rainbow(9))))
breaks=c(0,0.1,0.25,0.5,1,2,4,8,16,32,64); col=c("gray",rev(rainbow(9))); for (i in 1:length(col)) { ixx<-which(envtmp$obs_value>=breaks[i] & envtmp$obs_value<breaks[i+1]); points(envtmp$obs_x[ixx],envtmp$obs_y[ixx],pch=21,bg=col[i])}


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



#        print( paste( "j e of_par",j,e,of_par$par[1],of_par$par[2],of_par$par[3],
#                      round(range(getValues(of$u))[1],0),round(range(getValues(of$u))[2],0),
#                      round(range(getValues(of$v))[1],0),round(range(getValues(of$v))[2],0)))
      # Align smaller scales
      of_u <- of$u
      of_v <- of$v
      rm(of)
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

###
save(file=paste0("of_",j,"_",formatC(e,width=2,flag="0"),"_",jj,".rdata"),of_u,of_v)

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
      rm(rb_finer,rbmod_finer)
      r <- mrtree$raster[[j]]$r  
      r[] <- Ea[,e] * 2**jw
      dwt[[e]][[3*jw+1]][] <- as.matrix(r)
      rm(r)
###
save(file=paste0("dwt_",j,"_",formatC(e,width=2,flag="0"),".rdata"),dwt)
#      rm( rb_finer, of_u, of_v, rbmod_finer)
#e<-1;Eaj <- idwt.2d( dwt[[e]]);r <- mrtree$raster[[1]]$r;r[] <- array(data=as.matrix(Eaj),dim=c(sqrt(mrtree$m_dim[[1]]),sqrt(mrtree$m_dim[[1]])));r[r<0.1]<-0;image(r,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#s<-r;s[]<-mrbkg$data[[1]]$Eor[,e]; image(s,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#t<-r;t[]<-mrobs$idi[[1]]; image(t,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#png(file="test.png",width=1200,height=1200)
#image(r,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4)))) 

#save(file="tmp.rdata",r,env,dwt,mrbkg,envtmp,j,jw,mrtree,mrobs,argv,y_env)
#dev.off()
#q()
#    } # END - Loop over ensembles

q()


  # set dyadic domain
  xmn <- as.numeric(argv$grid_master.x1) - as.numeric(argv$grid_master.resx)/2
  xmx <- as.numeric(argv$grid_master.xn) + as.numeric(argv$grid_master.resx)/2
  ymn <- as.numeric(argv$grid_master.y1) - as.numeric(argv$grid_master.resy)/2
  ymx <- as.numeric(argv$grid_master.yn) + as.numeric(argv$grid_master.resy)/2
  nx <- ncol( env$rmaster)
  ny <- nrow( env$rmaster)

#  n <- ceiling( log( max(nx,ny), 2))
  n <- ceiling( log( min(nx,ny), 2))
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
###
argv$msa_ididense <- 0.1

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
  if ( !exists("jstop")) jstop <- jmax - 1 
  if ( exists( "robs" )) rm(robs)
  if ( exists( "ridi" )) rm(ridi)

###
save(file="multiresobs.rdata",argv,mrobs,mrtree)

  # safe-check, exit when no ok observations found 
  if (jstop == 0) return(NULL)

  dwt <- list() 
  envtmp$k_dim <- env$k_dim
  # Loop over spatial scales
  for (j in jstop:2) {
    cat( paste( "alligning spatial level", j))
    t0b <- Sys.time()
    jw <- j-1

    envtmp$Eb  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
    envtmp$D <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
    envtmp$Z  <- array( data=NA, dim=c( mrtree$m_dim[[j]], env$k_dim))
    envtmp$Y  <- array( data=NA, dim=c( mrobs$d_dim[[j]], env$k_dim))
    if (j == jstop) 
      Eor <- array( data=NA, dim=c( mrtree$m_dim[[1]], env$k_dim))

    # loop over ensemble members - prepare the background 
    for (e in 1:env$k_dim) {
      # first iteration, background is from the original ensemble
      if ( j == jstop ) {
        r <- resample( subset( fg_env$fg[[fg_env$ixf[fg_env$ixs[e]]]]$r_main, subset=fg_env$ixe[fg_env$ixs[e]]), rdyad, method="bilinear")
        if (!is.na(y_env$rain)) r[r<y_env$rain] <- 0
        Eor[,e] <- getValues(r)
      # 2nd,3rd,... iterations, background is from the previous iteration
      } else {
        aux <- idwt.2d( dwt[[e]])
        r <- s <- mrtree$raster[[1]]$r 
        r[] <- array(data=as.matrix(aux),dim=c(sqrt(mrtree$m_dim[[1]]),sqrt(mrtree$m_dim[[1]])))
        s[] <- Eor[,e]
        r <- clean_and_smooth_the_field( r, y_env$rain, s, 
                                         mrobs$val[[1]], mrobs$x[[1]], mrobs$y[[1]], 
                                         lambda=2, area_small_clumps=argv$area_small_clumps)
        rm(s)
      }
      #discrete wavelet transformation
      dwt[[e]] <- dwt.2d( as.matrix(r), wf=argv$msaensi_wf, J=jw, boundary=argv$msaensi_boundary)
      # father wavelet as the reconstructed field at the j-th spatial level
      r   <- mrtree$raster[[j]]$r
      r[] <- dwt[[e]][[3*jw+1]] / 2**jw
      # store the background in the multiresolution tree
      envtmp$Eb[,e]  <- getValues(r)
      envtmp$D[,e] <- mrobs$val[[j]] - extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]))
    } # end loop over ensemble members
    if ( j == jstop ) rm(rdyad) 
    rm(r)

    # background ensemble correlation matrices 
    Emean <- rowMeans(envtmp$Eb)
    Esd <- apply( envtmp$Eb, FUN=function(x){sd(x)}, MAR=1)
    for (e in 1:env$k_dim) { 
      # covariances
#      mrbkg$data[[j]]$X[,e] <- 1/sqrt(env$k_dim-1) * (envtmp$Eb[,e] - Emean)
      # correlations
      envtmp$Z[,e] <- 1/sqrt(env$k_dim-1) * (envtmp$Eb[,e] - Emean) / Esd
      envtmp$Z[,e][!is.finite(envtmp$Z[,e])] <- 1/sqrt(env$k_dim) 
      r   <- mrtree$raster[[j]]$r
      r[] <- envtmp$Z[,e]
      envtmp$Y[,e] <- extract( r, cbind( mrobs$x[[j]], mrobs$y[[j]]), method="simple")
    }
    rm( r, Emean, Esd)

###
save(file=paste0("envtmp_before_oi_",j,".rdata"),envtmp)

    # Spatial analysis with multiresolution background and observations
    envtmp$x <- mrtree$x[[j]]
    envtmp$y <- mrtree$y[[j]]
    envtmp$m_dim <- mrtree$m_dim[[j]]
    envtmp$obs_x <- mrobs$x[[j]]
    envtmp$obs_y <- mrobs$y[[j]]
    envtmp$eps2 <- rep( argv$msa_eps2, envtmp$m_dim) 
    envtmp$nn2 <- nn2( cbind(mrobs$x[[j]],mrobs$y[[j]]), 
                       query = cbind(mrtree$x[[j]],mrtree$y[[j]]), 
                       k = min( c(pmax,mrobs$d_dim[[j]])), 
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
                                           k_dim_corr=argv$k_dim_cor,
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
                                         k_dim_corr=argv$k_dim_cor,
                                         idi=F)))
    }
    Ea <- res[,1:env$k_dim]
    rm(res)
    if (any(is.na(Ea))) Ea[is.na(Ea)] <- 0
    envtmp$Y <- envtmp$Z <- envtmp$D <- envtmp$x <- envtmp$y <- NULL 
    envtmp$m_dim <- envtmp$obs_x <- envtmp$obs_y <- envtmp$eps2 <- NULL
    envtmp$nn2 <- NULL 

###
save(file=paste0("Ea_after_oi_",j,".rdata"),Ea)

    # Loop over ensembles
    ra <- mrtree$raster[[j]]$r
    rb <- mrtree$raster[[j]]$r
    for (e in 1:env$k_dim) {
      ra[] <- Ea[,e]
      rb[] <- envtmp$Eb[,e]
      of_nlevel <- min( c( floor(log2(nrow(rb))), floor(log2(ncol(rb)))))
      if ( floor(log2(nrow(rb))) == of_nlevel & floor(log2(ncol(rb))) == of_nlevel) 
        of_nlevel <- of_nlevel - 1
      of <- optical_flow_HS( rb, ra, nlevel=of_nlevel, niter=100, w1=100, w2=0, tol=0.00001)

###
save(file=paste0("of_",j,"_",formatC(e,width=2,flag="0"),".rdata"),of,ra,rb)

#        print( paste( "j e of_par",j,e,of_par$par[1],of_par$par[2],of_par$par[3],
#                      round(range(getValues(of$u))[1],0),round(range(getValues(of$u))[2],0),
#                      round(range(getValues(of$v))[1],0),round(range(getValues(of$v))[2],0)))
      # Align smaller scales
      of_u <- of$u
      of_v <- of$v
      rm(of)
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

###
save(file=paste0("of_",j,"_",formatC(e,width=2,flag="0"),"_",jj,".rdata"),of_u,of_v)

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
      rm(rb_finer,rbmod_finer)
      r <- mrtree$raster[[j]]$r  
      r[] <- Ea[,e] * 2**jw
      dwt[[e]][[3*jw+1]][] <- as.matrix(r)
      rm(r)
###
save(file=paste0("dwt_",j,"_",formatC(e,width=2,flag="0"),".rdata"),dwt)
#      rm( rb_finer, of_u, of_v, rbmod_finer)
#e<-1;Eaj <- idwt.2d( dwt[[e]]);r <- mrtree$raster[[1]]$r;r[] <- array(data=as.matrix(Eaj),dim=c(sqrt(mrtree$m_dim[[1]]),sqrt(mrtree$m_dim[[1]])));r[r<0.1]<-0;image(r,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#s<-r;s[]<-mrbkg$data[[1]]$Eor[,e]; image(s,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#t<-r;t[]<-mrobs$idi[[1]]; image(t,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4))))
#png(file="test.png",width=1200,height=1200)
#image(r,breaks=c(-100,0,0.1,1,2,4,8),col=c("beige","gray",rev(rainbow(4)))) 

#save(file="tmp.rdata",r,env,dwt,mrbkg,envtmp,j,jw,mrtree,mrobs,argv,y_env)
#dev.off()
#q()
    } # END - Loop over ensembles
    rm(Ea,ra,rb)
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
    r[] <- array(data=as.matrix(aux),dim=c(sqrt(mrtree$m_dim[[1]]),sqrt(mrtree$m_dim[[1]])))
    s[] <- Eor[,e]
    r <- clean_and_smooth_the_field( r, y_env$rain, s, 
                                     mrobs$val[[1]], mrobs$x[[1]], mrobs$y[[1]], 
                                     lambda=2, area_small_clumps=argv$area_small_clumps)
    rm(s)
    t <- resample( r, env$rmaster, method="ngb")
    rm(r)
    env$Xa[,e] <- getValues(t)
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a[,e] <- extract( t, cbind( y_env$yov$x, y_env$yov$y), method="simple")
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( t, cbind(  y_env$yo$x, y_env$yo$y), method="simple")
    rm(t)
  } # END loop over ensemble members

} # END FUNCTION
