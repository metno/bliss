#+ Multiscale Alignment as a pre-proc for Ensemble Statistical Interpolation
msa <- function( argv, y_env, fg_env, env) {
#------------------------------------------------------------------------------
#  library(waveslim)
#  library(smoothie)

  t0a <- Sys.time()

  cat( "-- MSA --\n")

  # number of background ensemble members
  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

  # observations
  r <- env$rmaster
  
  envtmp$m_dim <- ncell(env$rmaster)
  xy <- xyFromCell( env$rmaster, 1:envtmp$m_dim)
  envtmp$x <- xy[,1]
  envtmp$y <- xy[,2]
  cat( paste( "         number of grid points, m dim >", envtmp$m_dim, "\n"))
  cat( paste( "number of in-situ observations, p dim >", y_env$yo$n, "\n"))

  # Loop over ensembles
  envtmp$Eb  <- array( data=NA, dim=c( ncell(env$rmaster), env$k_dim))
  ra <- env$rmaster; ra[] <- NA
  rb <- env$rmaster; rb[] <- NA

argv$dh_vec_min <- 2500
argv$dh_vec_max <- 200000
argv$dh_vec_by  <- 5000
argv$dh_r_max  <- 10000
argv$dh_obs_fact <- 3
argv$delta_sbe_mean_r <- 0.1
argv$sbe_filter <- 3
argv$pmax <- 20

  # Initialization
  env$Xa <- array( data=NA, dim=c( envtmp$m_dim, env$k_dim))

  dh_vec <- seq( argv$dh_vec_min, argv$dh_vec_max, by=argv$dh_vec_by)
  n_dh_vec <- length(dh_vec)
  safe_loop <- 10
  for (e in 1:env$k_dim) {
    cat( paste("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"))
    first <- T
    sbe_mean_prev <- NA
    eps2 <- 0.1
    for (ii in 1:safe_loop) {
      cat( paste("-----------------------------------------------------------\n"))
      cat( paste("ensemble",e,"/iteration",ii,"\n"))
      if (first) {
        rb <- subset( fg_env$fg[[fg_env$ixf[fg_env$ixs[e]]]]$r_main, 
                      subset=fg_env$ixe[fg_env$ixs[e]])
        if (!is.na(y_env$rain)) r[r<y_env$rain] <- 0
        first <- F
      } else {
        rb <- rbmod
        rb[is.na(rb)] <- 0
        if (!is.na(y_env$rain)) rb[rb<y_env$rain] <- 0
        env$Xa[,e] <- getValues(rb)
      }
      sbe  <- array( data=NA, dim=n_dh_vec)
      yb <- extract( rb, cbind( y_env$yo$x, y_env$yo$y), method="simple")
      if (!is.na(y_env$rain)) yb[yb<y_env$rain] <- 0
      for (j in 1:n_dh_vec) {
        fact <- ceiling( dh_vec[j] / mean(res(env$rmaster)))
        if (fact == 1) {
          raux_agg <- env$rmaster 
        } else {
          raux_agg <- aggregate(env$rmaster, fact=fact, expand=T)
        }
        rauxo <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                           field=y_env$yo$value, fun=mean, na.rm=T )
        rcounto <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                              field=y_env$yo$value, fun="count", na.rm=T )
        ix <- which( !is.na(getValues(rauxo)) & 
                     getValues(rcounto) > median(getValues(rcounto), na.rm=T))
        vo <- getValues(rauxo)[ix]
        rauxb <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                            field=yb, fun=mean, na.rm=T )
        vb <- getValues(rauxb)[ix]
        corr <- ( (mean(vo)-mean(vb))**2 + sd(vo)**2 + sd(vb)**2 - mean((vo-vb)**2)) /
                (2 * sd(vo) * sd(vb))
        sbe[j] <- 1 - sqrt( (corr-1)**2 + 
                  ((sd(vo)-sd(vb))/(sd(vo)+sd(vb)))**2 + 
                  ((mean(vo)-mean(vb))/(mean(vo)+mean(vb)))**2)
      } # end loop over spatial scale to determine sbe
      sbe_smooth <- as.vector( filter( sbe, filter=rep(1/argv$sbe_filter,argv$sbe_filter)))
      sbe_mean <- mean( range( sbe_smooth, na.rm=T))
      delta_sbe_mean <- abs(sbe_mean - sbe_mean_prev) / sbe_mean_prev
      if ( !is.na(delta_sbe_mean) & delta_sbe_mean < delta_sbe_mean_r ) break
      sbe_mean_prev <- sbe_mean
      dh <- min( c( dh_vec[min( which( sbe_smooth >= sbe_mean))], argv$dh_r_max))
      fact <- ceiling( (dh / argv$dh_obs_fact) / mean(res(env$rmaster)))
      cat( paste("dh pmax fact delta_sbe_mean",dh,argv$pmax,fact,round(100*delta_sbe_mean),"\n")) 
      if (fact == 1) {
        raux_agg <- env$rmaster 
      } else {
        raux_agg <- aggregate(env$rmaster, fact=fact, expand=T)
      }
      raux <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                         field=y_env$yo$value, fun=max, na.rm=T )
      ix <- which( !is.na(getValues(raux)))
      xy <- xyFromCell( raux, 1:ncell(raux))
      envtmp$obs_value <- getValues(raux)[ix]
#    envtmp$obs_x <- xy[ix,1]
#    envtmp$obs_y <- xy[ix,2]
      raux_x <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                           field=y_env$yo$x, fun=mean, na.rm=T )
      envtmp$obs_x <- getValues(raux_x)[ix]
      raux_y <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
                           field=y_env$yo$y, fun=mean, na.rm=T )
      envtmp$obs_y <- getValues(raux_y)[ix]
      envtmp$obs_n <- length(ix)
      cat( paste( "number of rasterized observations, p dim >", envtmp$obs_n, "\n"))
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

#        cat( paste("-----------------------------------------------------------\n"))
#        cat( paste("m1 m2 range / m_dim",m1,m2,envtmp$m_dim,"/",m_dim,"\n"))
        t00 <- Sys.time()
        envtmp$nn2 <- nn2( cbind( envtmp$obs_x, envtmp$obs_y), 
                           query = cbind( envtmp$x, envtmp$y), 
                           k = argv$pmax, searchtype = "radius", 
                           radius = (7*dh))
        t11 <- Sys.time()
#        cat( paste( "nn2 total time", round(t11-t00,1), attr(t11-t00,"unit"), "\n"))
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
#        cat( paste( "oi total time", round(t11-t00,1), attr(t11-t00,"unit"), "\n"))
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
#        of <- optical_flow_HS( rb_norm, ra_norm, nlevel=of_nlevel, 
      of <- optical_flow_HS( rb, ra, nlevel=of_nlevel, 
                             niter=100, w1=100, w2=0, tol=0.00001)
      t11 <- Sys.time()
      cat( paste( "optflow of_nlevel", of_nlevel, "\n"))
      cat( paste( "optflow total time", round(t11-t00,1), attr(t11-t00,"unit"), "\n"))
      u <- env$rmaster
      v <- env$rmaster
#      u[] <- t( matrix( data=as.vector(of$u[,,1]), 
#                ncol=nrow(env$rmaster), nrow=ncol(env$rmaster)))
#      v[] <- t( matrix( data=as.vector(of$v[,,1]), 
#                ncol=nrow(env$rmaster), nrow=ncol(env$rmaster)))
#      rbmod <- warp( rb, -u, -v, method="simple")
      rbmod <- warp( rb, -of$u, -of$v, method="simple")
save(file=paste0("tmp_e",formatC(e,width=2,flag="0"),"_ii",ii,".rdata"),argv,raux,ra,rb,ra_norm,rb_norm,res,y_env,xa,envtmp,env,of,u,v,rbmod,robs_on_master,sbe,sbe_smooth,sbe_mean,dh_r)
    } # loop over ii
  }
q()
source("~/projects/bliss/functions/optflow/optflow_util.r")
###
save(file=paste0("of_",j,"_",formatC(e,width=2,flag="0"),".rdata"),of,ra,rb)

e<-"03"; load(paste0("tmp_e",e,"_j40.rdata")); rb1<-rb; load(paste0("tmp_e",e,"_j4.rdata"));image(rb1,breaks=c(0,0.1,0.25,0.5,1,2,4,8,16,32,64),col=c("gray",rev(rainbow(9))))
image(rbmod,breaks=c(0,0.1,0.25,0.5,1,2,4,8,16,32,64),col=c("gray",rev(rainbow(9))))
breaks=c(0,0.1,0.25,0.5,1,2,4,8,16,32,64); col=c("gray",rev(rainbow(9))); for (i in 1:length(col)) { ixx<-which(envtmp$obs_value>=breaks[i] & envtmp$obs_value<breaks[i+1]); points(envtmp$obs_x[ixx],envtmp$obs_y[ixx],pch=21,bg=col[i])}

  # Initialization
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$k_dim))
  if (env$cv_mode | env$cv_mode_random) {
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$k_dim))
  }

  # loop over ensemble members
  r <- env$rmaster
  for (e in 1:env$k_dim) {
    r[] <- env$Xa[,e]
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a[,e] <- extract( r, cbind( y_env$yov$x, y_env$yov$y), 
                                        method="simple")
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y), 
                                      method="simple")
  } # END loop over ensemble members

} # END FUNCTION


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


  # structure functions
#  raino <- y_env$yo$value; raino[] <- NA
#  rainb <- y_env$yo$value; rainb[] <- NA
#  raino[y_env$yo$value>=0.1] <- 1
#  raino[y_env$yo$value< 0.1] <- 0
#  dh_vec <- seq(2500,200000,by=5000)
#  rmse <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
#  corr <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
#  sbe  <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
#  nse  <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
#  kge  <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
#  ss   <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
#  fss   <- array(data=NA,dim=c(length(dh_vec),env$k_dim))
#  for (e in 1:env$k_dim) {
#    rb <- subset( fg_env$fg[[1]]$r_main, subset=e)
#    yb <- extract( rb, cbind(y_env$yo$x,y_env$yo$y))
#    rainb[yb>=0.1] <- 1
#    rainb[yb< 0.1] <- 0
##    dist <- sqrt( outer( y_env$yo$x,y_env$yo$x,FUN="-")**2 + outer( y_env$yo$y,y_env$yo$y,FUN="-")**2)
##save(file="tmp.rdata",yb,y_env,dist)
#    cat( paste("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"))
#    cat( paste("ensemble",e,"\n"))
#    for (j in 1:length(dh_vec)) {
##      cat( paste("===========================================================\n"))
#      fact <- ceiling( dh_vec[j] / mean(res(env$rmaster)))
##      cat( paste("dh pmax fact",dh,fact,"\n"))
#      if (fact == 1) {
#        raux_agg <- env$rmaster 
#      } else {
#        raux_agg <- aggregate(env$rmaster, fact=fact, expand=T)
#      }
#      rauxo <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
#                         field=y_env$yo$value, fun=mean, na.rm=T )
#      rcounto <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
#                            field=y_env$yo$value, fun="count", na.rm=T )
#      ix <- which( !is.na(getValues(rauxo)) & getValues(rcounto)>median(getValues(rcounto),na.rm=T))
#      vo <- getValues(rauxo)[ix]
#      rauxb <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=raux_agg, 
#                          field=yb, fun=mean, na.rm=T )
#      vb <- getValues(rauxb)[ix]
#      rmse[j,e] <- sqrt(mean((vo-vb)**2))
#      corr[j,e] <- ( (mean(vo)-mean(vb))**2 + sd(vo)**2 + sd(vb)**2 - mean((vo-vb)**2)) / (2 * sd(vo) * sd(vb))
#      sbe[j,e] <- 1 - sqrt( (corr[j,e]-1)**2 + ((sd(vo)-sd(vb))/(sd(vo)+sd(vb)))**2 + ((mean(vo)-mean(vb))/(mean(vo)+mean(vb)))**2)
#      ss[j,e] <-  2 * sd(vo) * sd(vb) * corr[j,e] / ((mean(vo)-mean(vb))**2 + sd(vo)**2 + sd(vb)**2)
#      nse[j,e] <- 1 - mean((vo-vb)**2)/sd(vo)**2
#      kge[j,e] <- 1 - sqrt( (corr[j,e]-1)**2 + (sd(vb)/sd(vo)-1)**2 + (mean(vb)/mean(vo)-1)**2)
#      fss[j,e] <- 1 - mean( (vo - vb)**2) / (mean(vo**2) + mean(vb**2))
#    }
#  }
#  sbe_smooth <- as.vector( filter( rowMeans(sbe), filter=rep(1/3,3)))
#  sbe_r <- mean( range( sbe_smooth, na.rm=T))
#  dh_r <- dh_vec[min( which( sbe_smooth >= sbe_r))]
##save(file="tmp.rdata",raino,rainb,rmse,dh_vec,corr,sbe,ss,nse,kge,fss)
## plot(dh_vec,rmse[,1],ylim=c(0,3),col="white"); for (i in 1:30) lines(dh_vec,rmse[,i])
## plot(dh_vec,rmse[,1],ylim=c(0.1,0.5),col="white"); for (i in 1:30) lines(dh_vec,rmse[,i])
##points(dh_vec,rowMeans(rmse),pch=21,bg="red")
## try with spatial correlaiton instead of rmse (see Casati's euqation)
##q()
