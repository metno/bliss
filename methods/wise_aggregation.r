#+
wise_aggregation <- function( argv, y_env, fg_env, u_env, env, 
                              plot=F,dir_plot=NA) {

    # set dyadic domain
    xmn <- as.numeric(argv$grid_master.x1) - as.numeric(argv$grid_master.resx)/2
    xmx <- as.numeric(argv$grid_master.xn) + as.numeric(argv$grid_master.resx)/2
    ymn <- as.numeric(argv$grid_master.y1) - as.numeric(argv$grid_master.resy)/2
    ymx <- as.numeric(argv$grid_master.yn) + as.numeric(argv$grid_master.resy)/2
    nx <- ncol( env$rmaster)
    ny <- nrow( env$rmaster)
    n <- ceiling( log( max(nx,ny), 2))
    rdyad <- raster( extent( xmn, xmx, ymn, ymx), ncol=2**n, nrow=2**n, crs=argv$grid_master.proj4)
    rdyad[] <- 0
    sqrt_m_dim <- sqrt( env$m_dim)

    # 
    env$Xa <- array( data=NA, dim=c( env$ngrid, env$a_dim))
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$a_dim)) 
    y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$a_dim)) 

    for (a in 1:env$a_dim) {
      t0 <- Sys.time()
      rdyad[] <- array(data=env$Xa_dyad[,a],dim=c(sqrt_m_dim,sqrt_m_dim))
      r <- resample( rdyad, env$rmaster, method="ngb")
      env$Xa[,a] <- getValues(r)
      if (env$cv_mode | env$cv_mode_random) 
        y_env$yov$value_a[,a] <- extract( rdyad, cbind( y_env$yov$x, y_env$yov$y))
      y_env$yo$value_a[,a]  <- extract( rdyad, cbind(  y_env$yo$x, y_env$yo$y))
#      val<-getValues(rdyad)
#      xy<-xyFromCell(rdyad,env$mask)
#      ix_wet <- which( val >= y_env$rain)
##      r<-rasterize(x=xy,y=env$rmaster,field=val[env$mask],fun=function(x,...){quantile(x,probs=0.9,na.rm=T)})
##      r<-rasterize(x=xy,y=env$rmaster,field=val[env$mask],fun=mean)
#      if (a==1) {
#        rfobs<-rdyad
#        nn2 <- nn2( xyFromCell( rdyad, 1:ncell(rdyad)), 
#                    query = xyFromCell( env$rmaster, 1:ncell(env$rmaster)), 
#                    k = 1000, 
#                    searchtype = "radius", radius = (res(rdyad)[1]*sqrt(2)))
#        mat <- nn2[[1]]
#        print(dim(mat))
#        c_xy <- which( ( aux <- rowSums( mat)) > 0 )
#        rm(nn2)
#        r<-env$rmaster; r[]<-NA
#      }  
#      mapply_quantile  <- function(i) { if ((i%%10000)==0) print(i);quantile( val[mat[c_xy[i],1:length(which(mat[c_xy[i],]!=0))]], probs=0.9) }
#      r[] <- 0
#      r[c_xy] <- t( mapply( mapply_quantile, 1:length(c_xy), SIMPLIFY = T))
      print(paste(a,Sys.time()-t0))
    }
}

