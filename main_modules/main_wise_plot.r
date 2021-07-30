main_wise_plot <- function( argv, y_env, fg_env, env, dir_plot) {

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

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))

  proj4.lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
  b_utm33<-readOGR("/home/cristianl/data/geoinfo/TM_WORLD_BORDERS_UTM33/TM_WORLD_BORDERS_UTM33-0.2.shp","TM_WORLD_BORDERS_UTM33-0.2",verbose=F)
  b<-spTransform(b_utm33,CRS(proj4.lcc))
  rm(b_utm33,proj4.lcc)

  br<-c(0,1,2,4,8,16,32,64,128)
  col<-c("gray",rev(rainbow(7)))

  j<-0
  for (i in fg_env$ixs) {
    cat(".")
    j <- j + 1
    # interpolate onto the dyadic grid
    s <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
    s[s<0] <- 0
    fout<-file.path(dir_plot,paste0("wise_main_xb_e",formatC(j,width=2,flag="0"),".png"))
    png(file=fout,width=1200,height=1200)
    image(s,breaks=br,col=col,main="background")
    plot(b,add=T)
    for (i in 1:length(col)) {
      ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
      cex<-1; col_i <- col[i]
      if (i==1) { cex<-0.5; col_i<-"beige" }
      points(y_env$yo$x[ix],y_env$yo$y[ix],pch=21,bg=col_i,col=col_i,cex=cex)
    }
    dev.off()
  }

}
