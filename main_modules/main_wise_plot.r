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

  # raster cell from observation XY values
  c_xy <- cellFromXY( rdyad, cbind(y_env$yo$x,y_env$yo$y))
  q   <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=rdyad, field=y_env$yo$value)
  q[is.na(q)] <- 0
  yo_check <- extract( q, c_xy)


  proj4.lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
  b_utm33<-readOGR("/home/cristianl/data/geoinfo/TM_WORLD_BORDERS_UTM33/TM_WORLD_BORDERS_UTM33-0.2.shp","TM_WORLD_BORDERS_UTM33-0.2",verbose=F)
  b<-spTransform(b_utm33,CRS(proj4.lcc))
  rm(b_utm33,proj4.lcc)

  br<-c(0,1,2,4,8,16,32,64,128)
  col<-c("gray",rev(rainbow(7)))
jump<-F
if (!jump) {
  cat("plotting background ")
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
  cat("\n")

  cat("plotting combo ")
  for (e in 1:env$k_dim) {
    cat(".")
    q<-rdyad
    q[]<-array(data=env$Xa[,e],dim=c(sqrt(length(env$Xa[,e])),sqrt(length(env$Xa[,e]))))
    fout<-file.path(dir_plot,paste0("wise_main_xcombo_e",formatC(e,width=2,flag="0"),".png"))
    fb<-file.path(dir_plot,paste0("wise_main_xb_e",formatC(e,width=2,flag="0"),".png"))
    fa<-file.path(dir_plot,paste0("wise_main_xa_e",formatC(e,width=2,flag="0"),".png"))
    png(file=fa,width=1200,height=1200)
    image(q,breaks=br,col=col,main="analysis")
    plot(b,add=T)
#    for (i in 1:length(col)) {
#      ix <- which(yo_check>=br[i] & yo_check<br[i+1])
#      cex<-1; col_i <- col[i]
#      if (i==1) { cex<-0.5; col_i<-"beige" }
#      points(y_env$yo$x[ix],y_env$yo$y[ix],pch=21,bg=col_i,col=col_i,cex=cex)
#    }
    dev.off()
    system( paste0("convert +append ",fb," ",fa," ",fout))
    system( paste0("rm ",fa," ",fb))

    ya <- extract( q, c_xy)
    i<-fg_env$ixs[e]
    s <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
    yb <- extract( s, c_xy) 
    yo <- yo_check  
    fout<-file.path(dir_plot,paste0("wise_main_ycombo_e",formatC(e,width=2,flag="0"),".png"))
    fb<-file.path(dir_plot,paste0("wise_yb_e",formatC(e,width=2,flag="0"),".png"))
    fa<-file.path(dir_plot,paste0("wise_ya_e",formatC(e,width=2,flag="0"),".png"))
    png(file=fb,width=1200,height=1200)
    plot(yo,yb,ylim=range(c(yo,ya,yb)),xlim=range(c(yo,ya,yb)),main="background vs observations",pch=21,bg="cornflowerblue",col="blue")
    lines(-10000:10000,-10000:10000)
    seq<-c(0,0.1,seq(0.5,100,by=0.5))
    for (ii in 1:(length(seq)-1)) {
      ix <- which( yo>=seq[ii] & yo<seq[ii+1])
      if (length(ix)<10) next
      points( (seq[ii]+seq[ii+1])/2, median(yb[ix]), cex=2, pch=21,bg="black")
      points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(yb[ix],probs=0.25)), cex=1, pch=21,bg="black")
      points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(yb[ix],probs=0.75)), cex=1, pch=21,bg="black")
    }
    dev.off()
    png(file=fa,width=1200,height=1200)
    plot(yo,ya,ylim=range(c(yo,ya,yb)),xlim=range(c(yo,ya,yb)),main="analysis vs observations",pch=21,bg="pink",col="red")
    lines(-10000:10000,-10000:10000)
    for (ii in 1:(length(seq)-1)) {
      ix <- which( yo>=seq[ii] & yo<seq[ii+1])
      if (length(ix)<10) next
      points( (seq[ii]+seq[ii+1])/2, median(ya[ix]), cex=2, pch=21,bg="black")
      points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(ya[ix],probs=0.25)), cex=1, pch=21,bg="black")
      points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(ya[ix],probs=0.75)), cex=1, pch=21,bg="black")
    }
    dev.off()
    system( paste0("convert +append ",fb," ",fa," ",fout))
    system( paste0("rm ",fa," ",fb))
print(fout)
    rm(q)
  }
  cat("\n")
}
  ya <- array( data=NA, dim=c( env$p_dim, env$n_dim))

  for (n in 1:env$n_dim) {
    q<-rdyad
    q[]<-array(data=env$Xa_large[,n],dim=c(sqrt(length(env$Xa_large[,n])),sqrt(length(env$Xa_large[,n]))))
    fout<-file.path(dir_plot,paste0("wise_xa_e",formatC(n,width=2,flag="0"),".png"))
    fm<-file.path(dir_plot,paste0("wise_xa_m_e",formatC(n,width=2,flag="0"),".png"))
    png(file=fm,width=1200,height=1200)
    image(q,breaks=br,col=col)
    plot(b,add=T)
    for (i in 1:length(col)) {
      ix <- which(yo_check>=br[i] & yo_check<br[i+1])
      cex<-1; col_i <- col[i]
      if (i==1) { cex<-0.5; col_i<-"beige" }
      points(y_env$yo$x[ix],y_env$yo$y[ix],pch=21,bg=col_i,col=col_i,cex=cex)
    }
    dev.off()

    ya[,n] <- extract( q, c_xy)
    fx<-file.path(dir_plot,paste0("wise_xa_x_e",formatC(n,width=2,flag="0"),".png"))
    png(file=fx,width=1200,height=1200)
    plot(y_env$yo$x,ya[,n],ylim=range(c(yo_check,ya[,n])),pch=21,bg="pink",col="red")
    points(y_env$yo$x,yo_check)
    abline(h=0.1)
    abline(h=1:100,lty=2,col="gray")
    abline(h=y_env$rain,col="green")
    dev.off()
    fy<-file.path(dir_plot,paste0("wise_xa_y_e",formatC(n,width=2,flag="0"),".png"))
    png(file=fy,width=1200,height=1200)
    plot(y_env$yo$y,ya[,n],ylim=range(c(yo_check,ya[,n])),pch=21,bg="cornflowerblue",col="blue")
    points(y_env$yo$y,yo_check)
    dev.off()
    system( paste0("convert +append ",fm," ",fx," ",fy," ",fout))
    system( paste0("rm ",fm," ",fx," ",fy))
    rm(q)
  }


}
