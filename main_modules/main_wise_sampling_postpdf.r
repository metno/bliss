#+
main_wise_sampling_postpdf <- function( argv, y_env, fg_env, env, seed=NA,
                                        plot=F, dir_plot=NA) {
#
#------------------------------------------------------------------------------
  
  if (plot) {
    proj4.lcc<-"+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +no_defs +R=6.371e+06"
    b_utm33<-readOGR("/home/cristianl/data/geoinfo/TM_WORLD_BORDERS_UTM33/TM_WORLD_BORDERS_UTM33-0.2.shp","TM_WORLD_BORDERS_UTM33-0.2",verbose=F)
    b<-spTransform(b_utm33,CRS(proj4.lcc))
    rm(b_utm33,proj4.lcc)
    br<-c(0,1,2,4,8,16,32,64,128)
    col<-c("gray",rev(rainbow(7)))
  }

  options( warn = 2)

  library(waveslim)

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
  
  # Initialization of the wavelet structures
  #  dwt_out. dwt.2d class.
  #  dwt_ix. vector.
  #  dwt_res. vector.
  dwt_out <- dwt.2d( as.matrix(rdyad), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
  for (i in 1:env$n_levs_mx) {
    dwt_out[[3*(i-1)+1]][] <- i*10 + 1 # LHi - wavelet coefficients
    dwt_out[[3*(i-1)+2]][] <- i*10 + 2 # HLi - wavelet coefficients
    dwt_out[[3*(i-1)+3]][] <- i*10 + 3 # HHi - wavelet coefficients
  }
  # scaling function coefficients (father wavelet)
  dwt_out[[3*(env$n_levs_mx-1)+4]][] <- env$n_levs_mx*10 + 4 # LL env$n_levs_mx
  dwt_ix <- as.vector( unlist( dwt_out))
  for (i in 1:length(dwt_out)) dwt_out[[i]][] <- 0
  dwt_res <- 2**as.numeric(floor(dwt_ix/10))
  dwt_ix[dwt_ix<(env$n_levs_mn*10)] <- 0
  dwt_res[dwt_ix==(env$n_levs_mx+1)] <- 2**env$n_levs_mx
  # total number of coefficients used. Dimension of the state variable
  env$m_dim <- length(dwt_ix)
  cat( paste( "state variable, dim >", env$m_dim, "\n"))

  # Observations
  if (plot) {
    q   <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=rdyad, field=y_env$yo$value)
    q[is.na(q)] <- 0
    c_xy <- cellFromXY( rdyad, cbind(y_env$yo$x,y_env$yo$y))
    yo_check <- extract( q, c_xy)
    rm(q)
  }

  cat("Reconstruct analysis ensemble members on the state vector \n")
  Xa1     <- array( data=NA, dim=c(     env$m_dim, env$k_dim))
  for (e in 1:env$k_dim) {
    q<-rdyad
    q[]<-array(data=env$Xa[,e],dim=c(sqrt(length(env$Xa[,e])),sqrt(length(env$Xa[,e]))))
    Xa1[,e] <- as.vector( unlist( dwt.2d( as.matrix(q), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)))
    rm(q)
  }

  #
  cat("Sample realizations from the posterior distribution ")
  xa1        <- rowMeans(Xa1)
  Aa         <- Xa1 - xa1
  var_a1     <- 1/(env$k_dim-1) * rowSums( Aa*Aa) * 1/dwt_res**2
  env$Xa_large <- array( data=NA, dim=c( env$m_dim, env$n_dim))
  ya <- array( data=NA, dim=c( env$p_dim, env$n_dim))
  if (!is.na(seed)) set.seed(seed)
  for (n in 1:env$n_dim) {
    cat(".")
    aux     <- rnorm( env$m_dim, mean=xa1, sd=sqrt(var_a1) )
    dwt_aux <- dwt_out
    for (i in env$n_levs_mn:env$n_levs_mx) {
      ix <- which( dwt_ix == (i*10 + 1))
      dwt_aux[[3*(i-1)+1]][] <- aux[ix]
      ix <- which( dwt_ix == (i*10 + 2))
      dwt_aux[[3*(i-1)+2]][] <- aux[ix]
      ix <- which( dwt_ix == (i*10 + 3))
      dwt_aux[[3*(i-1)+3]][] <- aux[ix]
    }
    ix <- which( dwt_ix == (env$n_levs_mx*10 + 4))
    dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- aux[ix]
    env$Xa_large[,n] <- idwt.2d( dwt_aux)
    env$Xa_large[,n][env$Xa_large[,n]<y_env$rain]<-0

    if (plot) {
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
  cat("\n")

save(file=file.path(dir_plot,"tmp3.rdata"), argv, y_env, fg_env, u_env, env,rdyad,b)

  q()
}
