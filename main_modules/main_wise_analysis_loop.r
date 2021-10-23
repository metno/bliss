#+
main_wise_analysis_loop <- function( argv, y_env, fg_env, env, seed=NA, obs_k_dim=1000, 
                                plot=F, dir_plot=NA) {
#
#------------------------------------------------------------------------------
dp <- "/home/cristianl/data/wise/pngs"
plot<-F

#+ position of coefficients of a given spatial resolution level into the vector of all wavelet coefficients
ijFromLev <- function( n, lev, fw=F) { 
# n = max number of levels (i.e. 2**n is the length of one dimension of the dyadic domain)
# lev = spatial resolution level (1=smallest; ... to largest)
# fw = are these the father wavelet coefficients?
#-------------------------------------------------------------------
  i <- c( NA, NA, NA)
  j <- c( NA, NA, NA)
  if ( lev == 1) { 
    i[1] <- 1
  } else {
    i[1] <- sum( 3 * ( 2**n / 2**(1:(lev-1)))**2) + 1
  }
  i[2] <- i[1] + ( 2**n / 2**lev)**2
  i[3] <- i[2] + ( 2**n / 2**lev)**2
  j[3] <- i[1] + 3 * ( 2**n / 2**lev)**2 - 1
  j[2] <- j[3] - ( 2**n / 2**lev)**2
  j[1] <- j[2] - ( 2**n / 2**lev)**2 
  if (fw) {
    i <- j[3] + 1
    j <- i + ( 2**n / 2**lev)**2 - 1
    return( c(i,j))
  }
  return( cbind(i,j))
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
  
  # raster cell from observation XY values

#   uncomment the following lines to use radar
#  radarval <- getValues( u_env$uo[[1]]$r_main)
#  ixrad <- which(!is.na(radarval))
#   
#  y_env$yo$x <- c( y_env$yo$x, xyFromCell(u_env$uo[[1]]$r_main,ixrad)[,1])
#  y_env$yo$y <- c( y_env$yo$y, xyFromCell(u_env$uo[[1]]$r_main,ixrad)[,2])
#  y_env$yo$value <- c( y_env$yo$value, radarval[ixrad])

#  c_xy <- cellFromXY( rdyad, cbind(y_env$yo$x,y_env$yo$y))
  c_xy <- fourCellsFromXY( rdyad, cbind(y_env$yo$x,y_env$yo$y))
  
  # Initialization of the wavelet structures
  #  dwt_out. dwt.2d class.
  #  dwt_ix. vector.
  #  dwt_res. vector.
env$n_levs_mx<-10
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

  # define constants
  # total number of coefficients used. Dimension of the state variable
  env$n_dim <- length(dwt_ix)
  # total number of coefficients used. Dimension of the state variable
  env$m_dim <- length( getValues(rdyad))
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( "state variable, wavelet coefficients, n dim >", env$n_dim, "\n"))
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))

  # total number of coefficients used. Dimension of the state variable
  env$m_dim <- length(dwt_ix)
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( "state variable, dim >", env$m_dim, "\n"))

  # ---~--------------
  # -- Observations --

  cat("Transform observations and compute error covariance matrix ")
 
  rasterize_with_buffer <- function( xy, r, l) {
    
  }

  # rasterize
  rfobs <- rdyad
  for (c in 1:4) rfobs[c_xy[,c]] <- y_env$yo$value
  rfidi <- rdyad
  for (c in 1:4) rfidi[c_xy[,c]] <- 1

  # transformation operator
  dwt     <- dwt.2d( as.matrix(rfobs), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
  dwt_idi <- dwt.2d( as.matrix(rfidi), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

  # unlist the result and compute squared energies
  yo1   <- vector( mode="numeric", length=env$n_dim); yo1[]<-NA
  yidi1 <- vector( mode="numeric", length=env$n_dim); yidi1[]<-NA
  env$En2idi <- vector( mode="numeric", length=env$n_dim); env$En2idi[]<-NA
  env$En2obs <- vector( mode="numeric", length=(env$n_levs_mx+1)); env$En2obs[]<-NA
  jj <- 0; ii <- 0
  for (l in 1:env$n_levs_mx) {
    lh <- dwt[[1+3*(l-1)]]; hl <- dwt[[2+3*(l-1)]]; hh <- dwt[[3+3*(l-1)]]
    lh_idi <- dwt_idi[[1+3*(l-1)]]; hl_idi <- dwt_idi[[2+3*(l-1)]]; hh_idi <- dwt_idi[[3+3*(l-1)]]
    ii <- jj + 1
    jj <- ii + length(lh) + length(hl) + length(hh) - 1
    yo1[ii:jj] <- c( lh, hl, hh)
    yidi1[ii:jj] <- c( lh_idi, hl_idi, hh_idi)
    env$En2idi[l] <- mean( (lh_idi / 2**l)**2) + mean( (hl_idi / 2**l)**2) + mean( (hh_idi / 2**l)**2)
    env$En2obs[l] <- mean( ( lh / 2**l)**2)     + mean( ( hl / 2**l)**2)     + mean( ( hh / 2**l)**2)
  }
  ll <- dwt[[4+3*(env$n_levs_mx-1)]] 
  yo1[(jj+1):env$n_dim] <- ll
  env$En2obs[env$n_levs_mx+1] <- mean( ( ll / 2**l)**2)
  ll_idi <- dwt_idi[[4+3*(env$n_levs_mx-1)]] 
  yidi1[(jj+1):env$n_dim] <- ll_idi
  env$En2idi[env$n_levs_mx+1] <- mean( (ll_idi/2**env$n_levs_mx)**2)
  rm( lh, hl, hh, ll)
  rm( lh_idi, hl_idi, hh_idi, ll_idi)
  cat("\n")

  # ---~------------
  # -- Background --

  t0 <- Sys.time()
  cat("Transform background")
  Xb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  Yb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  Yd1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  j <- 0
  env$En2xb <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2yb <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2in <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  fxb <- vector()
  for (i in fg_env$ixs) {
    cat(".")
    j <- j + 1

    # interpolate onto the dyadic grid
    # background at grid  points
    rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
    rfxb[rfxb<0] <- 0
    # background at observation points
    rfyb  <- rfxb; rfyb[] <- 0; for (c in 1:4) rfyb[c_xy[,c]] <- extract( rfxb, c_xy[,c])
    # innovation 
    rfin  <- rfxb; rfin[] <- 0; for (c in 1:4) rfin[c_xy[,c]] <- extract( rfobs, c_xy[,c]) - extract( rfxb, c_xy[,c])

if (plot) {
    fxb[j] <- file.path(dp,paste0("rfxb_",formatC(j,flag="0",width=2),".png"))
    br<-c(0,1,2,4,8,16,32,64,128)
    col<-c("lightgray",rev(rainbow(7)))
    png(file=fxb[j],width=1200,height=1200)
    image(rfxb,breaks=br,col=col,main="background")
    dev.off()
}

    # transformation operator
    dwtxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
    dwtyb <- dwt.2d( as.matrix(rfyb), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
    dwtin <- dwt.2d( as.matrix(rfin), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
    # unlist the result and compute squared energies
    jj <- 0; ii <- 0
    for (l in 1:env$n_levs_mx) {
      lh <- dwtxb[[1+3*(l-1)]]; hl <- dwtxb[[2+3*(l-1)]]; hh <- dwtxb[[3+3*(l-1)]]
      env$En2xb[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
      ii <- jj + 1
      jj <- ii + length(lh) + length(hl) + length(hh) - 1
      Xb1[ii:jj,j] <- c( lh, hl, hh)
      lh  <-  dwtyb[[1+3*(l-1)]]; hl  <-  dwtyb[[2+3*(l-1)]];  hh <-  dwtyb[[3+3*(l-1)]]
      env$En2yb[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
      Yb1[ii:jj,j] <- c( lh, hl, hh)
      lh  <-  dwtin[[1+3*(l-1)]]; hl  <-  dwtin[[2+3*(l-1)]];  hh <-  dwtin[[3+3*(l-1)]]
      env$En2in[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
      Yd1[ii:jj,j] <- c( lh, hl, hh)
    }
    ll <- dwtxb[[4+3*(env$n_levs_mx-1)]] 
    env$En2xb[env$n_levs_mx+1,j] <- mean( (ll/2**env$n_levs_mx)**2)
    Xb1[(jj+1):env$n_dim,j] <- ll
    ll <- dwtyb[[4+3*(env$n_levs_mx-1)]] 
    env$En2yb[env$n_levs_mx+1,j] <- mean( (ll/2**env$n_levs_mx)**2)
    Yb1[(jj+1):env$n_dim,j] <- ll
    ll <- dwtin[[4+3*(env$n_levs_mx-1)]] 
    env$En2in[env$n_levs_mx+1,j] <- mean( (ll/2**env$n_levs_mx)**2)
    Yd1[(jj+1):env$n_dim,j] <- ll
    rm( dwtxb, dwtyb, dwtin, lh, hl, hh, ll, ii, jj)

  } # end loop over background fields
  t1 <- Sys.time()
  cat( paste0("time=",round(t1-t0,1), attr(t1-t0,"unit"),"\n"))
if (plot) {
  for (j in 1:env$k_dim) {
    ylim<-range(c(env$En2xb[,j],env$En2obs,env$En2yb[,j]))
    yylim<-range(c(env$En2obs,env$En2yb[,j]))
    fout<-file.path(dp,paste0("en2_",formatC(0,width=2,flag="0"),"_",formatC(j,width=2,flag="0"),".png"))
    png(file=fout,width=800,height=800)
    plot(1:env$n_levs_mx,env$En2xb[1:env$n_levs_mx,j],type="l",col="red",ylim=ylim)
    points(env$n_levs_mx,env$En2xb[env$n_levs_mx+1,j],pch=21,bg="red",cex=2)
    lines(1:env$n_levs_mx,env$En2obs[1:env$n_levs_mx],col="blue")
    lines(1:env$n_levs_mx,env$En2obs[1:env$n_levs_mx],col="blue")
    lines(1:env$n_levs_mx,env$En2yb[1:env$n_levs_mx,j],col="pink4")
    points(env$n_levs_mx,env$En2obs[env$n_levs_mx+1],pch=21,bg="blue",cex=2)
    points(env$n_levs_mx,env$En2yb[env$n_levs_mx+1,j],pch=21,bg="pink4",cex=2)
    par(new=T)
    plot(1:env$n_levs_mx,env$En2in[1:env$n_levs_mx],col="gold",axes=F,type="l",lty=1,lwd=2)
    points(env$n_levs_mx,env$En2in[env$n_levs_mx+1],pch=21,bg="gold",cex=2)
    axis(4)
    par(new=T)
    plot(1:env$n_levs_mx,env$En2idi[1:env$n_levs_mx],col="black",axes=F,type="l",lty=2)
    dev.off()
    print(paste("written file",fout))
  }
}

  # Analysis
  cat("Analysis on the transformed space \n")

  costf <- vector()
  ya <- extract( rfxb, cbind( y_env$yo$x, y_env$yo$y))

  for (loop in 1:20) {
    t0 <- Sys.time()
    cat( paste(" loop", loop, " "))
    Xb1[abs(Xb1)<(1e-02)]<-0
    Yb1[abs(Yb1)<(1e-02)]<-0
#    yo1[abs(yo1)<(1e-02)]<-0
    # Background error covariance matrices
    Ab <- Xb1 - rowMeans(Xb1)
    Db <- Yb1 - rowMeans(Yb1)
#    D  <- Yd1 - rowMeans(Yd1)
    D  <- Yd1
    var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) #* 1/dwt_res**2 
    var_b1_yy <- 1/(env$k_dim-1) * rowSums( Db*Db) #* 1/dwt_res**2
    var_bd1_xy <- 1/(env$k_dim-1) * rowSums( Ab*D) #* 1/dwt_res**2 
    var_d1_yy  <- 1/(env$k_dim-1) * rowSums( D *D) #* 1/dwt_res**2
    rm(Ab,Db,D)
  
    costf[loop] <- mean( sqrt( var_d1_yy))
    print( paste("rmse", costf[loop]))
    
    coeff <- var_b1_xy / var_d1_yy
    coeff[!is.finite(coeff)] <- 0
    Xa1     <- array( data=NA, dim=c(     env$n_dim, env$k_dim))
    env$Xa  <- array( data=NA, dim=c( length(rdyad), env$k_dim))
    for (e in 1:env$k_dim) {
      cat(".")
      Xa1[,e] <- Xb1[,e] + coeff * ( yo1 - Yb1[,e])
      dwt_aux <- dwt_out
      for (i in env$n_levs_mn:env$n_levs_mx) {
        ij <- ijFromLev( env$n_levs_mx, i, F)
#        ix <- which( dwt_ix == (i*10 + 1))
        dwt_aux[[3*(i-1)+1]][] <- Xa1[ij[1,1]:ij[1,2],e]
#        ix <- which( dwt_ix == (i*10 + 2))
        dwt_aux[[3*(i-1)+2]][] <- Xa1[ij[2,1]:ij[2,2],e]
#        ix <- which( dwt_ix == (i*10 + 3))
        dwt_aux[[3*(i-1)+3]][] <- Xa1[ij[3,1]:ij[3,2],e]
      }
      ij <- ijFromLev( env$n_levs_mx, env$n_levs_mx, T)
#      ix <- which( dwt_ix == (env$n_levs_mx*10 + 4))
      dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- Xa1[ij[1]:ij[2],e]
      env$Xa[,e] <- idwt.2d( dwt_aux)
      rm(dwt_aux)
      env$Xa[,e][env$Xa[,e]<y_env$rain]<-0
    }
    t1 <- Sys.time()
    cat( paste0("time=",round(t1-t0,1), attr(t1-t0,"unit"),"\\"))

    Xb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    Yb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    Yd1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    for (e in 1:env$k_dim) {
      cat(".")

      # -- background at grid  points --
#      rfxb_to_plot[] <- array(data=env$Xa[,e],dim=c(sqrt(length(env$Xa[,e])),sqrt(length(env$Xa[,e]))))
      rfxb[] <- array(data=env$Xa[,e],dim=c(sqrt(length(env$Xa[,e])),sqrt(length(env$Xa[,e]))))
      rfxb[rfxb<0] <- 0
      ya <- extract( rfxb, cbind( y_env$yo$x, y_env$yo$y))
      # background at observation points
      rfyb  <- rfxb; rfyb[] <- 0; for (c in 1:4) rfyb[c_xy[,c]] <- extract( rfxb, c_xy[,c])
      # innovation 
      rfin  <- rfxb; rfin[] <- 0; for (c in 1:4) rfin[c_xy[,c]] <- extract( rfobs, c_xy[,c]) - extract( rfxb, c_xy[,c])

      # transformation operator
      dwtxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
      dwtyb <- dwt.2d( as.matrix(rfyb), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
      dwtin <- dwt.2d( as.matrix(rfin), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in 1:env$n_levs_mx) {
        lh <- dwtxb[[1+3*(l-1)]]; hl <- dwtxb[[2+3*(l-1)]]; hh <- dwtxb[[3+3*(l-1)]]
        env$En2xb[l,e] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        ii <- jj + 1
        jj <- ii + length(lh) + length(hl) + length(hh) - 1
        Xb1[ii:jj,e] <- c( lh, hl, hh)
        lh  <-  dwtyb[[1+3*(l-1)]]; hl  <-  dwtyb[[2+3*(l-1)]];  hh <-  dwtyb[[3+3*(l-1)]]
        env$En2yb[l,e] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        Yb1[ii:jj,e] <- c( lh, hl, hh)
        lh  <-  dwtin[[1+3*(l-1)]]; hl  <-  dwtin[[2+3*(l-1)]];  hh <-  dwtin[[3+3*(l-1)]]
        env$En2in[l,e] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        Yd1[ii:jj,e] <- c( lh, hl, hh)
      }
      ll <- dwtxb[[4+3*(env$n_levs_mx-1)]] 
      env$En2xb[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      Xb1[(jj+1):env$n_dim,e] <- ll
      ll <- dwtyb[[4+3*(env$n_levs_mx-1)]] 
      env$En2yb[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      Yb1[(jj+1):env$n_dim,e] <- ll
      ll <- dwtin[[4+3*(env$n_levs_mx-1)]] 
      env$En2in[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      Yd1[(jj+1):env$n_dim,e] <- ll
      rm( dwtxb, dwtyb, dwtin, lh, hl, hh, ll, ii, jj)

if (plot) {
      ylim<-range(c(env$En2xb[,e],env$En2obs,env$En2yb[,e]))
      yylim<-range(c(env$En2obs,env$En2yb[,e]))
      fout<-file.path(dp,paste0("en2_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
      png(file=fout,width=800,height=800)
      plot(1:env$n_levs_mx,env$En2xb[1:env$n_levs_mx,e],type="l",col="red",ylim=ylim)
      points(env$n_levs_mx,env$En2xb[env$n_levs_mx+1,e],pch=21,bg="red",cex=2)
      lines(1:env$n_levs_mx,env$En2obs[1:env$n_levs_mx],col="blue")
      lines(1:env$n_levs_mx,env$En2obs[1:env$n_levs_mx],col="blue")
      lines(1:env$n_levs_mx,env$En2yb[1:env$n_levs_mx,e],col="pink4")
      points(env$n_levs_mx,env$En2obs[env$n_levs_mx+1],pch=21,bg="blue",cex=2)
      points(env$n_levs_mx,env$En2yb[env$n_levs_mx+1,e],pch=21,bg="pink4",cex=2)
      par(new=T)
      plot(1:env$n_levs_mx,env$En2in[1:env$n_levs_mx],col="gold",axes=F,type="l",lty=1,lwd=2)
      points(env$n_levs_mx,env$En2in[env$n_levs_mx+1],pch=21,bg="gold",cex=2)
      dev.off()
      print(paste("written file",fout))

      br<-c(0,1,2,4,8,16,32,64,128)
      col<-c("lightgray",rev(rainbow(7)))
      fouta<-file.path(dp,paste0("rr1a.png"))
      png(file=fouta,width=1200,height=1200)
      image(rfxb,breaks=br,col=col,main="analysis")
      dev.off()
      foutb<-file.path(dp,paste0("rr1b.png"))
      png(file=foutb,width=1200,height=1200)
      image(rfxb,breaks=br,col=col,main="analysis")
      for (i in 1:length(col)) {
        ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
        points(y_env$yo$x[ix],y_env$yo$y[ix],cex=0.8,pch=21,bg=col[i],col="darkgray")
      }
      dev.off()
      fout<-file.path(dp,paste0("rr1_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
      system(paste0("convert +append ",fouta," ",foutb," ",fxb[e]," ",fout))
      system(paste0("rm -f ",fouta," ",foutb))
      print(paste("written file",fout))

      aux <- extract(rfxb,c_xy[,1])
      ylim=range(c(0,aux,y_env$yo$value))
      fouta<-file.path(dp,paste0("rr1x.png"))
      png(file=fouta,width=1200,height=1200)
      plot(y_env$yo$x,aux,pch=21,bg="lightgray",col="gray")
      points(y_env$yo$x,y_env$yo$value,pch=21,bg="black")
      dev.off()
      foutb<-file.path(dp,paste0("rr1y.png"))
      png(file=foutb,width=1200,height=1200)
      plot(y_env$yo$y,aux,pch=21,bg="lightgray",col="gray")
      points(y_env$yo$y,y_env$yo$value,pch=21,bg="black")
      dev.off()
      foutc<-file.path(dp,paste0("rr1scatt.png"))
      png(file=foutc,width=1200,height=1200)
      plot(y_env$yo$value,aux,pch=21,bg="lightgray",col="gray",xlim=ylim,ylim=ylim)
      dev.off()
      fout<-file.path(dp,paste0("rr1xy_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
      system(paste0("convert +append ",fouta," ",foutb," ",foutc," ",fout))
      system(paste0("rm -f ",fouta," ",foutb," ",foutc))
      print(paste("written file",fout))
}

    } # end loop over ensemble members
    cat( paste0("time=",round(t1-t0,1), attr(t1-t0,"unit"),"\n"))
  } # end main loop
}
