#+
main_wise_analysis_loop_alignment <- function( argv, y_env, fg_env, env, seed=NA, obs_k_dim=1000, 
                                               plot=F, dir_plot=NA) {
#
#------------------------------------------------------------------------------

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

  c_xy <- cellFromXY( rdyad, cbind(y_env$yo$x,y_env$yo$y))
  
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
  env$m_dim <- length(dwt_ix)
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( "state variable, wavelet coefficients, n dim >", env$n_dim, "\n"))
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))

  # ---~--------------
  # -- Observations --

  cat("Transform observations and compute error covariance matrix ")
  
  # rasterize
  # q   <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=rdyad, field=y_env$yo$value)
  # q[is.na(q)] <- 0
  # yo1 <- as.vector( unlist( dwt.2d( as.matrix(q), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)))
  rfobs <- rdyad
  rfobs[c_xy] <- y_env$yo$value
  rfidi <- rdyad
  rfidi[c_xy] <- 1

  # transformation operator
  dwt <- dwt.2d( as.matrix(rfobs), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
  dwt_idi <- dwt.2d( as.matrix(rfidi), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

  # unlist the result and compute squared energies
  yo1 <- vector( mode="numeric", length=env$n_dim); yo1[]<-NA
  yidi1 <- vector( mode="numeric", length=env$n_dim); yidi1[]<-NA
  env$En2idi <- vector( mode="numeric", length=(env$n_levs_mx+1)); env$En2idi[]<-NA
  env$En2obs <- vector( mode="numeric", length=(env$n_levs_mx+1)); env$En2obs[]<-NA
  jj <- 0; ii <- 0
  for (l in 1:env$n_levs_mx) {
    lh <- dwt[[1+3*(l-1)]]; hl <- dwt[[2+3*(l-1)]]; hh <- dwt[[3+3*(l-1)]]
    lh_idi <- dwt_idi[[1+3*(l-1)]]; hl_idi <- dwt_idi[[2+3*(l-1)]]; hh_idi <- dwt_idi[[3+3*(l-1)]]
    ii <- jj + 1
    jj <- ii + length(lh) + length(hl) + length(hh) - 1
    yo1[ii:jj] <- c( lh, hl, hh)
    yidi1[ii:jj] <- c( lh_idi, hl_idi, hh_idi)
    env$En2idi[l] <- mean( ( lh_idi / 2**l)**2) + mean( ( hl_idi / 2**l)**2) + mean( ( hh_idi / 2**l)**2)
    env$En2obs[l] <- mean( ( lh / 2**l)**2)     + mean( ( hl / 2**l)**2)     + mean( ( hh / 2**l)**2)
  }
  ll <- dwt[[4+3*(env$n_levs_mx-1)]] 
  yo1[(jj+1):env$n_dim] <- ll
  env$En2obs[env$n_levs_mx+1] <- mean( ( ll / 2**l)**2)
  ll_idi <- dwt_idi[[4+3*(env$n_levs_mx-1)]] 
  yidi1[(jj+1):env$n_dim] <- ll_idi
  env$En2idi[env$n_levs_mx+1] <- mean( ( ll_idi / 2**env$n_levs_mx)**2)
  rm( lh, hl, hh, ll)
  rm( lh_idi, hl_idi, hh_idi, ll_idi)
  cat("\n")

  # ---~------------
  # -- Background --

  t0 <- Sys.time()
  cat("Transform background")
  Xb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  Xb  <- array( data=NA, dim=c( env$m_dim, env$k_dim))
  Yb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  Yd1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  j <- 0
  env$En2xb <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2yb <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2in <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  for (i in fg_env$ixs) {
    cat(".")
    j <- j + 1

    # interpolate onto the dyadic grid
    # background at grid  points
    rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
    rfxb[rfxb<0] <- 0
    Xb[,j] <- getValues(rfxb)
    # background at observation points
    rfyb  <- rfxb; rfyb[] <- 0; rfyb[c_xy] <- extract( rfxb, c_xy)
    # innovation 
    rfin  <- rfxb; rfin[] <- 0; rfin[c_xy] <- extract( rfobs, c_xy) - extract( rfxb, c_xy)

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

  # Background error covariance matrices
  Ab <- Xb1 - rowMeans(Xb1)
  Db <- Yb1 - rowMeans(Yb1)
#  D  <- Yd1 - rowMeans(Yd1)
  D  <- Yd1
  var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) #* 1/dwt_res**2 
  var_b1_yy <- 1/(env$k_dim-1) * rowSums( Db*Db) #* 1/dwt_res**2
  var_bd1_xy <- 1/(env$k_dim-1) * rowSums( Ab*D) #* 1/dwt_res**2 
  var_d1_yy  <- 1/(env$k_dim-1) * rowSums( D *D) #* 1/dwt_res**2
  rm(Ab,Db,D)

  # Analysis coefficients
#  coeff <- var_b1_xy / var_b1_yy
  coeff <- var_b1_xy / var_d1_yy
  coeff[!is.finite(coeff)] <- 0
  Xa1     <- array( data=NA, dim=c(     env$n_dim, env$k_dim))
  env$Xa  <- array( data=NA, dim=c( length(rdyad), env$k_dim))
  for (e in 1:env$k_dim) {
    Xa1[,e] <- Xb1[,e] + coeff * ( yo1 - Yb1[,e])
    dwt_aux <- dwt_out
    for (i in env$n_levs_mn:env$n_levs_mx) {
      ix <- which( dwt_ix == (i*10 + 1))
      dwt_aux[[3*(i-1)+1]][] <- Xa1[ix,e]
      ix <- which( dwt_ix == (i*10 + 2))
      dwt_aux[[3*(i-1)+2]][] <- Xa1[ix,e]
      ix <- which( dwt_ix == (i*10 + 3))
      dwt_aux[[3*(i-1)+3]][] <- Xa1[ix,e]
    }
    ix <- which( dwt_ix == (env$n_levs_mx*10 + 4))
    dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- Xa1[ix,e]
    env$Xa[,e] <- idwt.2d( dwt_aux)
    rm(dwt_aux)
    env$Xa[,e][env$Xa[,e]<y_env$rain]<-0
q <- rdyad
      q[] <- array( data=env$Xa[,e], dim=c(2**n,2**n))
rfxa <- q
      q[] <- Xb[,e]
rfxb <- q
library(raster)
  br<-c(0,1,2,4,8,16,32,64,128)
  col<-c("lightgray",rev(rainbow(7)))
fouta<-paste0("testa_",e,".png")
    png(file=fouta,width=1200,height=1200)
    image(rfxa,breaks=br,col=col,main="background")
#    for (i in 1:length(col)) {
#      ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
#      points(y_env$yo$x[ix],y_env$yo$y[ix],cex=1,pch=21,bg=col[i])
#    }
dev.off()
foutb<-paste0("testb_",e,".png")
    png(file=foutb,width=1200,height=1200)
    image(rfxb,breaks=br,col=col,main="background")
    for (i in 1:length(col)) {
      ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
      points(y_env$yo$x[ix],y_env$yo$y[ix],cex=0.8,pch=21,bg=col[i],col="darkgray")
    }
dev.off()
fout<-paste0("test_",formatC(e,width=2,flag="0"),".png")
#system(paste0("convert +append ",fouta," ",foutb," test.png ",fout))
system(paste0("convert +append ",fouta," ",foutb," ",fout))
system(paste0("rm -f ",fouta," ",foutb))
print(paste("written file",fout))

  }

q()

###############

  env$En2xb <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2x  <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2yb <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2yo <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2in  <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
ylim<-0
  j <- 0
  for (i in fg_env$ixs) {
    cat(".")
    j <- j + 1

    # -- background at grid  points --

    # interpolate onto the dyadic grid
    rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
    rfxb[rfxb<0] <- 0
    Xs <- array( data=NA, dim=c( env$n_dim, 1))
    Ys <- array( data=NA, dim=c( env$n_dim, 1))
    Xs_ref <- array( data=NA, dim=c( env$n_dim, 1))
    q <- rdyad
  br<-c(0,1,2,4,8,16,32,64,128)
  col<-c("gray",rev(rainbow(7)))
fout<-paste0("test.png")
    png(file=fout,width=1200,height=1200)
    image(rfxb,breaks=br,col=col,main="background")
dev.off()
print(paste("written file",fout))

    for (s in env$n_levs_mx:1) {

#      n_levs_mx <- env$n_levs_mx 
      n_levs_mx <- s

      # transformation operator
      dwtxb  <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=n_levs_mx, boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in 1:n_levs_mx) {
        lh <- dwtxb[[1+3*(l-1)]]; hl <- dwtxb[[2+3*(l-1)]]; hh <- dwtxb[[3+3*(l-1)]]
        env$En2xb[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        ii <- jj + 1
        jj <- ii + length(lh) + length(hl) + length(hh) - 1
        Xs[ii:jj] <- c( lh, hl, hh)
      }
      ll <- dwtxb[[4+3*(n_levs_mx-1)]] 
      env$En2xb[n_levs_mx+1,j] <- mean( (ll/2**n_levs_mx)**2)
      Xs[(jj+1):env$n_dim] <- ll

      # -- background at observation points --

      # rasterize
      # s  <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=rdyad, field=yb, fun=mean)
      rfin  <- rfxb
      rfin[] <- 0
      rfin[c_xy] <- extract( rfobs, c_xy) - extract( rfxb, c_xy)

      rfyb  <- rfxb
      rfyb[] <- 0
      rfyb[c_xy] <- extract( rfxb, c_xy)

      rfx <- rfxb
      rfx[c_xy] <- extract( rfobs, c_xy)

      # transformation
      dwtyin  <- dwt.2d( as.matrix(rfin),  wf=env$wf, J=n_levs_mx, boundary=env$boundary)
      dwtyb   <- dwt.2d( as.matrix(rfyb),  wf=env$wf, J=n_levs_mx, boundary=env$boundary)
      dwtyo   <- dwt.2d( as.matrix(rfobs), wf=env$wf, J=n_levs_mx, boundary=env$boundary)
      dwtx    <- dwt.2d( as.matrix(rfx),   wf=env$wf, J=n_levs_mx, boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in 1:n_levs_mx) {
        lh  <-  dwtyin[[1+3*(l-1)]]; hl  <-  dwtyin[[2+3*(l-1)]];  hh <-  dwtyin[[3+3*(l-1)]]
#        env$En2yb[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        env$En2in[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        ii <- jj + 1
        jj <- ii + length(lh) + length(hl) + length(hh) - 1
        Ys[ii:jj] <- c(  lh,  hl, hh)
        lh  <-  dwtyb[[1+3*(l-1)]]; hl  <-  dwtyb[[2+3*(l-1)]];  hh <-  dwtyb[[3+3*(l-1)]]
        env$En2yb[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        lh  <-  dwtyo[[1+3*(l-1)]]; hl  <-  dwtyo[[2+3*(l-1)]];  hh <-  dwtyo[[3+3*(l-1)]]
        env$En2yo[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        lh  <-  dwtx[[1+3*(l-1)]]; hl  <-  dwtx[[2+3*(l-1)]];  hh <-  dwtx[[3+3*(l-1)]]
        env$En2x[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        Xs_ref[ii:jj] <- c(  lh,  hl, hh)
      }
      ll  <-  dwtyin[[4+3*(n_levs_mx-1)]] 
      env$En2in[n_levs_mx+1,j] <- mean( (ll/2**n_levs_mx)**2)
      Ys[(jj+1):env$n_dim] <- ll
      ll  <-  dwtyb[[4+3*(n_levs_mx-1)]] 
      env$En2yb[n_levs_mx+1,j] <- mean( (ll/2**n_levs_mx)**2)
      ll  <-  dwtyo[[4+3*(n_levs_mx-1)]] 
      env$En2yo[n_levs_mx+1,j] <- mean( (ll/2**n_levs_mx)**2)
      ll  <-  dwtx[[4+3*(n_levs_mx-1)]] 
      env$En2x[n_levs_mx+1,j] <- mean( (ll/2**n_levs_mx)**2)
      Xs_ref[(jj+1):env$n_dim] <- ll

#ylim<-range(c(ylim,env$En2xb[1:n_levs_mx,j],env$En2obs[1:n_levs_mx],env$En2yb[1:n_levs_mx,j]))
ylim<-range(c(ylim,env$En2xb[,j],env$En2x[,j],env$En2obs,env$En2yb[,j]))
yylim<-range(c(env$En2obs,env$En2yb[,j]))
#ylim<-range(c(ylim,env$En2xb[,j],env$En2in[,j]))
#print(env$En2xb[,j])
#print(env$En2in[,j])
fout<-paste0("en2_",formatC(s,width=2,flag="0"),".png")
png(file=fout,width=800,height=800)
plot(1:n_levs_mx,env$En2xb[1:n_levs_mx,j],type="l",col="red",ylim=ylim)
points(n_levs_mx,env$En2xb[n_levs_mx+1,j],pch=21,bg="red",cex=2)
lines(1:n_levs_mx,env$En2x[1:n_levs_mx,j],col="cornflowerblue")
points(n_levs_mx,env$En2x[n_levs_mx+1],pch=21,bg="cornflowerblue",cex=2)
lines(1:n_levs_mx,env$En2obs[1:n_levs_mx],col="blue")
lines(1:n_levs_mx,env$En2obs[1:n_levs_mx],col="blue")
lines(1:n_levs_mx,env$En2yb[1:n_levs_mx,j],col="pink4")
points(n_levs_mx,env$En2obs[n_levs_mx+1],pch=21,bg="blue",cex=2)
points(n_levs_mx,env$En2yb[n_levs_mx+1,j],pch=21,bg="pink4",cex=2)
par(new=T)
plot(1:n_levs_mx,env$En2in[1:n_levs_mx],col="gold",axes=F,type="l",lty=1,lwd=2)
points(n_levs_mx,env$En2in[n_levs_mx+1],pch=21,bg="gold",cex=2)
axis(4)
par(new=T)
plot(1:n_levs_mx,env$En2idi[1:n_levs_mx],col="black",axes=F,type="l",lty=2)
dev.off()
print(paste("written file",fout))

      ij <- ijFromLev( n, n_levs_mx, T) 
#      scaling <- env$En2yb[n_levs_mx+1,j] / env$En2xb[n_levs_mx+1,j]
#      scaling <- env$En2yb[n_levs_mx,j] / env$En2xb[n_levs_mx,j]
#print(scaling)
#print(dwtyin[[3*(n_levs_mx-1)+4]][])
#print(dwtxb[[3*(n_levs_mx-1)+4]][])
        scaling <- 1
#      dwtxb[[3*(n_levs_mx-1)+4]][] <- dwtxb[[3*(n_levs_mx-1)+4]][] + scaling * dwtyin[[3*(n_levs_mx-1)+4]][]
      dwtxb[[3*(n_levs_mx-1)+4]][] <- dwtx[[3*(n_levs_mx-1)+4]][]
      for (r in 1:s) {
#      for (r in s) {
        ij <- ijFromLev( n, r, F)
        scaling <- env$En2xb[r,j] / (env$En2yb[r,j]+env$En2yo[r,j])
        varin <- mean( Ys[ij[1,1]:ij[1,2]]**2)
        varxy <- mean( (Xs[ij[1,1]:ij[1,2]]-Xs_ref[ij[1,1]:ij[1,2]])**2)
#        if (varin>0) {
#          scaling <- (Xs[ij[1,1]:ij[1,2]]-Xs_ref[ij[1,1]:ij[1,2]])**2 / Ys[ij[1,1]:ij[1,2]]**2
#          ix <- which( is.finite(scaling))
          scaling <- 0.5
          dwtxb[[3*(r-1)+1]][] <- dwtxb[[3*(r-1)+1]][] + scaling[] * dwtyin[[3*(r-1)+1]][]
          dwtxb[[3*(r-1)+2]][] <- dwtxb[[3*(r-1)+2]][] + scaling[] * dwtyin[[3*(r-1)+2]][]
          dwtxb[[3*(r-1)+3]][] <- dwtxb[[3*(r-1)+3]][] + scaling[] * dwtyin[[3*(r-1)+3]][]
#        }
      }
      
      Xsn <- idwt.2d( dwtxb)
      q[] <- array( data=Xsn, dim=c(2**n,2**n))
      q[q<0]<-0
      rfxb <- q
library(raster)
  br<-c(0,1,2,4,8,16,32,64,128)
  col<-c("gray",rev(rainbow(7)))
fouta<-paste0("testa_",s,".png")
    png(file=fouta,width=1200,height=1200)
    image(rfxb,breaks=br,col=col,main="background")
#    for (i in 1:length(col)) {
#      ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
#      points(y_env$yo$x[ix],y_env$yo$y[ix],cex=1,pch=21,bg=col[i])
#    }
dev.off()
foutb<-paste0("testb_",s,".png")
    png(file=foutb,width=1200,height=1200)
    image(rfxb,breaks=br,col=col,main="background")
    for (i in 1:length(col)) {
      ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
      points(y_env$yo$x[ix],y_env$yo$y[ix],cex=0.5,pch=21,bg=col[i],col=col[i])
    }
dev.off()
fout<-paste0("test_",formatC(s,width=2,flag="0"),".png")
#system(paste0("convert +append ",fouta," ",foutb," test.png ",fout))
system(paste0("convert +append ",fouta," ",foutb," ",fout))
system(paste0("rm -f ",fouta," ",foutb))
print(paste("written file",fout))
#if (s==(env$n_levs_mx-2)) q()
#      deltaXs <- array( data=0, dim=c( env$n_dim, 1))
#      for (r in 1:s) {
#        ij <- ijFromLev( n, r, F) 
#        deltaXs[ij[1]:ij[2]] <- Xs[ij[1]:ij[2]] + 1/(1.1) * ( yo1[ij[1]:ij[2]] - Ys[ij[1]:ij[2]]) 
#      }
      
#      for (r in 1:s) {
##      aa<-s-1; if (aa==0) aa <-1
##      for (r in aa:s) {
#        ij <- ijFromLev( n, r, F)
#        if (env$En2yb[r,j]>0) {
##          scaling <- env$En2xb[r,j] / env$En2yb[r,j]
#          scaling <- 1
#          dwtxb[[3*(r-1)+1]][] <- dwtxb[[3*(r-1)+1]][] + scaling * 1/(1.1) * ( yo1[ij[1,1]:ij[1,2]] - Ys[ij[1,1]:ij[1,2]])
#          dwtxb[[3*(r-1)+2]][] <- dwtxb[[3*(r-1)+1]][] + scaling * 1/(1.1) * ( yo1[ij[2,1]:ij[2,2]] - Ys[ij[2,1]:ij[2,2]])
#          dwtxb[[3*(r-1)+3]][] <- dwtxb[[3*(r-1)+1]][] + scaling * 1/(1.1) * ( yo1[ij[3,1]:ij[3,2]] - Ys[ij[3,1]:ij[3,2]])
#        }
#      }
##      if ( s == n_levs_mx) {
##        ij <- ijFromLev( n, r, T) 
##        if (env$En2yb[n_levs_mx+1,j]>0) {
###          scaling <- env$En2xb[n_levs_mx+1,j] / env$En2yb[n_levs_mx+1,j]
##          scaling <- 1
##          dwtxb[[3*(n_levs_mx-1)+4]][] <- dwtxb[[3*(n_levs_mx-1)+4]][] + scaling * 1/(1.1) * ( yo1[ij[1]:ij[2]] - Ys[ij[1]:ij[2]])
##        }
##      }
#      Xsn <- idwt.2d( dwtxb)
#      q[] <- array( data=Xsn, dim=c(2**n,2**n))
#      q[q<0]<-0
#      rfxb <- q
#library(raster)
#  br<-c(0,1,2,4,8,16,32,64,128)
#  col<-c("gray",rev(rainbow(7)))
#fouta<-paste0("testa_",s,".png")
#    png(file=fouta,width=1200,height=1200)
#    image(rfxb,breaks=br,col=col,main="background")
##    for (i in 1:length(col)) {
##      ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
##      points(y_env$yo$x[ix],y_env$yo$y[ix],cex=1,pch=21,bg=col[i])
##    }
#dev.off()
#foutb<-paste0("testb_",s,".png")
#    png(file=foutb,width=1200,height=1200)
#    image(rfxb,breaks=br,col=col,main="background")
#    for (i in 1:length(col)) {
#      ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
#      points(y_env$yo$x[ix],y_env$yo$y[ix],cex=1,pch=21,bg=col[i])
#    }
#dev.off()
#fout<-paste0("test_",formatC(s,width=2,flag="0"),".png")
##system(paste0("convert +append ",fouta," ",foutb," test.png ",fout))
#system(paste0("convert +append ",fouta," ",foutb," ",fout))
#system(paste0("rm -f ",fouta," ",foutb))
#print(paste("written file",fout))
    } # end loop over spatial scales
q()

  } # end loop over ensemble members

  # ---~------------
  # -- Background --

  t0 <- Sys.time()
  cat("Transform background")
  Xb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  Yb1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  Yd1 <- array( data=NA, dim=c( env$n_dim, env$k_dim))
  j <- 0
  env$En2b <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  env$En2i <- array( data=NA, dim=c(env$n_levs_mx+1,nfg))
  for (i in fg_env$ixs) {
    cat(".")
    j <- j + 1

    # -- background at grid  points --

    # interpolate onto the dyadic grid
    s <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
    s[s<0] <- 0

    # transformation operator
    dwt  <- dwt.2d( as.matrix(s), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

    # unlist the result and compute squared energies
    jj <- 0; ii <- 0
    for (l in 1:env$n_levs_mx) {
      lh <- dwt[[1+3*(l-1)]]; hl <- dwt[[2+3*(l-1)]]; hh <- dwt[[3+3*(l-1)]]
      env$En2b[l,j] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
      ii <- jj + 1
      jj <- ii + length(lh) + length(hl) + length(hh) - 1
      Xb1[ii:jj,j] <- c( lh, hl, hh)
    }
    ll <- dwt[[4+3*(env$n_levs_mx-1)]] 
    env$En2b[env$n_levs_mx+1,j] <- mean( (ll/2**env$n_levs_mx)**2)
    Xb1[(jj+1):env$n_dim,j] <- ll
    rm( dwt, lh, hl, hh, ll, ii, jj)

    # -- background at observation points --

    # rasterize
    # s  <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=rdyad, field=yb, fun=mean)
    sb  <- s
    s[] <- 0
    s[c_xy] <- extract( sb, c_xy)
    rm(sb)

    # transformation
    dwt  <- dwt.2d( as.matrix(s), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
    dwti <- dwt.2d( as.matrix(so-s), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

    # unlist the result and compute squared energies
    jj <- 0; ii <- 0
    for (l in 1:env$n_levs_mx) {
      lh  <-  dwt[[1+3*(l-1)]]; hl  <-  dwt[[2+3*(l-1)]];  hh <-  dwt[[3+3*(l-1)]]
      lhi <- dwti[[1+3*(l-1)]]; hli <- dwti[[2+3*(l-1)]]; hhi <- dwti[[3+3*(l-1)]]
      env$En2i[l,j] <- mean( (lhi / 2**l)**2) + mean( (hli / 2**l)**2) + mean( (hhi / 2**l)**2)
      ii <- jj + 1
      jj <- ii + length(lh) + length(hl) + length(hh) - 1
      Yb1[ii:jj,j] <- c(  lh,  hl, hh)
      Yd1[ii:jj,j] <- c( lhi, hli, hhi)
    }
    lli <- dwti[[4+3*(env$n_levs_mx-1)]] 
    env$En2i[env$n_levs_mx+1,j] <- mean( (lli/2**env$n_levs_mx)**2)
    Yd1[(jj+1):env$n_dim,j] <- lli
    ll  <-  dwt[[4+3*(env$n_levs_mx-1)]] 
    Yb1[(jj+1):env$n_dim,j] <- ll
    rm( lh, hl, hh, ll, dwt, dwti, ii, jj)

  } # end loop over background fields
  t1 <- Sys.time()
  cat( paste0("time=",round(t1-t0,1), attr(t1-t0,"unit"),"\n"))
  rm(s)

  # Background error covariance matrices
  Ab <- Xb1 - rowMeans(Xb1)
  Db <- Yb1 - rowMeans(Yb1)
  D  <- Yd1 - rowMeans(Yd1)
  var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) #* 1/dwt_res**2 
  var_b1_yy <- 1/(env$k_dim-1) * rowSums( Db*Db) #* 1/dwt_res**2
  var_bd1_xy <- 1/(env$k_dim-1) * rowSums( Ab*D) #* 1/dwt_res**2 
  var_d1_yy  <- 1/(env$k_dim-1) * rowSums( D *D) #* 1/dwt_res**2
  rm(Ab,Db,D)

  # Analysis
  cat("Analysis on the transformed space ")
#  var_o1 <- 0.1 * var_b1_yy
  var_o1 <- 0
  coeff <- var_b1_xy / ( var_b1_yy + var_o1)
  coeff[!is.finite(coeff)] <- 0
  Xa1     <- array( data=NA, dim=c(     env$n_dim, env$k_dim))
  env$Xa  <- array( data=NA, dim=c( length(rdyad), env$k_dim))
  for (e in 1:env$k_dim) {
    Xa1[]
    for (l in env$n_levs_mx:1) {

    }
    cat(".")
    Xa1[,e] <- Xb1[,e] + coeff * ( yo1 - Yb1[,e])
    dwt_aux <- dwt_out
    for (i in env$n_levs_mn:env$n_levs_mx) {
      ix <- which( dwt_ix == (i*10 + 1))
      dwt_aux[[3*(i-1)+1]][] <- Xa1[ix,e]
      ix <- which( dwt_ix == (i*10 + 2))
      dwt_aux[[3*(i-1)+2]][] <- Xa1[ix,e]
      ix <- which( dwt_ix == (i*10 + 3))
      dwt_aux[[3*(i-1)+3]][] <- Xa1[ix,e]
    }
    ix <- which( dwt_ix == (env$n_levs_mx*10 + 4))
    dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- Xa1[ix,e]
    env$Xa[,e] <- idwt.2d( dwt_aux)
    rm(dwt_aux)
    env$Xa[,e][env$Xa[,e]<y_env$rain]<-0
  }
  cat("\n")
}
