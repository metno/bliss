#+
main_wise_analysis <- function( argv, y_env, fg_env, env, seed=NA, obs_k_dim=1000, 
                                plot=F, dir_plot=NA) {
#
#------------------------------------------------------------------------------
  
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
  c_xy <- cellFromXY( rdyad, cbind(y_env$yo$x,y_env$yo$y))
  
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

  # define constants
  # total number of coefficients used. Dimension of the state variable
  env$m_dim <- length(dwt_ix)
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( "state variable, dim >", env$m_dim, "\n"))

  # ---~--------------
  # -- Observations --

  cat("Transform observations and compute error covariance matrix ")
  
  # rasterize
  # q   <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=rdyad, field=y_env$yo$value)
  # q[is.na(q)] <- 0
  # yo1 <- as.vector( unlist( dwt.2d( as.matrix(q), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)))
  so <- rdyad
  so[c_xy] <- y_env$yo$value
  sidi <- rdyad
  sidi[c_xy] <- 1

  # transformation operator
  dwt <- dwt.2d( as.matrix(so), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
  dwt_idi <- dwt.2d( as.matrix(sidi), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

  # unlist the result and compute squared energies
  yo1 <- vector( mode="numeric", length=env$m_dim); yo1[]<-NA
  yidi1 <- vector( mode="numeric", length=env$m_dim); yidi1[]<-NA
  env$En2idi <- vector( mode="numeric", length=env$m_dim); env$En2idi[]<-NA
  jj <- 0; ii <- 0
  for (l in 1:env$n_levs_mx) {
    lh <- dwt[[1+3*(l-1)]]; hl <- dwt[[2+3*(l-1)]]; hh <- dwt[[3+3*(l-1)]]
    lh_idi <- dwt_idi[[1+3*(l-1)]]; hl_idi <- dwt_idi[[2+3*(l-1)]]; hh_idi <- dwt_idi[[3+3*(l-1)]]
    ii <- jj + 1
    jj <- ii + length(lh) + length(hl) + length(hh) - 1
    yo1[ii:jj] <- c( lh, hl, hh)
    yidi1[ii:jj] <- c( lh_idi, hl_idi, hh_idi)
    env$En2idi[l] <- mean( (lh_idi / 2**l)**2) + mean( (hl_idi / 2**l)**2) + mean( (hh_idi / 2**l)**2)
  }
  ll <- dwt[[4+3*(env$n_levs_mx-1)]] 
  yo1[(jj+1):env$m_dim] <- ll
  ll_idi <- dwt_idi[[4+3*(env$n_levs_mx-1)]] 
  yidi1[(jj+1):env$m_dim] <- ll_idi
  env$En2idi[env$n_levs_mx+1] <- mean( (ll_idi/2**env$n_levs_mx)**2)
  rm( lh, hl, hh, ll)
  rm( lh_idi, hl_idi, hh_idi, ll_idi)
#  if (plot) {
#    yo_check <- extract( q, c_xy)
#    fout<-file.path(dir_plot,paste0("wise_checkobs.png"))
#    fx<-file.path(dir_plot,paste0("wise_checkobs_x.png"))
#    png(file=fx,width=1200,height=1200)
#    plot(y_env$yo$x,yo_check,ylim=range(c(y_env$yo$value,yo_check)),pch=21,bg="pink",col="red")
#    points(y_env$yo$x,y_env$yo$value,lwd=2,col="gold")
#    abline(h=0.1)
#    abline(h=1:100,lty=2,col="gray")
#    abline(h=y_env$rain,col="green")
#    dev.off()
#    fy<-file.path(dir_plot,paste0("wise_xa_y_e",formatC(n,width=2,flag="0"),".png"))
#    png(file=fy,width=1200,height=1200)
#    plot(y_env$yo$y,yo_check,ylim=range(c(y_env$yo$value,yo_check)),pch=21,bg="cornflowerblue",col="blue")
#    points(y_env$yo$y,y_env$yo$value,lwd=2,col="gold")
#    dev.off()
#    system( paste0("convert +append ",fx," ",fy," ",fout))
#    system( paste0("rm ",fx," ",fy))
#  }

  # Observation error covariance matrix is assumed to be diagonal in multi-resolution space
  #  based on an ensemble of observation error vectors
#  env$var_o <- y_env$yo$value * 0.01
#  epso_dwt <- array(data=NA,dim=c(env$m_dim,obs_k_dim))
#  if (!is.na(seed)) set.seed(seed)
#  for (i in 1:obs_k_dim) {
#    cat(".")
#    q   <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=rdyad, 
#                      field=rnorm( env$p_dim, mean=0, sd=sqrt(env$var_o)))
#    q[is.na(q)] <- 0
#    epso_dwt[,i] <- as.vector( unlist( dwt.2d( as.matrix(q), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)))
#  }
  cat("\n")
#  Do <- epso_dwt - rowMeans(epso_dwt)
#  var_o1 <- 1/(obs_k_dim-1) * rowSums( Do*Do) * 1/dwt_res**2
#  rm(q,Do,epso_dwt)

  # ---~------------
  # -- Background --

  t0 <- Sys.time()
  cat("Transform background")
  Xb1 <- array( data=NA, dim=c( env$m_dim, env$k_dim))
  Yb1 <- array( data=NA, dim=c( env$m_dim, env$k_dim))
  Yd1 <- array( data=NA, dim=c( env$m_dim, env$k_dim))
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

#    if (plot) {
#      br<-c(0,1,2,4,8,16,32,64,128)
#      col<-c("gray",rev(rainbow(7)))
#      fout<-file.path(dir_plot,paste0("wise_main_xb_e",formatC(j,width=2,flag="0"),".png"))
#      png(file=fout,width=1200,height=1200)
#      image(s,breaks=br,col=col,main="background")
#      plot(b,add=T)
#      for (i in 1:length(col)) {
#        ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
#        cex<-1; col_i <- col[i]
#        if (i==1) { cex<-0.5; col_i<-"beige" }
#        points(y_env$yo$x[ix],y_env$yo$y[ix],pch=21,bg=col_i,col=col_i,cex=cex)
#      }
#      dev.off()
#    }

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
    Xb1[(jj+1):env$m_dim,j] <- ll
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
    Yd1[(jj+1):env$m_dim,j] <- lli
    ll  <-  dwt[[4+3*(env$n_levs_mx-1)]] 
    Yb1[(jj+1):env$m_dim,j] <- ll
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
#  coeff_res <- (sqrt(env$En2i[,1]) / sqrt(env$En2b[,1])) / max(sqrt(env$En2i[,1]) / sqrt(env$En2b[,1]))
#  coeff <- vector( mode="numeric", length=env$m_dim)
#  for (l in 1:env$n_levs_mx) {
#    ix <- which( dwt_res == (2**l))
#    coeff[ix] <- coeff_res[l]
#  }
  coeff_res <- vector( mode="numeric", length=env$n_levs_mx)
  for (l in 1:env$n_levs_mx) {
    ix <- which( dwt_res == (2**l))
    coeff_res[l] <- mean(coeff[ix])
  }
  Xa1     <- array( data=NA, dim=c(     env$m_dim, env$k_dim))
  env$Xa  <- array( data=NA, dim=c( length(rdyad), env$k_dim))
  for (e in 1:env$k_dim) {
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
#    if (plot) {
#      q<-rdyad
#      q[]<-array(data=env$Xa[,e],dim=c(sqrt(length(env$Xa[,e])),sqrt(length(env$Xa[,e]))))
#      fout<-file.path(dir_plot,paste0("wise_main_xcombo_e",formatC(e,width=2,flag="0"),".png"))
#      fb<-file.path(dir_plot,paste0("wise_main_xb_e",formatC(e,width=2,flag="0"),".png"))
#      fa<-file.path(dir_plot,paste0("wise_main_xa_e",formatC(e,width=2,flag="0"),".png"))
#      png(file=fa,width=1200,height=1200)
#      image(q,breaks=br,col=col,main="analysis")
#      plot(b,add=T)
#      for (i in 1:length(col)) {
#        ix <- which(yo_check>=br[i] & yo_check<br[i+1])
#        cex<-1; col_i <- col[i]
#        if (i==1) { cex<-0.5; col_i<-"beige" }
#        points(y_env$yo$x[ix],y_env$yo$y[ix],pch=21,bg=col_i,col=col_i,cex=cex)
#      }
#      dev.off()
#      system( paste0("convert +append ",fb," ",fa," ",fout))
#      system( paste0("rm ",fa," ",fb))
#
#      ya <- extract( q, c_xy)
#      i<-fg_env$ixs[e]
#      s <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
#      yb <- extract( s, c_xy) 
#      yo <- yo_check  
#      fout<-file.path(dir_plot,paste0("wise_main_ycombo_e",formatC(e,width=2,flag="0"),".png"))
#      fb<-file.path(dir_plot,paste0("wise_yb_e",formatC(e,width=2,flag="0"),".png"))
#      fa<-file.path(dir_plot,paste0("wise_ya_e",formatC(e,width=2,flag="0"),".png"))
#      png(file=fb,width=1200,height=1200)
#      plot(yo,yb,ylim=range(c(yo,ya,yb)),xlim=range(c(yo,ya,yb)),main="background vs observations",pch=21,bg="cornflowerblue",col="blue")
#      lines(-10000:10000,-10000:10000)
#      seq<-c(0,0.1,seq(0.5,100,by=0.5))
#      for (ii in 1:(length(seq)-1)) {
#        ix <- which( yo>=seq[ii] & yo<seq[ii+1])
#        if (length(ix)<10) next
#        points( (seq[ii]+seq[ii+1])/2, median(yb[ix]), cex=2, pch=21,bg="black")
#        points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(yb[ix],probs=0.25)), cex=1, pch=21,bg="black")
#        points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(yb[ix],probs=0.75)), cex=1, pch=21,bg="black")
#      }
#      dev.off()
#      png(file=fa,width=1200,height=1200)
#      plot(yo,ya,ylim=range(c(yo,ya,yb)),xlim=range(c(yo,ya,yb)),main="analysis vs observations",pch=21,bg="pink",col="red")
#      lines(-10000:10000,-10000:10000)
#      for (ii in 1:(length(seq)-1)) {
#        ix <- which( yo>=seq[ii] & yo<seq[ii+1])
#        if (length(ix)<10) next
#        points( (seq[ii]+seq[ii+1])/2, median(ya[ix]), cex=2, pch=21,bg="black")
#        points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(ya[ix],probs=0.25)), cex=1, pch=21,bg="black")
#        points( (seq[ii]+seq[ii+1])/2, as.numeric(quantile(ya[ix],probs=0.75)), cex=1, pch=21,bg="black")
#      }
#      dev.off()
#      system( paste0("convert +append ",fb," ",fa," ",fout))
#      system( paste0("rm ",fa," ",fb))
#
#      rm(q)
#    }

  }
  cat("\n")
}
