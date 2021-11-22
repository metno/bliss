#+
main_wise_analysis_loop <- function( argv, y_env, fg_env, env, seed=NA, obs_k_dim=1000, 
                                plot=F, dir_plot=NA) {
#
#------------------------------------------------------------------------------

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
  rfxb <- rdyad
  rfyb <- rdyad
  rfin <- rdyad
  rfxb_alt <- rdyad

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))
  
  # ---~--------------
  # Superobbing

#   uncomment the following lines to use radar
#  radarval <- getValues( u_env$uo[[1]]$r_main)
#  ixrad <- which(!is.na(radarval))
#   
#  y_env$yo$x <- c( y_env$yo$x, xyFromCell(u_env$uo[[1]]$r_main,ixrad)[,1])
#  y_env$yo$y <- c( y_env$yo$y, xyFromCell(u_env$uo[[1]]$r_main,ixrad)[,2])
#  y_env$yo$value <- c( y_env$yo$value, radarval[ixrad])

  t0 <- Sys.time()
  rfobs<-rdyad
  nn2 <- nn2( cbind(y_env$yo$x, y_env$yo$y), 
              query = xyFromCell( rdyad, 1:ncell(rdyad)), 
              k = 50, 
              searchtype = "radius", radius = 1500)

  mat <- nn2[[1]]
  c_xy <- which( ( aux <- rowSums( mat)) > 0 )

  mapply_quantile  <- function(i) { quantile( y_env$yo$value[mat[c_xy[i],1:length(which(mat[c_xy[i],]!=0))]], probs=0.99) }
  rfobs[] <- 0
  rfobs[c_xy] <- t( mapply( mapply_quantile, 1:length(c_xy), SIMPLIFY = T))
  t1 <- Sys.time()
  print(paste("mapply",t1-t0))

  # ---~--------------
  # Initialization of the wavelet structures
  #  dwt_out. dwt.2d class.
  #    for the i-th level (i=1,...,env$n_levs_mx) 
  #      dwt_out[[3*(i-1)+1]] LHi - wavelet coefficients
  #      dwt_out[[3*(i-1)+2]] HLi - wavelet coefficients
  #      dwt_out[[3*(i-1)+3]] HHi - wavelet coefficients
  #  resolution of the dyadic tree. coefficients = 2**i; base = 2**(i-1)
  #  resolution is the number of original grid points (in each direction) within the box a dydadic tree at level i
  # then i=1 is the finest resolution and i=n_levs_mx the coarser
  dwt_out <- dwt.2d( as.matrix(rdyad), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
  for (i in 1:length(dwt_out)) dwt_out[[i]][] <- 0

  # ---~--------------
  # define constants
  # total number of coefficients used. Dimension of the state variable
  env$n_dim <- length( as.vector( unlist( dwt_out)))
  # total number of coefficients used. Dimension of the state variable
  env$m_dim <- length( getValues(rdyad))
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( "state variable, wavelet coefficients, n dim >", env$n_dim, "\n"))
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))
  cat( paste( "number of observations, superob >", length(c_xy), "\n"))

  # ---~--------------
  # -- Observations --

  cat("Transform observations and compute error covariance matrix ")
 
  # transformation operator
  dwt     <- dwt.2d( as.matrix(rfobs), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

  # unlist the result and compute squared energies
  vo   <- vector( mode="numeric", length=env$n_dim); vo[]<-NA
  En2vo <- vector( mode="numeric", length=(env$n_levs_mx+1)); En2vo[]<-NA
  jj <- 0; ii <- 0
  for (l in 1:env$n_levs_mx) {
    lh <- dwt[[1+3*(l-1)]]; hl <- dwt[[2+3*(l-1)]]; hh <- dwt[[3+3*(l-1)]]
    ii <- jj + 1
    jj <- ii + length(lh) + length(hl) + length(hh) - 1
    vo[ii:jj] <- c( lh, hl, hh)
    En2vo[l] <- mean( ( lh / 2**l)**2) + mean( ( hl / 2**l)**2) + mean( ( hh / 2**l)**2)
  }
  ll <- dwt[[4+3*(env$n_levs_mx-1)]] 
  vo[(jj+1):env$n_dim] <- ll
  En2vo[env$n_levs_mx+1] <- mean( ( ll / 2**l)**2)
  rm( lh, hl, hh, ll)
  cat("\n")

  #--------------------------------------------------------    
  # Analysis
  cat("Analysis on the transformed space \n")

  costf <- vector()
#  ya <- extract( rfxb, cbind( y_env$yo$x, y_env$yo$y))

  En2 <- array( data=NA, dim=c(env$n_levs_mx+1,env$k_dim))
  En2_prev <- array( data=NA, dim=c(env$n_levs_mx+1,env$k_dim))
  rho <- array( data=NA, dim=c(env$n_levs_mx+1,env$k_dim))

#  En2_prev <- En2Ub
#  aux <- En2_prev - rowMeans(En2_prev) 
#  var_En2_prev <- 1/(env$k_dim-1) * rowSums( aux * aux)
#  rm(aux)

  #--------------------------------------------------------    
  # MAIN LOOP
  for (loop in 1:20) {

    t0 <- Sys.time()
    cat( paste(" loop", loop, " "))

    #--------------------------------------------------------    
    # set the background

    Ub     <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    Ub_alt <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    Vb     <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    vo_Vb  <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    En2Ub     <- array( data=NA, dim=c( env$n_levs_mx+1, env$k_dim))
    En2Ub_alt <- array( data=NA, dim=c( env$n_levs_mx+1, env$k_dim))
    En2Vb     <- array( data=NA, dim=c( env$n_levs_mx+1, env$k_dim))
    En2in     <- array( data=NA, dim=c( env$n_levs_mx+1, env$k_dim))

    if (plot & loop == 1) fxb <- vector()
    
    for (e in 1:env$k_dim) {

      # first iteration - background from ensemble 
      if ( loop == 1) {
        i <- fg_env$ixs[e]
        # interpolate onto the dyadic grid
        # background at grid  points
        rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
        rfxb[rfxb<0] <- 0

      # second iteration onwards - background from previous iteration 
      } else {
        # -- background at grid  points --
  #      rfxb_to_plot[] <- array(data=Xa[,e],dim=c(sqrt(length(Xa[,e])),sqrt(length(Xa[,e]))))
        rfxb[] <- array(data=Xa[,e],dim=c(sqrt(length(Xa[,e])),sqrt(length(Xa[,e]))))
        rfxb[rfxb<0] <- 0
#        ya <- extract( rfxb, cbind( y_env$yo$x, y_env$yo$y))
      }

      # background at observation points
      rfyb <- rdyad; rfyb[c_xy] <- rfxb[c_xy]
      # innovation 
      rfin <- rdyad; rfin[c_xy] <- rfobs[c_xy] - rfxb[c_xy]
      # background alternative (obs+backg)
      rfxb_alt <- rfxb; rfxb_alt[c_xy] <- rfobs[c_xy]

      # plot the background
      if (plot & loop==1) {
        fxb[e] <- file.path(dir_plot,paste0("rfxb_",formatC(e,flag="0",width=2),".png"))
        br<-c(0,1,2,4,8,16,32,64,128)
        col<-c("lightgray",rev(rainbow(7)))
        png(file=fxb[e],width=1200,height=1200)
        image(rfxb,breaks=br,col=col,main="background")
        dev.off()
      }

      # transformation operator
      dwtxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
      dwtxb_alt <- dwt.2d( as.matrix(rfxb_alt), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
      dwtyb <- dwt.2d( as.matrix(rfyb), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
      dwtin <- dwt.2d( as.matrix(rfin), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in env$n_levs_mn:env$n_levs_mx) {
        lh <- dwtxb[[1+3*(l-1)]]; hl <- dwtxb[[2+3*(l-1)]]; hh <- dwtxb[[3+3*(l-1)]]
        En2Ub[l,e] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        ii <- jj + 1
        jj <- ii + length(lh) + length(hl) + length(hh) - 1
        Ub[ii:jj,e] <- c( lh, hl, hh)
        lh  <-  dwtxb_alt[[1+3*(l-1)]]; hl  <-  dwtxb_alt[[2+3*(l-1)]];  hh <-  dwtxb_alt[[3+3*(l-1)]]
        En2Ub_alt[l,e] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        Ub_alt[ii:jj,e] <- c( lh, hl, hh)
        lh  <-  dwtyb[[1+3*(l-1)]]; hl  <-  dwtyb[[2+3*(l-1)]];  hh <-  dwtyb[[3+3*(l-1)]]
        En2Vb[l,e] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        Vb[ii:jj,e] <- c( lh, hl, hh)
        lh  <-  dwtin[[1+3*(l-1)]]; hl  <-  dwtin[[2+3*(l-1)]];  hh <-  dwtin[[3+3*(l-1)]]
        En2in[l,e] <- mean( (lh / 2**l)**2) + mean( (hl / 2**l)**2) + mean( (hh / 2**l)**2)
        vo_Vb[ii:jj,e] <- c( lh, hl, hh)
      }
      ll <- dwtxb[[4+3*(env$n_levs_mx-1)]] 
      En2Ub[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      Ub[(jj+1):env$n_dim,e] <- ll
      ll <- dwtxb[[4+3*(env$n_levs_mx-1)]] 
      En2Ub_alt[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      Ub_alt[(jj+1):env$n_dim,e] <- ll
      ll <- dwtyb[[4+3*(env$n_levs_mx-1)]] 
      En2Vb[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      Vb[(jj+1):env$n_dim,e] <- ll
      ll <- dwtin[[4+3*(env$n_levs_mx-1)]] 
      En2in[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      vo_Vb[(jj+1):env$n_dim,e] <- ll

      rm( dwtxb, dwtyb, dwtin, lh, hl, hh, ll, ii, jj)

    } # end loop over ensemble members

#    En2_prev <- En2Ub_alt
    En2_prev <- En2Ub
    aux <- En2_prev - rowMeans(En2_prev) 
    var_En2_prev <- 1/(env$k_dim-1) * rowSums( aux * aux)
    rm(aux)
    
    # set the background: end
    #--------------------------------------------------------    

    #--------------------------------------------------------    
    # Background error covariance matrices
    Ab <- Ub - rowMeans(Ub)
    Db <- Vb - rowMeans(Vb)
    D  <- vo_Vb
    var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) #* 1/dwt_res**2 
    var_b1_yy <- 1/(env$k_dim-1) * rowSums( Db*Db) #* 1/dwt_res**2
    var_bd1_xy <- 1/(env$k_dim-1) * rowSums( Ab*D) #* 1/dwt_res**2 
    var_d1_yy  <- 1/(env$k_dim-1) * rowSums( D *D) #* 1/dwt_res**2

    D_alt  <- vo_Vb - rowMeans(vo_Vb)
    Ab_alt <- Ub_alt - rowMeans(Ub_alt)
    var_bd1alt_xy <- 1/(env$k_dim-1) * rowSums( Ab_alt*D_alt)
    var_d1alt_yy  <- 1/(env$k_dim-1) * rowSums( D_alt *D_alt)

    rm( Ab, Db, D, D_alt, Ab_alt)

    #--------------------------------------------------------    
    #   
    costf[loop] <- mean( sqrt( var_d1_yy))
    print( paste("rmse", costf[loop]))

    #--------------------------------------------------------    
    # Analysis
    
    coeff <- var_b1_xy / var_d1_yy
#    coeff <- var_bd1alt_xy / var_d1alt_yy
    coeff[!is.finite(coeff)] <- 0
    Ua  <- array( data=NA, dim=c(     env$n_dim, env$k_dim))
    Xa  <- array( data=NA, dim=c( length(rdyad), env$k_dim))
    for (e in 1:env$k_dim) {
      cat(".")
      Ua[,e] <- Ub[,e] + coeff * ( vo - Vb[,e])
      dwt_aux <- dwt_out
      for (i in env$n_levs_mn:env$n_levs_mx) {
        ij <- ijFromLev( env$n_levs_mx, i, F)
        dwt_aux[[3*(i-1)+1]][] <- Ua[ij[1,1]:ij[1,2],e]
        dwt_aux[[3*(i-1)+2]][] <- Ua[ij[2,1]:ij[2,2],e]
        dwt_aux[[3*(i-1)+3]][] <- Ua[ij[3,1]:ij[3,2],e]
      }
      ij <- ijFromLev( env$n_levs_mx, env$n_levs_mx, T)
      dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- Ua[ij[1]:ij[2],e]
      # reconstruct analysis
      Xa[,e] <- idwt.2d( dwt_aux)
      Xa[,e][Xa[,e]<y_env$rain]<-0

      # compute energies
      rfxb[] <- array(data=Xa[,e],dim=c(sqrt(length(Xa[,e])),sqrt(length(Xa[,e]))))
      dwtxa <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
      jj <- 0; ii <- 0
      for (i in env$n_levs_mn:env$n_levs_mx) {
        lh <- dwtxa[[1+3*(i-1)]]
        hl <- dwtxa[[2+3*(i-1)]]
        hh <- dwtxa[[3+3*(i-1)]]
        En2[i,e] <- mean( (lh / 2**i)**2) + 
                    mean( (hl / 2**i)**2) + 
                    mean( (hh / 2**i)**2)
        ii <- jj + 1
        jj <- ii + length(lh) + length(hl) + length(hh) - 1
        Ua[ii:jj,e] <- c( lh, hl, hh)
      }
      ll <- dwtxa[[4+3*(env$n_levs_mx-1)]] 
      En2[env$n_levs_mx+1,e] <- mean( (ll/2**env$n_levs_mx)**2)
      Ua[(jj+1):env$n_dim,e] <- ll
    }
    t1 <- Sys.time()
    cat( paste0("time=",round(t1-t0,1), attr(t1-t0,"unit"),"\\"))

    # Constraint on energies
    for (e in 1:env$k_dim) {
      rho[,e] <- exp( -0.5 * (En2[,e] - En2_prev[,e])**2 / var_En2_prev)
#      Ua[,e] <- Ub[,e] + coeff * ( vo - Vb[,e])
      dwt_aux <- dwt_out
      for (i in env$n_levs_mn:env$n_levs_mx) {
        ij <- ijFromLev( env$n_levs_mx, i, F)
        dwt_aux[[3*(i-1)+1]][] <- Ub[ij[1,1]:ij[1,2],e] + rho[i,e] * coeff[ij[1,1]:ij[1,2]] * ( vo[ij[1,1]:ij[1,2]] - Vb[ij[1,1]:ij[1,2],e])
        dwt_aux[[3*(i-1)+2]][] <- Ub[ij[2,1]:ij[2,2],e] + rho[i,e] * coeff[ij[2,1]:ij[2,2]] * ( vo[ij[2,1]:ij[2,2]] - Vb[ij[2,1]:ij[2,2],e])
        dwt_aux[[3*(i-1)+3]][] <- Ub[ij[3,1]:ij[3,2],e] + rho[i,e] * coeff[ij[3,1]:ij[3,2]] * ( vo[ij[3,1]:ij[3,2]] - Vb[ij[3,1]:ij[3,2],e])
#        dwt_aux[[3*(i-1)+1]][] <- rho[i,e] * Ua[ij[1,1]:ij[1,2],e]
#        dwt_aux[[3*(i-1)+2]][] <- rho[i,e] * Ua[ij[2,1]:ij[2,2],e]
#        dwt_aux[[3*(i-1)+3]][] <- rho[i,e] * Ua[ij[3,1]:ij[3,2],e]
        En2[i,e] <- mean( ( dwt_aux[[3*(i-1)+1]][] / 2**i)**2) + 
                    mean( ( dwt_aux[[3*(i-1)+2]][] / 2**i)**2) +
                    mean( ( dwt_aux[[3*(i-1)+3]][] / 2**i)**2)
      }
      ij <- ijFromLev( env$n_levs_mx, env$n_levs_mx, T)
      dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- Ub[ij[1]:ij[2],e] + rho[i,e] * coeff[ij[1]:ij[2]] * ( vo[ij[1]:ij[2]] - Vb[ij[1]:ij[2],e])
#      dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- rho[env$n_levs_mx+1,e] * Ua[ij[1]:ij[2],e]
      En2[env$n_levs_mx+1,e] <- mean( ( dwt_aux[[3*(env$n_levs_mx-1)+4]][] / 2**i)**2)
      Xa[,e] <- idwt.2d( dwt_aux)
      Xa[,e][Xa[,e]<y_env$rain]<-0
    }
 
    if (plot) {
save(file="tmp.rdata",env,En2,var_En2_prev,En2_prev)
      for (e in 1:env$k_dim) {
        rfxb[] <- array(data=Xa[,e],dim=c(sqrt(length(Xa[,e])),sqrt(length(Xa[,e]))))
        if (loop==1) ylim<-range(En2Ub)
        fout<-file.path(dir_plot,paste0("en2_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
        png(file=fout,width=800,height=800)
#        plot(1:env$n_levs_mx,En2Ub[1:env$n_levs_mx,e],type="l",col="red",ylim=ylim)
        plot(1:env$n_levs_mx,En2[1:env$n_levs_mx,e],type="l",col="red",ylim=ylim)
        points(env$n_levs_mx,En2[env$n_levs_mx+1,e],pch=21,bg="red",cex=2)
        lines(1:env$n_levs_mx,En2_prev[1:env$n_levs_mx,e],col="pink")
        points(env$n_levs_mx,En2_prev[env$n_levs_mx+1,e],pch=21,bg="pink",cex=2)
        polygon(c(1:env$n_levs_mx,env$n_levs_mx:1),c(En2_prev[1:env$n_levs_mx,e]-var_En2_prev[1:env$n_levs_mx],En2_prev[env$n_levs_mx:1,e]+var_En2_prev[env$n_levs_mx:1]),col="pink")
        lines(1:env$n_levs_mx,En2[1:env$n_levs_mx,e],type="l",col="red",ylim=ylim)
        points(env$n_levs_mx,En2[env$n_levs_mx+1,e],pch=21,bg="red",cex=2)
        lines(1:env$n_levs_mx,En2vo[1:env$n_levs_mx],col="blue")
        lines(1:env$n_levs_mx,En2vo[1:env$n_levs_mx],col="blue")
        lines(1:env$n_levs_mx,En2Vb[1:env$n_levs_mx,e],col="pink4")
        points(env$n_levs_mx,En2vo[env$n_levs_mx+1],pch=21,bg="blue",cex=2)
        points(env$n_levs_mx,En2Vb[env$n_levs_mx+1,e],pch=21,bg="pink4",cex=2)
        par(new=T)
        plot(1:env$n_levs_mx,rho[1:env$n_levs_mx,e],col="gold",axes=F,type="l",lty=1,lwd=2)
        axis(4)
        dev.off()
        print(paste("written file",fout))

        br<-c(0,1,2,4,8,16,32,64,128)
        col<-c("lightgray",rev(rainbow(7)))
        fouta<-file.path(dir_plot,paste0("rr1a.png"))
        png(file=fouta,width=1200,height=1200)
        image(rfxb,breaks=br,col=col,main="analysis")
        dev.off()
        foutb<-file.path(dir_plot,paste0("rr1b.png"))
        png(file=foutb,width=1200,height=1200)
        image(rfxb,breaks=br,col=col,main="analysis")
        for (i in 1:length(col)) {
          ix <- which(y_env$yo$value>=br[i] & y_env$yo$value<br[i+1])
          points(y_env$yo$x[ix],y_env$yo$y[ix],cex=0.8,pch=21,bg=col[i],col="darkgray")
        }
        dev.off()
        fout<-file.path(dir_plot,paste0("rr1_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
        system(paste0("convert +append ",fouta," ",foutb," ",fxb[e]," ",fout))
        system(paste0("rm -f ",fouta," ",foutb))
        print(paste("written file",fout))

#        aux <- getValues(rfxb,c_xy)
#        ylim=range(c(0,aux,y_env$yo$value))
#        fouta<-file.path(dir_plot,paste0("rr1x.png"))
#        png(file=fouta,width=1200,height=1200)
#        plot(y_env$yo$x,aux,pch=21,bg="lightgray",col="gray")
#        points(y_env$yo$x,y_env$yo$value,pch=21,bg="black")
#        dev.off()
#        foutb<-file.path(dir_plot,paste0("rr1y.png"))
#        png(file=foutb,width=1200,height=1200)
#        plot(y_env$yo$y,aux,pch=21,bg="lightgray",col="gray")
#        points(y_env$yo$y,y_env$yo$value,pch=21,bg="black")
#        dev.off()
#        foutc<-file.path(dir_plot,paste0("rr1scatt.png"))
#        png(file=foutc,width=1200,height=1200)
#        plot(y_env$yo$value,aux,pch=21,bg="lightgray",col="gray",xlim=ylim,ylim=ylim)
#        dev.off()
#        fout<-file.path(dir_plot,paste0("rr1xy_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
#        system(paste0("convert +append ",fouta," ",foutb," ",foutc," ",fout))
#        system(paste0("rm -f ",fouta," ",foutb," ",foutc))
#        print(paste("written file",fout))
      }
    }

    cat( paste0("time=",round(t1-t0,1), attr(t1-t0,"unit"),"\n"))
  } # end main loop

}
