#+
wise_analysis_loop <- function( argv, y_env, fg_env, env,
                                supob_nobs=50, supob_radius=1500, supob_q=0.99,
                                max_it=100, opttol=0.02,
                                En2_adj_fun="Gaussian",
                                En2_adj_min=0,
                                rescale_min_obs=100, rescale_min_cells=100,
                                plot=F, dir_plot=NA) {
#
#------------------------------------------------------------------------------

  t0a <- Sys.time()

  cat( "-- wise analysis --\n")

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

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))

  # ---~--------------
  # observations on the dyadic grid
  rfobs <- resample( env$mergeobs$r, rdyad, method="bilinear")
  rfobs[rfobs<y_env$rain] <- 0

  vobs <- getValues(rfobs)
  vobs_val <- !is.na(getValues(rfobs))
  vobs_nan <- is.na(getValues(rfobs))

  c_xy <- which( !is.na( getValues(rfobs)))
  cneg_xy <- which( is.na( getValues(rfobs)))
  env$rfobs <- rfobs

  rfobsidi <- rfobs
  rfobsidi[cneg_xy] <- 0
  rfobsidi[c_xy] <- 1

  rfnoidi <- rfobs
  rfnoidi[cneg_xy] <- 1
  rfnoidi[c_xy] <- 0

  rfobsidiwet <- rfobs
  rfobsidiwet[which(vobs  < y_env$rain & vobs_val)] <- 0
  rfobsidiwet[which(vobs >= y_env$rain & vobs_val)] <- 1
  rfobsidiwet[cneg_xy] <- 0

  rfobsididry <- rfobs
  rfobsididry[which(vobs  < y_env$rain & vobs_val)] <- 1
  rfobsididry[which(vobs >= y_env$rain & vobs_val)] <- 0
  rfobsididry[cneg_xy] <- 0

  rfobs[cneg_xy] <- 0
  
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
  sqrt_m_dim <- sqrt( env$m_dim)
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
  dwtobs  <- dwt.2d( as.matrix(rfobs),    wf=env$wf, J=n, boundary=env$boundary)
  dwtidi  <- dwt.2d( as.matrix(rfobsidi), wf=env$wf, J=n, boundary=env$boundary)
  dwtnoidi   <- dwt.2d( as.matrix(rfnoidi), wf=env$wf, J=n, boundary=env$boundary)
  dwtidiwet  <- dwt.2d( as.matrix(rfobsidiwet), wf=env$wf, J=n, boundary=env$boundary)
  dwtididry  <- dwt.2d( as.matrix(rfobsididry), wf=env$wf, J=n, boundary=env$boundary)
  # multires representations
  mrobs <- list()
  mrnoidi <- list()
  mridi <- list()
  mridiwet <- list()
  mrididry <- list()
  # unlist the result and compute squared energies
  vo    <- vector( mode="numeric", length=env$n_dim); vo[]<-NA
  nn <- 0
  for (l in 1:n) {
    ij <- ijFromLev( n, l, F)
    vo[ij[1,1]:ij[3,2]] <- c( dwtobs[[3*l-2]], dwtobs[[3*l-1]], dwtobs[[3*l]])
    mrobs[[l]] <- dwt.2d( as.matrix(rfobs), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    mrnoidi[[l]] <- dwt.2d( as.matrix(rfnoidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    mridi[[l]] <- dwt.2d( as.matrix(rfobsidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    mridiwet[[l]] <- dwt.2d( as.matrix(rfobsidiwet), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    mrididry[[l]] <- dwt.2d( as.matrix(rfobsididry), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    if ( any( mridi[[l]] > mrnoidi[[l]])) nn <- l
  }
  ij <- ijFromLev( n, n, T)
  vo[ij[1]:ij[2]] <- dwtobs[[3*n+1]]

  if (nn == 0) return(NULL)
  env$nn_dim <- max( ijFromLev( n, (nn+1), F))
  env$nn1_dim <- sqrt( diff(ijFromLev( n, nn, T)) + 1 )
  cat("\n")
#
#moveDown <- function(matin) {
#  dimout <- 2*dim(matin)[1]
#  matout <- array( data=NA, dim=c(dimout,dimout))
#  for (i in 1:dimout) {
#    matout[2*i,]   <- rep( matin[i,], each=2)
#    matout[2*i-1,] <- rep( matin[i,], each=2)
#  }
#}
  #--------------------------------------------------------    
  # Analysis
  cat("Analysis on the transformed space \n")

  Xb_dyad <- array( data=NA, dim=c( env$m_dim, env$k_dim))
  Xb_dyad_original <- array( data=NA, dim=c( env$m_dim, env$k_dim))

  #--------------------------------------------------------    
  # MAIN LOOP
  rfxb_acc <- vector( mode="numeric", length=length(getValues(rdyad))); rfxb_acc[] <- 0
  for (loop in 1:max_it) {

    t0 <- Sys.time()
    cat( paste(" loop", loop, " "))

    #--------------------------------------------------------    
    # set the background

    Ub     <- array( data=NA, dim=c( env$nn_dim, env$k_dim))
    Vb     <- array( data=NA, dim=c( env$nn_dim, env$k_dim))
    vo_Vb  <- array( data=NA, dim=c( env$nn_dim, env$k_dim))

    if (plot & loop == 1) fxb <- vector()
    
    for (e in 1:env$k_dim) {

      # first iteration - background from ensemble 
      if ( loop == 1) {
        i <- fg_env$ixs[e]
        # interpolate onto the dyadic grid
        # background at grid  points
        rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
        rfxb[rfxb<y_env$rain] <- 0
        Xb_dyad_original[,e]  <- getValues(rfxb)
        Xb_dyad[,e]  <- getValues(rfxb)
        rfxb_acc <- rfxb_acc + getValues(rfxb)
      # second iteration onwards - background from previous iteration 
      } else {
        # -- background at grid  points --
        rfxb[] <- array( data=env$Xa_dyad[,e], dim=c(sqrt_m_dim,sqrt_m_dim))
        rfxb[rfxb<y_env$rain] <- 0
        Xb_dyad[,e]  <- getValues(rfxb)
      }

      # background at observation points
      rfyb <- rdyad; rfyb[c_xy] <- rfxb[c_xy]
      # innovation 
      rfin <- rdyad; rfin[c_xy] <- rfobs[c_xy] - rfxb[c_xy]

      # plot the background
      if (plot & loop==1) {
        fxb[e] <- file.path(dir_plot,paste0("rfxb_",formatC(e,flag="0",width=2),".png"))
        br<-c(0,0.1,1,2,4,8,16,32,64,128)
        col<-c("lightgray",rev(rainbow(8)))
        png(file=fxb[e],width=1200,height=1200)
        image(rfxb,breaks=br,col=col,main="background")
        abline(h=seq(-1000000,10000000,by=100000),v=seq(-1000000,10000000,by=100000),lty=2)
        dev.off()
      }

      # transformation operator
      dwtxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtyb <- dwt.2d( as.matrix(rfyb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtin <- dwt.2d( as.matrix(rfin), wf=env$wf, J=(nn+1), boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in 1:(nn+1)) {
        ij <- ijFromLev( n, l, F)
        Ub[ij[1,1]:ij[3,2],e] <- c( dwtxb[[3*l-2]], dwtxb[[3*l-1]], dwtxb[[3*l]])
        Vb[ij[1,1]:ij[3,2],e] <- c( dwtyb[[3*l-2]], dwtyb[[3*l-1]], dwtyb[[3*l]])
        vo_Vb[ij[1,1]:ij[3,2],e] <- c( dwtin[[3*l-2]], dwtin[[3*l-1]], dwtin[[3*l]])
      }
    } # end loop over ensemble members

    # set the background: end
    #--------------------------------------------------------    

    #--------------------------------------------------------    
    # Background error covariance matrices
    Ab <- Ub - rowMeans(Ub)
    Db <- Vb - rowMeans(Vb)
    var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) 
    var_d1_yy <- 1/(env$k_dim-1) * rowSums( vo_Vb *vo_Vb) 
    # scaling coefficient used in spatial analysis
    coeff <- var_b1_xy / var_d1_yy
    coeff[!is.finite(coeff)] <- 0
    rm( Ab, Db)

#    #--------------------------------------------------------    
#    # sort of root mean squared error of the decomposed innovation
#    env$costf[loop] <- mean( sqrt( var_d1_yy))
#    cat( paste("rmse", round(env$costf[loop],5)), "\n")

    #--------------------------------------------------------    
    # Rescale the different components
    
    Ua  <- array( data=NA, dim=c(     env$nn_dim, env$k_dim))
    env$Xa_dyad  <- array( data=NA, dim=c( env$m_dim, env$k_dim))
    for (e in 1:env$k_dim) {
      cat(".")
      Ua[,e] <- Ub[,e] + coeff * vo_Vb[,e]
      # FROM HERE - FROM HERE

      rfxb[] <- Xb_dyad[,e]
      r <- rdyad
      dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      for (i in 1:(length(dwt)-1)) dwt[[i]][] <- 0
      ij <- ijFromLev( n, (nn+1), F)
      dwt[[3*(nn+1)-2]][] <- Ua[ij[1,1]:ij[1,2],e]
      dwt[[3*(nn+1)-1]][] <- Ua[ij[2,1]:ij[2,2],e]
      dwt[[3*(nn+1)]][]   <- Ua[ij[3,1]:ij[3,2],e]
      xbadj <- idwt.2d( dwt) 
      r[] <- array(data=as.matrix(xbadj),dim=c(sqrt_m_dim,sqrt_m_dim))
      dwtxbadj <- dwt.2d( as.matrix(r), wf=env$wf, J=nn, boundary=env$boundary)
      mrxbadj <- dwtxbadj[[3*nn+1]]
      mrxbadj[mrididry[[nn]]>mridiwet[[nn]]] <- 0
      mrxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=nn, boundary=env$boundary)[[3*nn+1]]
      mrxa <- mrnoidi[[nn]] * mrxb + mridi[[nn]] * mrxbadj


      for (l in nn:2) {
        r[] <- mrxa
        dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=l, boundary=env$boundary)
        for (i in 1:(length(dwt))) dwt[[i]][] <- 0
        dwt[[3*l+1]][] <- mrxa
        ij <- ijFromLev( n, l, F)
        dwt[[3*l-2]][] <- Ua[ij[1,1]:ij[1,2],e]
        dwt[[3*l-1]][] <- Ua[ij[2,1]:ij[2,2],e]
        dwt[[3*l]][]   <- Ua[ij[3,1]:ij[3,2],e]
        xbadj <- idwt.2d( dwt) 
        r[] <- array(data=as.matrix(xbadj),dim=c(sqrt_m_dim,sqrt_m_dim))
        dwtxbadj <- dwt.2d( as.matrix(r), wf=env$wf, J=(l-1), boundary=env$boundary)
        mrxbadj <- dwtxbadj[[3*(l-1)+1]]
        mrxbadj[mrididry[[l-1]]>mridiwet[[l-1]]] <- 0
        mrxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(l-1), boundary=env$boundary)[[3*(l-1)+1]]
        mrxa <- mrnoidi[[l-1]] * mrxb + mridi[[l-1]] * mrxbadj
      }
save(file="tmp.rdata",coeff,env,mrobs,mrnoidi,mridi,mridiwet,mrididry,mrxa,n,nn,rdyad,Xb_dyad,Ua,Ub,vo_Vb,y_env)
q()
      # FROM HERE

      mra <- mrxb[1:env$nn1_dim,1:env$nn1_dim,e] 
      for (l in nn:1) {
        for (i in 1:length(dwt_aux)) dwt_aux[[i]][] <- 0
        dwt_aux[[3*l+1]][] <- mra
        ij <- ijFromLev( n, l, F)
        dwt_aux[[3*l-2]][] <- Ua[ij[1,1]:ij[1,2],e]
        dwt_aux[[3*l-1]][] <- Ua[ij[2,1]:ij[2,2],e]
        dwt_aux[[3*l]][]   <- Ua[ij[3,1]:ij[3,2],e]
        dx_l <- idwt.2d( dwt_aux) 
        mrdx_l <- dwt.2d( as.matrix(dx_l), wf=env$wf, J=(l-1), boundary=env$boundary)[[3*l+1]]
      }
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
      env$Xa_dyad[,e] <- idwt.2d( dwt_aux)
      if (!is.na(y_env$rain)) env$Xa_dyad[,e][env$Xa_dyad[,e]<y_env$rain] <- 0
      # compute energies
      rfxb[] <- array(data=env$Xa_dyad[,e],dim=c(sqrt_m_dim,sqrt_m_dim))
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
    } # end of rescaling cycle

    # Constraint on energies
    for (e in 1:env$k_dim) {
      # Gaussian
      if ( En2_adj_fun == "Gaussian") {
        rho[,e] <- exp( -0.5 * (En2[,e] - En2_prev[,e])**2 / var_En2_prev)
      # SOAR
      } else if ( En2_adj_fun == "SOAR") {
        rho[,e] <- (1 + abs(En2[,e] - En2_prev[,e]) / sqrt(var_En2_prev)) * exp( -abs(En2[,e] - En2_prev[,e]) / sqrt(var_En2_prev) )
      # powerlaw
      } else if ( En2_adj_fun == "powerlaw") {
        rho[,e] <- 1 / ( 1 + 0.5*(En2[,e] - En2_prev[,e])**2 / var_En2_prev)
      # TOAR
      } else if ( En2_adj_fun == "TOAR") {
        rho[,e] <- (1 + abs(En2[,e] - En2_prev[,e]) / sqrt(var_En2_prev) + (En2[,e] - En2_prev[,e])**2 / (3*var_En2_prev)) * exp( -abs(En2[,e] - En2_prev[,e]) / sqrt(var_En2_prev) )
      }
      rho[,e][!is.finite(rho[,e])] <- 0
      rho[,e] <- En2_adj_min + (1-En2_adj_min) * rho[,e]

      dwt_aux <- dwt_out
      for (i in env$n_levs_mn:env$n_levs_mx) {
        ij <- ijFromLev( env$n_levs_mx, i, F)
        dwt_aux[[3*(i-1)+1]][] <- Ub[ij[1,1]:ij[1,2],e] + rho[i,e] * coeff[ij[1,1]:ij[1,2]] * ( vo[ij[1,1]:ij[1,2]] - Vb[ij[1,1]:ij[1,2],e])
        dwt_aux[[3*(i-1)+2]][] <- Ub[ij[2,1]:ij[2,2],e] + rho[i,e] * coeff[ij[2,1]:ij[2,2]] * ( vo[ij[2,1]:ij[2,2]] - Vb[ij[2,1]:ij[2,2],e])
        dwt_aux[[3*(i-1)+3]][] <- Ub[ij[3,1]:ij[3,2],e] + rho[i,e] * coeff[ij[3,1]:ij[3,2]] * ( vo[ij[3,1]:ij[3,2]] - Vb[ij[3,1]:ij[3,2],e])
        En2[i,e] <- mean( ( dwt_aux[[3*(i-1)+1]][] / 2**i)**2) + 
                    mean( ( dwt_aux[[3*(i-1)+2]][] / 2**i)**2) +
                    mean( ( dwt_aux[[3*(i-1)+3]][] / 2**i)**2)
      }
      ij <- ijFromLev( env$n_levs_mx, env$n_levs_mx, T)
      dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- Ub[ij[1]:ij[2],e] + rho[env$n_levs_mx+1,e] * coeff[ij[1]:ij[2]] * ( vo[ij[1]:ij[2]] - Vb[ij[1]:ij[2],e])
      En2[env$n_levs_mx+1,e] <- mean( ( dwt_aux[[3*(env$n_levs_mx-1)+4]][] / 2**i)**2)
      env$Xa_dyad[,e] <- idwt.2d( dwt_aux)
      # adjust for negative precip values
      if (!is.na(y_env$rain)) env$Xa_dyad[,e][env$Xa_dyad[,e]<y_env$rain] <- 0
      # spatial analysis is useful only for precip events that includes observations 
      rfxb[] <- array( data=env$Xa_dyad[,e], dim=c(sqrt_m_dim,sqrt_m_dim))
      aclump <- clump(rfxb)
      oclump <- aclump[c_xy]
      fro <- as.data.frame( table( oclump), stringsAsFactors=F) 
      if ( length( ixo <- which( !is.na(fro[,1]) & fro[,2] >= 100)) > 0) {
        oclump_ok <- as.integer(fro[ixo,1])
      } else {
        oclump_ok <- integer(0)
      }
      fra <- freq(aclump)
      # remove clumps of YESprec cells less than 100 cells of the dyadic domain or not including enough obs
      ixa <- which( !is.na(fra[,1]) & !is.na(fra[,2]) & 
                   ( (fra[,2]<=100) | !(fra[,1] %in% oclump_ok)) )
      env$Xa_dyad[which(getValues(t(aclump)) %in% fra[ixa,1]),e] <- 0
    } # end of Constraint on energies

    if (plot) {
      for (e in 1:env$k_dim) {
        rfxb[] <- array(data=env$Xa_dyad[,e],dim=c(sqrt_m_dim,sqrt_m_dim))
        if (loop==1) ylim_en2<-range(c(En2,En2_prev+var_En2_prev))

        # En2
        fout<-file.path(dir_plot,paste0("en2_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
        sd_En2_prev <- sqrt(var_En2_prev)
        png(file=fout,width=800,height=800)
        par(mar=c(5,5,1,3))
#        plot(1:env$n_levs_mx,En2Ub[1:env$n_levs_mx,e],type="l",col="red",ylim=ylim)
        plot(1:env$n_levs_mx,En2[1:env$n_levs_mx,e],type="l",col="red",ylim=ylim_en2,xlab="level",ylab="Energy^2")
        abline(h=0,col="black",lwd=2)
        points(env$n_levs_mx,En2[env$n_levs_mx+1,e],pch=21,bg="red",cex=2)
        polygon( c(1:env$n_levs_mx,env$n_levs_mx:1),
                 c(En2_prev[1:env$n_levs_mx,e]-2*sd_En2_prev[1:env$n_levs_mx],
                   En2_prev[env$n_levs_mx:1,e]+2*sd_En2_prev[env$n_levs_mx:1]),
                 density=10,border="pink", col="pink",angle=45,lwd=3)
        polygon( c(1:env$n_levs_mx,env$n_levs_mx:1),
                 c(En2_prev[1:env$n_levs_mx,e]-sd_En2_prev[1:env$n_levs_mx],
                   En2_prev[env$n_levs_mx:1,e]+sd_En2_prev[env$n_levs_mx:1]),
                 border="pink", col="pink", density=10, angle=-45,lwd=3)
        lines(1:env$n_levs_mx,En2_prev[1:env$n_levs_mx,e],col="pink4",lwd=6)
        points(env$n_levs_mx,En2_prev[env$n_levs_mx+1,e],pch=21,bg="pink4",cex=2)
        lines(1:env$n_levs_mx,En2[1:env$n_levs_mx,e],type="l",col="red",lwd=6)
        points(env$n_levs_mx,En2[env$n_levs_mx+1,e],pch=21,bg="red",cex=2)
        lines(1:env$n_levs_mx,En2vo[1:env$n_levs_mx],col="blue",lwd=2)
        lines(1:env$n_levs_mx,En2Vb[1:env$n_levs_mx,e],col="pink4",lwd=2)
        points(env$n_levs_mx,En2vo[env$n_levs_mx+1],pch=21,bg="blue",cex=2)
        points(env$n_levs_mx,En2Vb[env$n_levs_mx+1,e],pch=21,bg="pink4",cex=2)
        par(new=T)
        plot(1:env$n_levs_mx,rho[1:env$n_levs_mx,e],ylim=c(0,1),col="gold",axes=F,type="l",lty=1,lwd=6,xlab="",ylab="")
        abline(h=0,col="black",lwd=2,lty=2)
        axis(4)
        dev.off()
        print(paste("written file",fout))

        # rr1 map
        br<-c(0,0.1,1,2,4,8,16,32,64,128)
        col<-c("lightgray",rev(rainbow(8)))
        fouta<-file.path(dir_plot,paste0("rr1a.png"))
        png(file=fouta,width=1200,height=1200)
        image(rfxb,breaks=br,col=col,main="analysis")
        abline(h=seq(-1000000,10000000,by=100000),v=seq(-1000000,10000000,by=100000),lty=2)
        dev.off()
        xy<-xyFromCell( rfobs, c_xy)
        yso_x <- xy[,1]
        yso_y <- xy[,2]
        yso_val <- getValues( rfobs)[c_xy]
        foutb<-file.path(dir_plot,paste0("rr1b.png"))
        png(file=foutb,width=1200,height=1200)
        image(rfxb,breaks=c(0,1),col="white",main="observations")
        for (i in 1:length(col)) {
          ix <- which(yso_val>=br[i] & yso_val<br[i+1])
          points(yso_x[ix],yso_y[ix],cex=0.8,pch=21,bg=col[i],col=col[i])
        }
        abline(h=seq(-1000000,10000000,by=100000),v=seq(-1000000,10000000,by=100000),lty=2)
        dev.off()
        fout<-file.path(dir_plot,paste0("rr1_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
        system(paste0("convert +append ",fouta," ",foutb," ",fxb[e]," ",fout))
        system(paste0("rm -f ",fouta," ",foutb))
        print(paste("written file",fout))

        # rr1 graph 
        aux <- extract( rfxb, cbind(yso_x,yso_y))
        ylim_rr1 <- range(c(0,yso_val))
        fouta<-file.path(dir_plot,paste0("rr1x.png"))
        png(file=fouta,width=1200,height=1200)
        plot(yso_x,aux,pch=21,bg="lightgray",col="gray")
        points(yso_x,yso_val,pch=21,bg="black")
        dev.off()
        foutb<-file.path(dir_plot,paste0("rr1y.png"))
        png(file=foutb,width=1200,height=1200)
        plot(yso_y,aux,pch=21,bg="lightgray",col="gray")
        points(yso_y,yso_val,pch=21,bg="black")
        dev.off()
        foutc<-file.path(dir_plot,paste0("rr1scatt.png"))
        png(file=foutc,width=1200,height=1200)
        plot(yso_val,aux,pch=21,bg="lightgray",col="gray",xlim=ylim_rr1,ylim=ylim_rr1)
        dev.off()
        fout<-file.path(dir_plot,paste0("rr1xy_",formatC(loop,width=2,flag="0"),"_",formatC(e,width=2,flag="0"),".png"))
        system(paste0("convert +append ",fouta," ",foutb," ",foutc," ",fout))
        system(paste0("rm -f ",fouta," ",foutb," ",foutc))
        print(paste("written file",fout))
      }
    }

    #
    #--------------------------------------------------------    
    # DeltaEn2idx - index related to the mean variation of the energy over the spatial scales (0=no variation, 1=mean variation of the order of the variability)
    # strict condition for convergence
    env$costf[loop] <- 1 - min( rho, na.rm=T)
    # less strict condition for convergence
#    env$costf[loop] <- 1 - mean( rho, na.rm=T)
    t1 <- Sys.time()
    cat( paste( "time=", round(t1-t0,1), attr(t1-t0,"unit"), "costf - Delta En2 index =", round(env$costf[loop],5), "\n"))
    # break out of the main loop early if variations  
    if ( env$costf[loop] < opttol) break

  } # end main loop

#save(file="tmp.rdata",env,Xb_dyad,sqrt_m_dim,rfxb,c_xy)
  cat( paste( "Add the background data where no observations are available\n"))
  for (e in 1:env$k_dim) {
    cat(".")
    rfxb[] <- Xb_dyad[,e]
    bclump <- clump(rfxb)
    oclump <- bclump[c_xy]
    fro <- as.data.frame( table( oclump), stringsAsFactors=F) 
    if ( length( ixo <- which( !is.na(fro[,1]) & fro[,2] > 100)) > 0) {
      oclump_ok <- as.integer(fro[ixo,1])
    } else {
      oclump_ok <- integer(0)
    }
    Xb_dyad[which(getValues(bclump) %in% oclump_ok),e] <- 0
    rfxb[] <- Xb_dyad[,e]
    env$Xa_dyad[,e] <- env$Xa_dyad[,e] + getValues(t(rfxb))
  }
  cat( paste( "\n", "End of analysis loop time=", round(t1-t0,1), attr(t1-t0,"unit"), "\n"))
}
