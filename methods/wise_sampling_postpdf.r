#+
wise_sampling_postpdf <- function( argv, y_env, fg_env, u_env, env, 
                                   resample=F,
                                   seed=NA,
                                   plot=F,dir_plot=NA) {
#
#------------------------------------------------------------------------------
  
  options( warn = 2)
  
  if (resample) {

    library(waveslim)

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

    cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))

    #
    # compute score for each analysis field
    cat (" compute score wrt u_env for each analysis field \n")
    # reference observed field
    # interpolate onto the dyadic grid
    uo <- getValues( resample( u_env$uo[[1]]$r_main, rdyad, method="bilinear"))
    uo[uo<u_env$rain] <- 0
    uo[uo>=u_env$rain] <- 1
    score <- array( data=NA, dim=c(1,env$k_dim))
    for (e in 1:env$k_dim) {
      r <- env$Xa_dyad[,e]
      r[r<u_env$rain] <- 0
      r[r>=u_env$rain] <- 1
      a <- as.numeric( length( which(r==1 & uo==1)))
      b <- as.numeric( length( which(r==1 & uo==0)))
      c <- as.numeric( length( which(r==0 & uo==1)))
      d <- as.numeric( length( which(r==0 & uo==0)))
      if (argv$wise_align_mode == "ets") {
        a_random <- (a+c)*(a+b) / (a+b+c+d)
        score[1,e] <- (a-a_random) / (a+c+b-a_random)
      } else if (argv$wise_align_mode == "maxoverlap") {
        score[1,e] <- a/length(r)
      }
      cat( ".")
    }
    rm(uo,r)
    weights <- score / sum(score)

    # 
    # Initialization of the wavelet structures
    dwt_out <- dwt.2d( as.matrix(rdyad), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
    for (i in 1:length(dwt_out)) dwt_out[[i]][] <- 0

    #
    cat("Reconstruct analysis ensemble members on the state vector \n")
    Ua  <- array( data=NA, dim=c( env$n_dim, env$k_dim))
    for (e in 1:env$k_dim) {
      dwt <- dwt.2d( as.matrix( array( data=env$Xa_dyad[,e], dim=c( sqrt_m_dim, sqrt_m_dim))), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in env$n_levs_mn:env$n_levs_mx) {
        lh  <-  dwt[[1+3*(l-1)]]; hl  <-  dwt[[2+3*(l-1)]];  hh <-  dwt[[3+3*(l-1)]]
        ii <- jj + 1
        jj <- ii + length(lh) + length(hl) + length(hh) - 1
        Ua[ii:jj,e] <- c(  lh,  hl, hh)
      }
      ll <- dwt[[4+3*(env$n_levs_mx-1)]] 
      Ua[(jj+1):env$m_dim,e] <- ll
    }
    rm( lh, hl, hh, ll, dwt, ii, jj)

    #
    cat("Sample realizations from the posterior distribution ")
    Ua_mean_dyad <- as.vector( t( tcrossprod( weights, Ua)))
    Ua_var_dyad <- tcrossprod( weights, (Ua - Ua_mean_dyad)**2)
    env$Xa_dyad <- array( data=NA, dim=c( env$m_dim, env$a_dim))
    if (!is.na(seed)) set.seed(seed)
    for (a in 1:env$a_dim) {
      cat(".")
#      aux     <- rnorm( env$m_dim, mean=Ua_mean_dyad, sd=sqrt(Ua_var_dyad) )
      aux     <- rnorm( env$n_dim, mean=Ua_mean_dyad, sd=sqrt(Ua_var_dyad) )
      dwt_aux <- dwt_out
      for (i in env$n_levs_mn:env$n_levs_mx) {
        ij <- ijFromLev( env$n_levs_mx, i, F)
        dwt_aux[[3*(i-1)+1]][] <- aux[ij[1,1]:ij[1,2]]
        dwt_aux[[3*(i-1)+2]][] <- aux[ij[2,1]:ij[2,2]]
        dwt_aux[[3*(i-1)+3]][] <- aux[ij[3,1]:ij[3,2]]
      }
      ij <- ijFromLev( env$n_levs_mx, env$n_levs_mx, T)
      dwt_aux[[3*(env$n_levs_mx-1)+4]][] <- aux[ij[1]:ij[2]]
      # reconstruct analysis
      env$Xa_dyad[,a] <- idwt.2d( dwt_aux)
      if (!is.na(y_env$rain)) env$Xa_dyad[,a][env$Xa_dyad[,a]<y_env$rain] <- 0
    }
  }

  cat("\n")

  if (plot) {
    for (a in 1:env$a_dim) {
      rdyad[] <- array(data=env$Xa_dyad[,a],dim=c(sqrt_m_dim,sqrt_m_dim))
      br<-c(0,1,2,4,8,16,32,64,128)
      col<-c("lightgray",rev(rainbow(7)))
      fouta<-file.path(dir_plot,paste0("rr1apost_",formatC(a,width=2,flag="0"),".png"))
      png(file=fouta,width=1200,height=1200)
      image(rdyad,breaks=br,col=col,main="analysis")
      abline(h=seq(-1000000,10000000,by=100000),v=seq(-1000000,10000000,by=100000),lty=2)
      dev.off()
      print(paste("written file",fouta))
    }
  }


}
