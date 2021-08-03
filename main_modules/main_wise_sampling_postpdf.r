#+
main_wise_sampling_postpdf <- function( argv, y_env, fg_env, env, seed=NA,
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

  cat("Reconstruct analysis ensemble members on the state vector \n")
  Xa1     <- array( data=NA, dim=c(     env$m_dim, env$k_dim))
  for (e in 1:env$k_dim) {
    q<-rdyad
    q[]<-array(data=env$Xa[,e],dim=c(sqrt(length(env$Xa[,e])),sqrt(length(env$Xa[,e]))))

#    Xa1[,e] <- as.vector( unlist( dwt.2d( as.matrix(q), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)))
    dwt <- dwt.2d( as.matrix(q), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)

    # unlist the result and compute squared energies
    jj <- 0; ii <- 0
    for (l in env$n_levs_mn:env$n_levs_mx) {
      lh  <-  dwt[[1+3*(l-1)]]; hl  <-  dwt[[2+3*(l-1)]];  hh <-  dwt[[3+3*(l-1)]]
      ii <- jj + 1
      jj <- ii + length(lh) + length(hl) + length(hh) - 1
      Xa1[ii:jj,e] <- c(  lh,  hl, hh)
    }
    ll <- dwt[[4+3*(env$n_levs_mx-1)]] 
    Xa1[(jj+1):env$m_dim,e] <- ll
    rm( lh, hl, hh, ll, dwt, ii, jj)

    rm(q)
  }

  #
  cat("Sample realizations from the posterior distribution ")
  xa1        <- rowMeans(Xa1)
  Aa         <- Xa1 - xa1
  var_a1     <- 1/(env$k_dim-1) * rowSums( Aa*Aa) #* 1/dwt_res**2
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

  }
  cat("\n")

}
