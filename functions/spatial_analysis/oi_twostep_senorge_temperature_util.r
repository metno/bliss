#+
vertical_profile_at_centroid_senorge2018 <- function( i) {

  if ( length( ix <- which( envtmp$nn2$nn.idx[i,] != 0)) == 0) 
    return(c( NA, NA, NA, NA, NA, NA))
 
  z_opt <- y_env$super_yo$z[envtmp$nn2$nn.idx[i,ix]]
  v_opt <- y_env$super_yo$value[envtmp$nn2$nn.idx[i,ix]]

  dz <- as.numeric( diff( quantile( z_opt, probs=c(0.05,0.95))))

  if (dz < argv$oi2step.bg_vertprof_dzmin) {
    par  <- mean( v_opt, na.rm=T)
    opt  <- optimize( f=vert_prof_basic_opt, 
                      interval=c(argv$oi2step.bg_vertprof_vmin, argv$oi2step.bg_vertprof_vmax),
                      z=z_opt, v=v_opt, gamma=argv$oi2step.bg_vertprof_gamma)
    tpar <- c( 1, opt$minimum, argv$oi2step.bg_vertprof_gamma, NA, NA, NA)
  } else {
    ui <- matrix( ncol=5, nrow=10, data=0)
    ci <- vector( length=10)
    q  <- as.numeric( quantile( z_opt, probs=c( 0.1, 0.25, 0.75)))
    v_range <- range(v_opt,na.rm=T)
    for (i in 1:5) {
      ui[(2*i-1),i] <- 1
      ui[(2*i),i]   <- (-1)
    }
    ci[1:2] <- c( v_range[1], -v_range[2])
    aux <- sort( c( 0.5*argv$oi2step.bg_vertprof_gamma, 2*argv$oi2step.bg_vertprof_gamma))
    ci[3:4] <- c( aux[1], -aux[2])
    ci[5:6] <- c(      0, -20)
    if ( q[1] == q[2]) {
      aux <- q[2] + mean(z_opt,na.rm=T) / 10
    } else {
      aux <- q[2]
    }
    ci[7:8] <- c( -10, -aux)
    aux <- q[3] - q[2] 
    ci[9:10] <- c( 0, -aux)
    par <- c( mean( v_range),
              argv$oi2step.bg_vertprof_gamma,
              5,
              q[1],
             -ci[10] / 2)
    opt <- constrOptim( theta=par,
                        f=vert_prof_Frei_opt,
                        ui=ui,
                        ci=ci,
                        grad=NULL,
                        z=z_opt, v=v_opt)
    tpar <- c( 2, opt$par[1], opt$par[2], opt$par[3], opt$par[4], opt$par[5])
  } 
 
  return(tpar)

}

#+
blend_vertical_profiles_senorge2018 <- function( i, 
                                                 corr = "soar",
                                                 dh = 10000) {
  if (i==1) cat(".")

  rloc <- corr1d( envtmp$nn2$nn.dist[i,], dh, corr)
  values <- vector( mode="numeric", length=y_env$centroids$n)
  if ( envtmp$n1 > 0) 
    values[envtmp$ix1] <- vert_prof_basic( rep( envtmp$z[i], envtmp$n1), 
                           y_env$centroids$vert_prof[envtmp$ix1,2], 
                           y_env$centroids$vert_prof[envtmp$ix1,3]) 
  if ( envtmp$n2 > 0) 
    values[envtmp$ix2] <- vert_prof_Frei_1( rep( envtmp$z[i], envtmp$n2), 
                           y_env$centroids$vert_prof[envtmp$ix2,2], 
                           y_env$centroids$vert_prof[envtmp$ix2,3],
                           y_env$centroids$vert_prof[envtmp$ix2,4],
                           y_env$centroids$vert_prof[envtmp$ix2,5],
                           y_env$centroids$vert_prof[envtmp$ix2,6])

  sum(rloc * values) / sum(rloc)

}

#+ vertical profile of temperature (linear)
vert_prof_basic <- function( z, t0, gamma) {
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  return( t0 + gamma * z)
}


#+ vertical profile of temperature (Frei, 2014)
vert_prof_Frei<-function( z, t0, gamma, a, h0, h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region 
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
#  a= numeric. Temperature contrast or inversion strength [degC] 
#  h0= numeric. z where inversion starts [m]
#  h1i= numeric. h0+h1i is z where inversion stops [m]
#       (Frei uses h1 directly, I use an increment to h0 so to avoid ending
#        up with h1<=h0 during the optimization)
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  t   <- z
  t[] <- NA
  h1  <- h0 + h1i
  if ( length( (z.le.h0 <- which( z <= h0))) > 0) 
    t[z.le.h0] <- t0 + gamma * z[z.le.h0] - a 
  if ( length( (z.ge.h1 <- which( z >= h1))) > 0) 
    t[z.ge.h1] <- t0 + gamma * z[z.ge.h1]
  if ( length( (z.in <- which( z > h0 & z < h1)))>0) 
    t[z.in] <- t0 + gamma * z[z.in] - a/2 * (1 + cos(pi*(z[z.in]-h0)/h1i))
  return(t)
}

#+ vertical profile of temperature (Frei, 2014)
vert_prof_Frei_1<-function( z, t0, gamma, a, h0, h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region 
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
#  a= numeric. Temperature contrast or inversion strength [degC] 
#  h0= numeric. z where inversion starts [m]
#  h1i= numeric. h0+h1i is z where inversion stops [m]
#       (Frei uses h1 directly, I use an increment to h0 so to avoid ending
#        up with h1<=h0 during the optimization)
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  t   <- z
  t[] <- NA
  h1  <- h0 + h1i
  if ( length( (z.le.h0 <- which( z <= h0))) > 0) 
    t[z.le.h0] <- t0[z.le.h0] + gamma[z.le.h0] * z[z.le.h0] - a[z.le.h0] 
  if ( length( (z.ge.h1 <- which( z >= h1))) > 0) 
    t[z.ge.h1] <- t0[z.ge.h1] + gamma[z.ge.h1] * z[z.ge.h1]
  if ( length( (z.in <- which( z > h0 & z < h1)))>0) 
    t[z.in] <- t0[z.in] + gamma[z.in] * z[z.in] - a[z.in]/2 * (1 + cos(pi*(z[z.in]-h0[z.in])/h1i[z.in]))
  return(t)
}

#+ cost function used for optimization of tvertprof parameter
vert_prof_Frei_opt<-function(par, z, v) {
  return( log( sqrt( mean( (vert_prof_Frei(z=z,t0=par[1],gamma=par[2],a=par[3],h0=par[4],h1i=par[5]) - v)**2))))
}


#+ cost function used for optimization of tvertprof parameter
vert_prof_basic_opt <- function( par, z, v, gamma) {
#  te <- vert_prof_basic( z=zopt, t0=par[1], gamma=gamma)
  return( log( sqrt( mean( (vert_prof_basic( z=z, t0=par, gamma=gamma) - v)**2))))
}

#+
Theil_Sen_regression <-function( x, y) {
# y = res[1] + res[2] * x
#------------------------------------------------------------------------------
# REF: Wilks (2019) p. 283
  bnum <- outer( y, y, FUN="-")
  bden <- outer( x, x, FUN="-")
  bset <- bnum / bden
  b <- median( bset[row(bset)>col(bset) & is.finite(bset)], na.rm=T)
  if ( !(!is.finite(b) | is.na(b) | is.null(b))) {
    residuals <- y - b * x
    a   <- median( residuals)
    res <- c( a, b)
  } else {
    res <- c( NA, NA)
  }
  res
}
