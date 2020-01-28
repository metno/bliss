#+ Hybrid Ensemble Optimal Interpolation
henoi <- function( i,
                   mode = "analysis",
                   pmax = 200,
                   rloc_min = 0,
                   nens = 10,
                   rr_inf = 0.1,
                   statcov_backg = "gauss") { 
#------------------------------------------------------------------------------
# remember (wtf!) ‘t(x) %*% y’ (‘crossprod’) or ‘x %*% t(y)’ (‘tcrossprod’)
# mode = analysis/cvanalysis
# statcov_backg = gauss/exp
# global variables
#  obs_x. observation coordinates
#  grid_x. grid point coordinates
#  yo. observations
#  yb. background at observation locations
#  xb. background at grid points
#  HAf. forecast perturbations at observation locations
#  Af. forecast perturbations at grid points
#  Pfdiag. diagonal of the forecast perturbation cov mat
#  henoi_Dh_loc
#  henoi_Dh
#  henoi_eps2
#  var_o_coeff
#  pmax
#  rloc_min
#  nens
# mode == "cvanalysis"
# global variables
#  obs_x. observation coordinates
#  grid_x. observation coordinates 
#  yo. observations
#  yb. background at observation locations
#  xb. background at observation locations
#  HAf. forecast perturbations at observation locations
#  Af. forecast perturbations at observation locations
#  Pfdiag. diagonal of the forecast perturbation cov mat at obs loc
#------------------------------------------------------------------------------
  dist2 <- ( obs_x - grid_x[i])**2 + ( obs_y - grid_y[i])**2  # nobs_tot 
  rloc  <- exp( -0.5 * ( dist2 / henoi_Dh_loc[i]**2))         # nobs_tot  
  cond <- rloc > rloc_min
  if ( mode == "cvanalysis") cond[i] <- F 
  sel <- which( cond)
  p_i <- length( sel)
  # no observations in the surroundings
  if ( p_i == 0 ) {
    xa    <- xb[i]
    xidi  <- 0
    var_a <- Pfdiag[i]
    alpha <- 0
    zero  <- 0
  # found something ...
  } else {
    if ( p_i > pmax) sel <- order( rloc, decreasing = T)[1:pmax]
    p_i <- length( sel)
    # all observations say no-prec, then take a shortcut
    if ( !any( is.prec[sel])) {
      xa    <- NA
      xidi  <- NA
      var_a <- NA
      alpha <- NA
      zero  <- 1
    # now do the math ...
    } else {
      d.i     <- yo[sel] - yb[sel]                                 # p_i x 1
      HAf.i   <- HAf[sel,,drop=F]                                  # p_i x k
      obs_x.i <- obs_x[sel]                                        # p_i
      obs_y.i <- obs_y[sel]                                        # p_i
      dist2_obs <- outer( obs_x.i, obs_x.i, FUN="-")**2 + 
                   outer( obs_y.i, obs_y.i, FUN="-")**2            # p_i x p_i
      if ( statcov_backg == "gauss" ) {
        S.i <- exp( -0.5 *  dist2_obs / henoi_Dh[i]**2)            # p_i x p_i
        G.i <- exp( -0.5 * dist2[sel] / henoi_Dh[i]**2)            # p_i
      } else if ( statcov_backg == "exp" ) {
        S.i <- exp( -sqrt(  dist2_obs) / henoi_Dh[i])              # p_i x p_i
        G.i <- exp( -sqrt( dist2[sel]) / henoi_Dh[i])              # p_i
      }
      Sloc.i <- exp( -0.5 *  dist2_obs / henoi_Dh_loc[i]**2)       # p_i x p_i
      Gloc.i <- exp( -0.5 * dist2[sel] / henoi_Dh_loc[i]**2)       # p_i
      relw_Gloc.i <- Gloc.i / sum( Gloc.i)
      # Sakov and Bertino (2011) Eq.(8)
      Gf.i <- Gloc.i * 1/(nens-1) * tcrossprod(Af[i,,drop=F],HAf.i) # 1 x p_i
      # Sakov and Bertino (2011) Eq.(9)
      Sf.i <- Sloc.i * 1/(nens-1) * tcrossprod(HAf.i,HAf.i)        # p_i x p_i
      var_o_coeff.i <- var_o_coeff[sel]                            # p_i
      var_o_coeff_mean.i <- sum( var_o_coeff.i * relw_Gloc.i)
      alpha <- max( 0, 
               sum( d.i * d.i * relw_Gloc.i) / ( 1 + henoi_eps2[i]))
      Sb.i<- Sf.i + alpha * S.i                                    # p_i x p_i
      Gb.i<- Gf.i + alpha * G.i                                    # 1 x p_i 
      diag_R.i <- alpha * henoi_eps2[i] * var_o_coeff.i / sum(var_o_coeff.i*relw_Gloc.i)
      SbRinv.i <- try( chol2inv( chol( (Sb.i+diag(diag_R.i)))))    # p_i x p_i
      # slower alternative
      if ( !is.null( attr( SbRinv.i, "class"))) SbRinv.i <- solve( Sb.i + diag( diag_R.i))
      K.i   <- tcrossprod( Gb.i, SbRinv.i) # 1 x p_i 
      xa    <- xb[i] + tcrossprod( K.i, d.i)
      xidi  <- rowSums( K.i)
      var_a <- Pfdiag[i] + alpha - tcrossprod( K.i, Gb.i)
      zero  <- as.integer( xa < rr_inf)
    }
  } 
  if ( mode == "analysis" | mode == "cvanalysis" ) {
    return( c( xa, var_a, xidi, alpha, zero))
  } else {
    return( NULL)
  }
}

