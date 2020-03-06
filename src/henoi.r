#+ Hybrid Ensemble Optimal Interpolation
henoi <- function( i,
                   mode = "analysis",
                   pmax = 200,
                   rloc_min = 0,
                   nens = 10,
                   statcov_backg = "gauss") { 
#------------------------------------------------------------------------------
# remember (wtf!) ‘t(x) %*% y’ (‘crossprod’) or ‘x %*% t(y)’ (‘tcrossprod’)
# mode = analysis/cvanalysis
# statcov_backg = gauss/exp
# global variables
#  obs_x. observation coordinates
#  grid_x. grid point coordinates
#  yo. observations
#  yb_best. background at observation locations
#  xb_best. background at grid points
#  HAf. forecast perturbations at observation locations
#  Af. forecast perturbations at grid points
#  Pfdiag. diagonal of the forecast perturbation cov mat
#  henoi_Dh_loc
#  henoi_Dh
#  henoi_eps2
#  henoi_alpha
#  var_o_coeff
#  pmax
#  rloc_min
#  nens
# mode == "cvanalysis"
# global variables
#  obs_x. observation coordinates
#  grid_x. observation coordinates 
#  yo. observations
#  yb_best. background at observation locations
#  xb_best. background at observation locations
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
    xa    <- xb_best[i]
    xidi  <- 0
    var_a <- Pfdiag[i]
    var_u <- 0
  # found something ...
  } else {
    if ( p_i > pmax) sel <- order( rloc, decreasing = T)[1:pmax]
    p_i <- length( sel)
    HAf.i   <- HAf[sel,,drop=F]                                  # p_i x k
    obs_x.i <- obs_x[sel]                                        # p_i
    obs_y.i <- obs_y[sel]                                        # p_i
    dist2_obs <- outer( obs_x.i, obs_x.i, FUN="-")**2 + 
                 outer( obs_y.i, obs_y.i, FUN="-")**2            # p_i x p_i
    Sloc.i <- exp( -0.5 *  dist2_obs / henoi_Dh_loc[i]**2)       # p_i x p_i
    Gloc.i <- exp( -0.5 * dist2[sel] / henoi_Dh_loc[i]**2)       # p_i
    relw_Gloc.i <- Gloc.i / sum( Gloc.i)
    # Sakov and Bertino (2011) Eq.(8)
    Gf.i <- Gloc.i * 1/(nens-1) * tcrossprod(Af[i,,drop=F],HAf.i) # 1 x p_i
    # Sakov and Bertino (2011) Eq.(9)
    Sf.i <- Sloc.i * 1/(nens-1) * tcrossprod(HAf.i,HAf.i)        # p_i x p_i
    var_f <- henoi_alpha[i] * sum( diag(Sf.i) * relw_Gloc.i)
    var_o_coeff.i <- var_o_coeff[sel]                            # p_i
    d.i     <- yo[sel] - yb_best[sel]                            # p_i x 1
    var_ob_emp <- henoi_alpha[i] * sum( d.i * d.i * relw_Gloc.i)
#    if ( (var_ob_emp <=  (henoi_eps2[i] * var_f)) & (var_f == 0) ) { 
    if ( ( (var_ob_emp/(1+henoi_eps2[i])) <=  var_f) & (var_f == 0) ) { 
      xa    <- xb_best[i]
      xidi  <- NA
      var_a <- NA
      var_u <- NA
    } else {
#      if ( var_ob_emp <=  (henoi_eps2[i] * var_f) ) { 
      if ( (var_ob_emp/(1+henoi_eps2[i])) <=  var_f ) { 
        var_u <- 0
        Sb.i<- Sf.i                                     # p_i x p_i
        Gb.i<- Gf.i                                     # 1 x p_i
      } else {
        if ( statcov_backg == "gauss" ) {
          S.i <- exp( -0.5 *  dist2_obs / henoi_Dh[i]**2)            # p_i x p_i
          G.i <- exp( -0.5 * dist2[sel] / henoi_Dh[i]**2)            # p_i
        } else if ( statcov_backg == "exp" ) {
          S.i <- exp( -sqrt(  dist2_obs) / henoi_Dh[i])              # p_i x p_i
          G.i <- exp( -sqrt( dist2[sel]) / henoi_Dh[i])              # p_i
        }
#        var_u <- ( var_ob_emp - henoi_eps2[i] * var_f ) / ( 1 + henoi_eps2[i])
        var_u <- var_ob_emp / ( 1 + henoi_eps2[i] ) - var_f
        Sb.i<- Sf.i + var_u * S.i                                    # p_i x p_i
        Gb.i<- Gf.i + var_u * G.i                                    # 1 x p_i
      }
      var_b <- var_f + var_u 
      diag_R.i <- var_b * henoi_eps2[i] * var_o_coeff.i
      SbRinv.i <- try( chol2inv( chol( (Sb.i + diag( diag_R.i, nrow=p_i, ncol=p_i)))))    # p_i x p_i
      # slower alternative
      if ( !is.null( attr( SbRinv.i, "class"))) 
        SbRinv.i <- solve( Sb.i + diag( diag_R.i, nrow=p_i, ncol=p_i))
      K.i   <- tcrossprod( Gb.i, SbRinv.i) # 1 x p_i 
      xa    <- xb_best[i]
      if ( any( d.i != 0)) xa <- xa + tcrossprod( K.i, d.i)
      xidi  <- rowSums( K.i)
      var_a <- Pfdiag[i] + var_u - tcrossprod( K.i, Gb.i)
    }
  } 
  if ( mode == "analysis" | mode == "cvanalysis" ) {
    return( c( xa, var_a, xidi, var_u))
  } else {
    return( NULL)
  }
}

