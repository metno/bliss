#+ Hybrid Ensemble Optimal Interpolation
henoi_alt<-function(i,
                    mode="analysis",
                    pmax=200,
                    rloc_min=0,
                    var_f_thresh=0.001,
                    nens=10,
                    nrand=1000) {
#------------------------------------------------------------------------------
# remember (wtf!) ‘t(x) %*% y’ (‘crossprod’) or ‘x %*% t(y)’ (‘tcrossprod’)
# mode == "analysis"
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
#  is.prec
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
  dist2<-(obs_x-grid_x[i])**2 + (obs_y-grid_y[i])**2  # nobs_tot 
  rloc<-exp(-0.5*(dist2/henoi_Dh_loc[i]**2)) # nobs_tot  
  cond<-rloc>rloc_min
  if (mode=="cvanalysis") cond[i]<-F 
  sel<-which(cond)
  p_i<-length(sel)
  # no observations in the surroundings
  if (p_i==0) {
    xa<-xb[i]
    xidi<-0
    var_a<-Pfdiag[i]
#    alpha<-0
    zero<-F
  # found something ...
  } else {
    if (p_i>pmax) sel<-order(rloc,decreasing=T)[1:pmax]
    p_i<-length(sel)
    # all observations say no-prec, then take a shortcut
    if (!any(is.prec[sel])) {
      xa<-NA
      xidi<-NA
      var_a<-NA
#      alpha<-NA
      zero<-T
    # now do the math ...
    } else {
      obs_x.i<-obs_x[sel]                                           # p_i
      obs_y.i<-obs_y[sel]                                           # p_i
      dist2_obs<-outer(obs_x.i,obs_x.i,FUN="-")**2 + 
                 outer(obs_y.i,obs_y.i,FUN="-")**2                  # p_i x p_i
      Sb.i<-exp(-0.5*dist2_obs/henoi_Dh[i]**2)                      # p_i x p_i
      Gb.i<-exp(-0.5*dist2[sel]/henoi_Dh[i]**2)                     # p_i
      var_o_coeff.i<-var_o_coeff[sel]                               # p_i
      diag_R.i<- henoi_eps2[i] * var_o_coeff.i                      # p_i
      SbRinv.i<-try(chol2inv(chol( (Sb.i+diag(diag_R.i)) )))        # p_i x p_i
      # slower alternative
      if (!is.null(attr(SbRinv.i,"class"))) SbRinv.i<-solve( Sb.i+diag(diag_R.i) )
      K.i<-tcrossprod(Gb.i,SbRinv.i)                                # 1 x p_i 
      ens<-rbind(HAf[sel,,drop=F],Af[i,,drop=F])                    # p_i+1 x k
      Pens<-1/(nens-1)*tcrossprod(ens,ens)+diag(0.0001,dim(ens)[1]) # p_i+1 x p_i+1 
      eig <- eigen(Pens,symmetric=T) 
      Pens_sqrt <- tcrossprod( tcrossprod( eig$vectors,             # p_i+1 x p_i+1 
                                           diag(sqrt(eig$values))), 
                               eig$vectors )
      z<-array(rnorm((p_i+1)*nrand),dim=c((p_i+1),nrand))           # p_i+1 x nrand
      yxb.i<-tcrossprod(Pens_sqrt,t(z))+c(yb[sel],xb[i])            # p_i+1 x nrand 
      yxb.i[yxb.i<zerot]<-zerot
      d.i<-yo[sel]-yxb.i[1:p_i,]                                    # p_i+1 x nrand
      xa_real<-yxb.i[(p_i+1),]+t(tcrossprod(K.i,t(d.i)))            # nrand
      xa<-mean(xa_real)
      xidi<- rowSums(K.i)
      var_a<-var(xa_real) 
      zero<-F
    }
  } 
  if ( mode=="analysis" | mode=="cvanalysis" ) {
    return(c(xa,var_a,xidi,NA,zero))
  } else {
    return(NULL)
  }
}

