#+ Gaussian anamorphosis, from Gamma to Gaussian
gamma_anamorphosis<-function(x,
                             shape,
                             rate=0.1,
                             small_const=NA,
                             gauss_mean=0,
                             gauss_sd=1) {
# x. vector with the Gamma values to be transformed into Gaussians
# shape. vector. Gamma shape
# rate. vector. Gamma rate
# small_const. small constant to be added to the Gamma values to ensure a 
#              regular transformation (useful for 0s)
# gauss_mean. Gaussian's mean value
# gauss_sd. Gaussian's standard deviation
#------------------------------------------------------------------------------ 
  if (!is.na(small_const)) x<-x+small_const
  qnorm(pgamma(x,shape=shape,rate=rate),mean=gauss_mean,sd=gauss_sd)
}

#+ Gaussian anamorphosis, back-transformation. from Gaussian to Gamma. Expected value
inv_gamma_anamorphosis<-function(values,
                                 sigma=NA,
                                 shape,
                                 rate,
                                 small_const=NA,
                                 threshold0=0,
                                 gauss_mean=0,
                                 gauss_sd=1) {
# values. vector. Gaussian values to be transformed into Gamma
# sigma. vector. Gaussian standard deviations (NAs if direct back-transform)
# shape. vector. Gamma shape
# rate. vector. Gamma rate
# small_const. small constant to be removed to the Gamma values
# gauss_mean. Gaussian's mean value
# gauss_sd. Gaussian's standard deviation
#------------------------------------------------------------------------------ 
  # direct back-transformation
  if (any(is.na(sigma))) {
    res<-qgamma( pnorm(values,mean=gauss_mean,sd=gauss_sd), shape=shape, rate=rate)
  # PDF-based back-transformation
  } else {
    res<-apply( cbind( values, sqrt(sigma), shape, rate ),
                MAR=1,
                FUN=function(x){
                  y<-x[1]; sd<-x[2]; s<-x[3]; r<-x[4]
                  if (sd==0) return(y)
##                  seq <- seq( (y-3*sd), (y+3*sd), length=100)
##                  dy  <- abs(seq[2]-seq[1])
##                  aux <- inv_gamma_anamorphosis(x=seq, shape=s, rate=r)
##                  ix  <- which(is.finite(aux))
##                  sum( aux[ix] * dnorm(seq[ix],mean=y,sd=sd) * dy) 
                  seq <- qnorm(seq(0.001,0.999,length=400),mean=y,sd=sd)
                  aux <- inv_gamma_anamorphosis(values=seq, shape=s, rate=r)
                  ix  <- which(is.finite(aux))
                  n<-length(ix)
                  xgamma <- aux[ix]
                  wnorm  <- dnorm(seq[ix],mean=y,sd=sd)
                  xnorm  <- seq[ix]
                  sum( xgamma[2:n] * wnorm[2:n] * diff(xnorm)) 
#                  mean(inv_gamma_anamorphosis(x= qnorm(seq(0.001,0.999,length=400),mean=y,sd=sd), shape=s, rate=r))
                }
              ) 
  }
  # remove small constant
  if (!is.na(small_const)) res<-res-small_const
  # filter out unplausible values
  if (!is.na(threshold0)) res[res<threshold0]<-0
  res
}

#+ Gaussian anamorphosis, back-transformation. from Gaussian to Gamma. Variance
inv_gamma_anamorphosis_var<-function(values,
                                     sigma,
                                     shape,
                                     rate,
                                     xstar) {
# x. vector with the Gaussian values to be transformed into Gamma
# sigma. vector. Gaussian standard deviations (NAs if direct back-transform)
# shape. vector. Gamma shape
# rate. vector. Gamma rate
# xstar. Back-transformed Gamma expected values
#------------------------------------------------------------------------------ 
  return( apply( cbind(values,sqrt(sigma),shape,rate,xstar),
                 MAR=1,
                 FUN=function(x){
                   y<-x[1]; sd<-x[2]; s<-x[3]; r<-x[4]; z<-x[5]
                   if (sd==0) return(y)
##                   seq<-seq( (y-3*sd), (y+3*sd), length=100)
##                   dy<-abs(seq[2]-seq[1])
##                   aux<-inv_gamma_anamorphosis(x=seq, shape=s, rate=r)
##                   ix<-which(is.finite(aux))
##                   sum( aux[ix]**2 * dnorm(seq[ix],mean=y,sd=sd) * dy) - z**2
                  seq <- qnorm(seq(0.001,0.999,length=400),mean=y,sd=sd)
                  aux <- inv_gamma_anamorphosis(values=seq, shape=s, rate=r)
                  ix  <- which(is.finite(aux))
                  n<-length(ix)
                  xgamma <- aux[ix]
                  wnorm  <- dnorm(seq[ix],mean=y,sd=sd)
                  xnorm  <- seq[ix]
                  sum( xgamma[2:n]**2 * wnorm[2:n] * diff(xnorm)) - z**2
#                  mean(inv_gamma_anamorphosis(x= (qnorm(seq(0.001,0.999,length=400)**2),mean=y,sd=sd), shape=s, rate=r)) -z**2
                 }
               ) 
        )
}


inv_gamma_anamorphosis_constrOptim <- function( i,
                                                mean_ga=0,
                                                sd_ga=1,
#                                                prob_min=0.0001,
#                                                prob_max=0.9999,
#                                                prob_length=400,
                                                prob_min=0.001,
                                                prob_max=0.999,
                                                prob_length=100,
                                                rr_inf = 0.1,
                                                cv_mode=F ) {
#
#------------------------------------------------------------------------------
  # initialization
  if (cv_mode) {
    mean     <- yav_henoi_mean[i]
    sd       <- sqrt( yav_henoi_var[i])
    shape_ga <- shape[i]
    rate_ga  <- rate[i]
  } else {
    mean     <- xa_henoi_mean[i]
    sd       <- sqrt( xa_henoi_var[i])
    shape_ga <- shape[i]
    rate_ga  <- rate[i]
  }
  # inverse transformation
  # the scalar one
  if ( is.na( sd)) {
    xpv <- qgamma( pnorm( mean, mean=mean_ga, sd=sd_ga), 
                   shape=shape_ga, rate=rate_ga)
    var <- NA
    shape_final <- NA
    rate_final  <- NA
  } else {
    # the pdf one
    #  transform from normal to gamma
    qseq0   <- seq( prob_min, prob_max, length=prob_length)
    qnorm0  <- qnorm( qseq0, mean = mean, sd = sd)
    qgamma0 <- qgamma( pnorm( qnorm0, mean=mean_ga, sd=sd_ga), 
                       shape=shape_ga, rate=rate_ga)
    ix <- which( is.finite( qgamma0))
    qgamma <- qgamma0[ix]
    qseq   <- qseq0[ix]
    #  first guesses
    xpv <- mean( qgamma)
    var <- mean( qgamma**2) - xpv**2
    sdev <- ifelse( var<0, 0, sqrt(var) )
    #
    if ( xpv < rr_inf & sdev < (rr_inf/2) ) {
#    if ( xpv < rr_inf & sdev < rr_inf ) {
      xpv <- 0
      var <- NA
      shape_final <- NA
      rate_final  <- NA
    } else if ( sdev < (rr_inf/2) ) {
#    } else if ( sdev < rr_inf ) {
      xpv <- xpv 
      var <- NA
      shape_final <- NA
      rate_final  <- NA
    } else {
      shape_start <- xpv**2 / var
      rate_start  <- xpv    / var
      res <- try( constrOptim( theta = c( shape_start, rate_start), 
              f = function(x,q,qseq) {
                qgammafun <- qgamma( qseq, shape=x[1], rate=x[2])
                ix <- which( is.finite(qgammafun) )
                log( sum( ( qgammafun[ix] - q[ix])**2)) 
                },
                               grad  = NULL, 
                               ui    = rbind( c(1,0),    # shape > 0
                                              c(0,1) ),  # rate  > 0
                               ci    = c(0,0),            # the thresholds
                               q     = qgamma,
                               qseq  = qseq ), silent=T)
      if ( class(res) == "try-error") {
        shape_final <- shape_start 
        rate_final  <- rate_start
      } else {
        shape_final <- res$par[1]
        rate_final  <- res$par[2]
        xpv <- shape_final / rate_final
        var <- shape_final / rate_final**2
      }
    }
  }
  return( c( shape_final, rate_final, xpv, var)) 
}
