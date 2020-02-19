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
                                                prob_min=0.001,
                                                prob_max=0.999,
                                                prob_length=100,
                                                cv_mode=F ) {
#
#------------------------------------------------------------------------------
  if (cv_mode) {
    mean     <- yav_henoi_mean[i]
    sd       <- sqrt(yav_henoi_var[i])
    shape_ga <- shape[i]
    rate_ga  <- rate[i]
  } else {
    mean     <- xa_henoi_mean[i]
    sd       <- sqrt(xa_henoi_var[i])
    shape_ga <- shape[i]
    rate_ga  <- rate[i]
  }
  #
  qseq   <- seq( prob_min, prob_max, length=prob_length)
  qnorm  <- qnorm( qseq, mean = mean, sd = sd)
  qgamma <- qgamma( pnorm( qnorm, mean=mean_ga, sd=sd_ga), 
                    shape=shape_ga, rate=rate_ga)
  xpv <- mean( qgamma)
  var <- mean( qgamma**2) - xpv**2
  shape_start <- xpv**2 / var
  rate_start  <- xpv    / var
  #
  res <- constrOptim( theta = c( shape_start, rate_start), 
                      f     = function(x,q,qseq)
                       { log( sum( ( qgamma( qseq, shape=x[1], rate=x[2]) - q)**2)) },
                      grad  = NULL, 
                      ui    = rbind( c(1,0),    # shape > 0
                                     c(0,1) ),  # rate  > 0
                      ci    = c(0,0),            # the thresholds
                      q     = qgamma,
                      qseq  = qseq ) 
  return(c(res$par[1],res$par[2],mean,sd,shape_ga,rate_ga))
}
