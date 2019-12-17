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
inv_gamma_anamorphosis<-function(x,
                                 sigma=NA,
                                 shape,
                                 rate,
                                 small_const=NA,
                                 threshold0=0,
                                 gauss_mean=0,
                                 gauss_sd=1) {
# x. vector. Gaussian values to be transformed into Gamma
# sigma. vector. Gaussian standard deviations (NAs if direct back-transform)
# shape. vector. Gamma shape
# rate. vector. Gamma rate
# small_const. small constant to be removed to the Gamma values
# gauss_mean. Gaussian's mean value
# gauss_sd. Gaussian's standard deviation
#------------------------------------------------------------------------------ 
  # direct back-transformation
  if (any(is.na(sigma))) {
    res<-qgamma( pnorm(x,mean=gauss_mean,sd=gauss_sd), shape=shape, rate=rate)
  # PDF-based back-transformation
  } else {
    res<-apply( cbind( x, sqrt(sigma), shape, rate ),
                MAR=1,
                FUN=function(x){
                  y<-x[1]; sd<-x[2]; s<-x[3]; r<-x[4]
                  if (sd==0) return(y)
                  seq<-seq( (y-3*sd), (y+3*sd), length=100)
                  dy<-abs(seq[2]-seq[1])
                  aux<-inv_gamma_anamorphosis(x=seq, shape=s, rate=r)
                  ix<-which(is.finite(aux))
                  sum( aux[ix] * dnorm(seq[ix],mean=y,sd=sd) * dy) 
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
inv_gamma_anamorphosis_var<-function(x,
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
  return( apply( cbind(x,sqrt(sigma),shape,rate,xstar),
                 MAR=1,
                 FUN=function(x){
                   y<-x[1]; sd<-x[2]; s<-x[3]; r<-x[4]; z<-x[5]
                   if (sd==0) return(y)
                   seq<-seq( (y-3*sd), (y+3*sd), length=100)
                   dy<-abs(seq[2]-seq[1])
                   aux<-inv_gamma_anamorphosis(x=seq, shape=s, rate=r)
                   ix<-which(is.finite(aux))
                   sum( aux[ix]**2 * dnorm(seq[ix],mean=y,sd=sd) * dy) - z**2
                 }
               ) 
        )
}

