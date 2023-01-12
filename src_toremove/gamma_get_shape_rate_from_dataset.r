#+ estimation of the gamma-pdf parameters
gamma_get_shape_rate_from_dataset<-function(x,by=0.01) {
  df <- approxfun(density(x))
  seq <-seq(min(x,na.rm=T),max(x,na.rm=T),by=by)
  mean<-0
  var<-0
  for (xx in seq) {
    if (!is.na(df(xx))) {
      mean<- mean + xx    * df(xx) * by
      var <- var  + xx**2 * df(xx) * by
    }
  }
  var<-var-mean**2
  return(list(shape=(mean**2/var),rate=(mean/var)))
}

#+ estimation of the gamma-pdf parameters
gamma_get_shape_rate_from_dataset_constrOptim<-function( x, 
                                                         small_const=10**(-10)) {
# Wilks(2019) p. 127. It formulates MLE for shape and scale (=1/rate)
#------------------------------------------------------------------------------
  x[which(x==0)]<-small_const
  par_start <- c( mean(x)**2/var(x), var(x)/mean(x))
  gamma_est_mle<-function(par,x) {
    n <- length(x)
    mat <- array( data=NA, dim=c(2,2))
#    mat[1,1] <- -n*( gamma(par[1])*(trigamma(par[1])+digamma(par[1])**2))
    mat[1,1] <- -n*trigamma(par[1])
    mat[1,2] <- -n/par[2]
    mat[2,1] <- mat[1,2]
    mat[2,2] <- n*par[1]/par[2]**2 - 2*sum(x)/par[2]**3
    vec <- vector(mode="numeric", length=2)
#    vec[1] <- sum(log(x)) - n*log(par[2]) - n*gamma(par[1])*digamma(par[1])
    vec[1] <- sum(log(x)) - n*log(par[2]) - n*digamma(par[1])
    vec[2] <- sum(x)/par[2]**2 - n*par[1]/par[2]
    par_new <- par - solve(mat) %*% vec
    mean(abs(par-par_new)/abs(par))
  }
  res <- constrOptim( theta = par_start, 
                      f     = gamma_est_mle,
                      grad  = NULL, 
                      ui    = rbind( c(1,0),    # shape > 0
                                     c(0,1) ),  # scale > 0
                      ci    = c(0,0),           # the thresholds
                      x     = x)
  return( list( shape=res$par[1], rate=1/res$par[2]))
}
