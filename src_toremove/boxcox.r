#
boxcox<-function(x,lambda) {
  if (lambda==0) {
    return(log(x))
  } else {
    return((x**lambda-1)/lambda)
  }
}


#+ Box-Cox back-transformation function for precipitation values
tboxcox<-function(x,
                  x_mean_sd_in_xbctmp=F,
                  lambda=0.5,
                  threshold=NA,
                  distribution=F,
                  statistics=NULL,
                  method=NULL,
                  p=NULL) {
#------------------------------------------------------------------------------
# x. numeric. It can be either an array of values to be back-transformed or
#             a distribution of Box-Cox transformed values to be used to 
#             estimate the mean (expected) value of the back-transformed
#             non-Gaussian distribution
# lambda. numeric. Box-Cox transformation parameter
# threshold. numeric. Box-Cox transformed rain/no-rain threshold (used only if
#                     distribution is TRUE)
# distribution. logical. TRUE=x values define an empirical Gaussian distribution
#                        FALSE=x values are an array of (independent) numbers
# statistics character. distribution statistics to be computed
#                       "mean_sd" mean and (pseuso) standard-deviation
#                       "quantile" quantile corresponding to probability p
# method. character. Approach to compute the back-transformed mean for a 
#                    distribution or values: "quantile" or "Taylor"
# p. numeric. probability (0-1) for quantile statistics
#==============================================================================
  if (is.null(x)) return(NULL)
  # DITRIBUTION=F
  if (!distribution) {
    if (lambda==0) {
      return(exp(x))
    } else {
      return((1+lambda*x)**(1./lambda))
    }
  # DITRIBUTION=T
  } else {
    if (is.null(statistics)) {
      print("ERROR in tboxcox: distribution=T and statistics is NULL")
      return(NULL)
    }
    if (!(statistics %in% c("mean_sd","quantile"))) {
      print("ERROR in tboxcox: distribution=T and statistics is not defined")
      return(NULL)
    }
    #
    if (!x_mean_sd_in_xbctmp) {
      if (length(x)<=1) return(c(NA,NA))
      mean<-mean(x,na.rm=T)
      if (is.na(mean)) return(c(NA,NA))
      if (!is.na(threshold)) {
       n0<-length(which(x<threshold))
       if (n0>(0.5*length(x))) return(c(0,0))
      }
      sd<-sd(x,na.rm=T)
    } else {
      mean<-xbctmp[x,1]
      sd<-xbctmp[x,2]
    }
    if (is.na(mean)) return(c(NA,NA))
    if (mean<threshold) return(c(0,NA))
    if (is.na(sd)) return(c(tboxcox(mean,lambda=lambda,distribution=F),NA))
    if (sd==0) return(c(tboxcox(mean,lambda=lambda,distribution=F),0))
    #
    if (statistics=="mean_sd") {
      if (is.null(method)) {
        print("ERROR in tboxcox: distribution=T and statistics=\"mean_sd\" and method is NULL")
        return(NULL)
      }
      if (!(method %in% c("quantile","Taylor"))) return(NULL)
      if (method=="quantile") {
        qn<-qnorm(mean=mean,sd=sd,p=seq(1/400,(1-1/400),length=399))
        if (lambda!=0) qn[qn<(-1/lambda)]<-(-1/lambda)
        tqn<-tboxcox(qn,lambda=lambda,distribution=F)
        tmean<-mean(tqn)
        # robust estimator of the (pseudo) standard deviation for a non-Gaussian distribution (Lanzante)
        tsd<-diff(as.vector(quantile(tqn,probs=c(0.25,0.75))))/1.349
      } else if (method=="Taylor") {
        if (lambda==0) {
          f<- exp(mean)
          f1<-exp(mean)
          f2<-exp(mean)
        } else {
          f<- (lambda*mean+1)**(1/lambda)
          f1<-(lambda*mean+1)**(1/lambda-1) # first derivative
          f2<-(1-lambda)*(lambda*mean+1)**(1/lambda-2) # second der
        }
        # mean and standard deviation in back-transformed space (2nd order approx)
        tmean<-f+0.5*sd*sd*f2
        tsd<-sqrt(f1**2 * sd**2 + 0.5 * f2**2 * sd**4)
      } # end if over methods
      return(c(tmean,tsd))
    } else if (statistics=="quantile") {
      if (is.null(p)) {
        print("ERROR in tboxcox: distribution=T and statistics=\"quantile\" and p is NULL")
        return(NULL)
      }
      return(tboxcox(qnorm(mean=mean,sd=sd,p=p),lambda=lambda,distribution=F))
    } # end if statistics
  } # end if over array/distribution
} # end function

