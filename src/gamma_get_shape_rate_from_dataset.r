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
  #
# Wilks p. 127
#  n<-length(x)
#  x[which(x==0)]<-0.001
#  sum_log_x<-sum(log(x))
#  sum_x<-sum(x)
#  shape_i<-mean**2/var
#  beta_i<-var/mean
#  mat<-matrix(data=NA,dim=c(2,2))
#  vec<-vector(length=2,mode="numeric"); vec[]<-NA
#  for (i in 1:1000) {
#    mat[1,1]<-(-n)*
#    mat[1,2]<-(-n)*beta_i
#    mat[2,1]<-mat[1,2]
#    mat[2,2]<-n*alpha_i/beta_i**2-2*sum(x)/beta_i**3
#    vec[1]<-
#    shape_i<-shape_i-
#  }
}
