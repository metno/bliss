#+ Define correlations
corr2d <- function( values, par, label) {
# correlations between pair of points
#------------------------------------------------------------------------------
  if ( length( par == 1)) {
    x  <- values[,1]
    y  <- values[,2]
    dh <- par
    if (corr=="gaussian") {
      S<-exp(-0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
    } else if (corr=="soar")  {
      distnorm<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2) / dh 
      S<-(1+distnorm)*exp(-distnorm)
    } else if (corr=="powerlaw")  {
      S<-1 / (1 + 0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
    } else if (corr=="toar")  {
      dist<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)
      S<- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
    }
  }
  res
}
