#+ observation operator, precipitation
obsop_precip<-function(method="bilinear",inf=0) {
# method. "simple" nearest neighbour. "bilinear" bilinear
# inf. lower limit for the plausible values
# - Global variables
# rmaster. raster with grid parameters
# n0. number of observations
# VecX, VecY. observation coordinates
# nens. number of ensemble members
# aix. indexes. the ensemble members (xb) are defined only over aix cells
# xb. matrix. ensemble members (aix cells, nens ensembles)
# - Returned values
# yb. matrix. ensemble at observation locations (ix observation loc, nens)
# ix. vector. indeces over the observation locations.
#------------------------------------------------------------------------------
  r<-rmaster
  yb0<-array(data=NA,dim=c(n0,nens))
  for (e in 1:nens) {
    r[]<-NA
    r[aix]<-xb[,e]
    yb0[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
    auxx<-which(yb0[,e]<0 | is.na(yb0[,e]) | is.nan(yb0[,e]))
    if (length(auxx)>0) yb0[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
    rm(auxx)
  }
  rm(r)
  # index over yb with all the ensemble members finite and not NAs
  ix<-which(apply(yb0,MAR=1,
                  FUN=function(x){length(which(!is.na(x) & !is.nan(x)))})
            ==nens)
  n1<-length(ix)
  yb<-array(data=NA,dim=c(n1,nens))
  for (e in 1:nens) yb[,e]<-yb0[ix,e]
  rm(yb0)
  yb[yb<inf]<-inf
  return(list(yb=yb,ix=ix))
}

