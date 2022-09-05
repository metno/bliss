#alignment techniques
#+ compute optical flow in observed space using hierarchical HS algorithm
optical_flow_HS <- function(r1, r2, nlevel,
                            niter=100, w1=100, w2=0) {
# Adapted From Yue (Michael) Ying repository on github
# https://github.com/myying/QG_Multiscale_DA/blob/master/demo_optical_flow.py
# Main modifications are for making this work in R and for using raster
# GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
# r1 and r2 are rasters
#=============================================================================
  u <- v <- r1
  u[] <- v[] <- 0
  ni <- nrow(r1)
  nj <- ncol(r1)
  res_x <- res(r1)[1]
  res_y <- res(r1)[2]
  for (lev in nlevel:0) {
    cat(paste(" ** lev=",lev," "))
#    r1warp <- warp( r1, -u, -v, method="simple")
    r1warp <- warp( r1, -u, -v, method="bilinear")
    if (any(is.na(getValues(r1warp)))) r1warp[is.na(r1warp)] <- 0
    r1c <- coarsen( r1warp, lev)
    r2c <- coarsen(r2, lev)
    Im1c <- as.matrix(r1c)
    Im2c <- as.matrix(r2c)

    Ix <- 0.5*(deriv_x(Im1c) + deriv_x(Im2c))
    Iy <- 0.5*(deriv_y(Im1c) + deriv_y(Im2c))
    It <- Im2c - Im1c
    du <- dv <- Ix
    du[] <- dv[] <- 0
    for (k in 1:niter) {
      if ((k %% 20)==0) cat(".")
      ubar2 <- laplacian(du) + du
      vbar2 <- laplacian(dv) + dv
      ubar1 <- deriv_xx(du) + du
      vbar1 <- deriv_yy(dv) + dv
      uxy <- deriv_xy(du)
      vxy <- deriv_xy(dv)
      du <- (w1*ubar2 + w2*(ubar1+vxy))/(w1+w2) - Ix*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
      dv <- (w1*vbar2 + w2*(vbar1+uxy))/(w1+w2) - Iy*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
    }
    rdu <- rdv <- r1c
    rdu[] <- du * 2**lev * res_x
    rdv[] <- dv * 2**lev * res_y
    u <- u + sharpen(rdu, r1, method="bilinear")
    v <- v + sharpen(rdv, r1, method="bilinear")
#    u <- u + sharpen(rdu, r1, method="ngb")
    if (any(is.na(getValues(u)))) u[is.na(u)] <- 0
#    v <- v + sharpen(rdv, r1, method="ngb")
    if (any(is.na(getValues(v)))) v[is.na(v)] <- 0
    cat( paste( "u", round(range(getValues(u))[1]/res_x),
                     round(range(getValues(u))[2]/res_x),
                "v", round(range(getValues(v))[1]/res_y),
                     round(range(getValues(v))[2]/res_y)))
#    print(u)
#    print(v)
  }
  cat("\n")
  return( list( u=u, v=v))
}

#source("~/projects/bliss/src/optflow_msa.r")
#source("~/projects/bliss/src/optflow_util.r")
#of <- optical_flow_HS(rb,ra,nlevel=8)
#rbmod<-warp(rb,-of$u,-of$v)
#if (any(is.na(getValues(rbmod)))) rbmod[is.na(rbmod)] <- 0
#image(rbmod,breaks=c(0,seq(0.1,45,length=10)),col=c("gray",rev(rainbow(9))))
#image(ra-rb,breaks=seq(-31,31,length=10),col=c(rev(rainbow(9))))
#image(ra-rbmod,breaks=seq(-31,31,length=10),col=c( rev(rainbow(9))[1:4],"gray",rev(rainbow(9))[6:9]))

#r<-env$rmaster
#r[]<-env$Xa[,1]
#mat<-as.matrix(r)
#dy<-deriv_y(mat)
#u[] <- array(data=dy, dim=c(env$ny,env$nx))
#
#  ni, nj = Im1.shape
#  u = np.zeros((ni, nj))
#  v = np.zeros((ni, nj))
#  for lev in range(nlevel, -1, -1):
#    Im1warp = util.warp(Im1, -u, -v)
#    Im1c = util.coarsen(Im1warp, lev)
#    Im2c = util.coarsen(Im2, lev)
#    
#    niter = 100
#    w1 = 100
#    w2 = 0
#    Ix = 0.5*(util.deriv_x(Im1c) + util.deriv_x(Im2c))
#    Iy = 0.5*(util.deriv_y(Im1c) + util.deriv_y(Im2c))
#    It = Im2c - Im1c
#    du = np.zeros(Ix.shape)
#    dv = np.zeros(Ix.shape)
#    for k in range(niter):
#      ubar2 = util.laplacian(du) + du
#      vbar2 = util.laplacian(dv) + dv
#      ubar1 = util.deriv_xx(du) + du
#      vbar1 = util.deriv_yy(dv) + dv
#      uxy = util.deriv_xy(du)
#      vxy = util.deriv_xy(dv)
#      du = (w1*ubar2 + w2*(ubar1+vxy))/(w1+w2) - Ix*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
#      dv = (w1*vbar2 + w2*(vbar1+uxy))/(w1+w2) - Iy*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
#
#    u += util.sharpen(du*2**lev, lev)
#    v += util.sharpen(dv*2**lev, lev)
#  return u, 
